# Storage & index

The store has two components under one directory: a sparse Zarr store
holding pixel data, and a tree of JSON markers recording what that data is
and how it was obtained.

```
{config.tmp_dir}/sentinel2_cube/
├── index/            # the ledger: coverage rects, scenes, searches (JSON markers)
└── cube.zarr/
    ├── 2023-12-18/   # one group per solar day
    │   ├── nbart_red         # one array per band on the full global grid
    │   ├── nbart_green
    │   ├── ...
    │   └── oa_fmask
    └── 2024-01-22/ ...
```

## The Zarr store

Each solar-day group holds one array per band, logically global
(3 473 920 × 1 465 344 px on the [fixed grid](grid.md)) but physically
sparse: Zarr materialises only chunks that have been written. Band
arrays are chunked at 256 × 256 px, so the on-disk chunk coincides with
the grid's unit of deduplication.

- Reflectance bands are `int16` with nodata −999 (DEA ARD convention);
  `oa_fmask` is `uint8` with nodata 0. The nodata value is recorded as
  an array attribute and doubles as the Zarr fill value, so reading an
  unwritten region yields nodata — the same value as a region the
  satellite never sensed, and the [cleaning pipeline](cleaning.md)
  treats the two identically.
- Because fills are [fmask-first](architecture.md#the-fill-algorithm),
  a day group always holds the fmask array, while reflectance arrays
  exist only where at least one chunk passed the download screen; the
  read path returns nodata for reflectance bands that were screened
  out.
- Grouping by **solar day** (the UTC acquisition time shifted by the
  scene-centre longitude in degrees, $t_{solar} = t_{UTC} + \lambda / 15$ hours)
  merges the two Sentinel-2 satellites and adjacent swath tiles into one
  temporal layer per overpass day.
- Storage cost scales with *observed area × observed days*, not with
  the global grid: the example window (4 chunks × 12 days × 11 bands)
  occupies ≈ 14 MB.

## The marker index

Three marker trees under `index/`, none holding pixel data
(`pysentinel2/index.py`), each file committed by write-to-temp + atomic
rename (`troi.ledger.Markers`):

| tree | one file per | holds |
|---|---|---|
| `scenes/<YYYY>/<YYYY-MM-DD>/<item_id>.json` | STAC item ever seen | solar day, `eo:cloud_cover`, the full item JSON |
| `coverage/<YYYY-MM-DD>/<uuid>.json` | written pixel rect | `row0, row1, col0, col1` on the global grid |
| `searches/<uuid>.json` | STAC search ever run | EPSG:6933 bbox and date range |

- **scenes** lets day selection and re-fills work offline; cloud
  filtering is applied at read time, so a laxer threshold later needs
  no re-search.
- **coverage** is the dedup ledger at pixel exactness: a day's missing
  work is `grid.rect_subtract(tight_window, covered_rects(day))`.
- **searches** distinguishes *"searched, no scenes exist"* (a valid,
  cacheable answer) from *"never asked"*.

Why files and not SQLite: the cube is filled by many PBS jobs on many
Gadi nodes against one store on Lustre, which is mounted `localflock`,
so file locks are node-local and SQLite is unsafe in any journal mode.
What Lustre does serialise across nodes is `mkdir` and `rename`; the
index and the claims are built from those alone
(`troi/docs/ledger.md`).

## Crash-safety semantics

The ordering of operations ensures that an interrupted fill cannot
corrupt the store:

```mermaid
sequenceDiagram
    participant F as fill()
    participant S3 as DEA S3
    participant Z as cube.zarr
    participant C as claims/s2-<day>
    participant IX as index/coverage

    F->>S3: fetch missing rects (odc.stac.load)
    S3-->>F: pixel window
    F->>C: mkdir (one writer per day, across nodes)
    F->>Z: write the rect (zarr commits each chunk by rename)
    Note over Z: a crash here leaves pixels written<br/>but unrecorded; they are re-fetched later
    F->>IX: mark_rect(day, rect) — one file, atomic rename
    Note over IX: only now is the rect "done"
    F->>C: rmdir
```

- Pixels are written before the marker, and the marker is one atomic
  rename. A crash at any point leaves rects unmarked, and the next run
  re-downloads and overwrites them; re-writing a rect reproduces the
  same bytes rather than accumulating inconsistency.
- The inverse failure — a marker without its pixels — cannot occur,
  because `mark_rect` runs only after the Zarr writes return.
- Rect writes are partial chunks (read-modify-write), so a day is
  written under a claim directory; a claim abandoned by a dead job
  expires after its lease and is taken over. Readers never block.

## Where the store lives

Locations derive from the shared lab `Config`
(`pysentinel2/paths.py`). The store is keyed by data root, not by
troi: every troi on the machine reads and fills the same cube.

```python
from pysentinel2.paths import Paths
paths = Paths()          # from the default Config
paths.store              # .../sentinel2_cube/cube.zarr
paths.index              # .../sentinel2_cube/index
```

The directory is a cache: deleting it loses no state that cannot be
rebuilt, and subsequent queries re-download exactly what they need.
