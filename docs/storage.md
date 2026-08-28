# Storage & index

The store has two components under one directory: a sparse Zarr store
holding pixel data, and a SQLite database recording what that data is
and how it was obtained.

```
{config.tmp_dir}/sentinel2_cube/
├── index.db          # the ledger (SQLite, WAL mode)
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

## The SQLite ledger

Three tables, none holding pixel data (`pysentinel2/index.py`):

```mermaid
erDiagram
    scenes {
        TEXT item_id PK "STAC item id"
        TEXT solar_day "YYYY-MM-DD"
        REAL cloud_cover "eo:cloud_cover, nullable"
        TEXT item_json "full STAC item"
    }
    chunks {
        TEXT solar_day PK
        INTEGER cy PK "chunk row"
        INTEGER cx PK "chunk column"
        TEXT written_at
    }
    searches {
        REAL x0 "searched extent, EPSG:6933"
        REAL y0 "searched extent"
        REAL x1 "searched extent"
        REAL y1 "searched extent"
        TEXT start "date range searched"
        TEXT end "date range searched"
        TEXT searched_at
    }
```

**`scenes`** caches every STAC item ever returned, including its full
JSON. Day selection, cloud filtering and re-fills therefore work
offline; the STAC API is only contacted for regions/ranges never seen
before. Because *all* scenes are recorded regardless of cloud cover,
changing `max_cloud_cover` later is a pure read-side change.

**`coverage`** is the deduplication ledger, at pixel exactness: a row
`(solar_day, row0, row1, col0, col1)` records one written window. A day's
missing work is its requested window minus the union of its rects
(`grid.rect_subtract`); fills append one rect per day they complete.
Legacy stores that used the fixed-chunk `chunks` table migrate
automatically on first open (each chunk becomes a rect; x-adjacent runs
coalesce). The old semantics for reference: a row `(solar_day, cy, cx)`
asserts that this cell of the global grid is completely written for
that day. `fill()` computes `wanted − done` per day and downloads only
the difference.

**`searches`** records the spatial extent and date range of every STAC
search. `search_covered()` answers "is this request contained in some
past search?" — which is what distinguishes *"searched, no scenes
exist"* (a valid, cacheable answer) from *"never asked"*.

## Crash-safety semantics

The ordering of operations ensures that an interrupted fill cannot
corrupt the store:

```mermaid
sequenceDiagram
    participant F as fill()
    participant S3 as DEA S3
    participant Z as cube.zarr
    participant DB as index.db (WAL)

    F->>S3: fetch missing cells (odc.stac.load)
    S3-->>F: pixel window
    F->>Z: write whole chunks
    Note over Z: a crash here leaves chunks written<br/>but unrecorded; they are re-fetched later
    F->>DB: mark_chunks(day, cells) — one transaction
    Note over DB: only now is the cell "done"
```

- Pixels are written before the ledger row, and the ledger commit is a
  single WAL transaction. A crash at any point leaves cells unmarked,
  and the next run re-downloads and overwrites them; because writes are
  whole chunks at fixed grid positions, re-writing a cell reproduces
  the same bytes rather than accumulating inconsistency.
- The inverse failure — a ledger row without its pixels — cannot occur,
  because `mark_chunks` runs only after the Zarr writes return.
- WAL mode allows concurrent readers while a fill is in progress.

## Where the store lives

Locations derive from the shared lab `Config`
(`pysentinel2/paths.py`). The store is keyed by data root, not by
query: every query on the machine reads and fills the same cube.

```python
from pysentinel2.paths import Paths
paths = Paths()          # from the default Config
paths.store              # .../sentinel2_cube/cube.zarr
paths.index_db           # .../sentinel2_cube/index.db
```

The directory is a cache: deleting it loses no state that cannot be
rebuilt, and subsequent queries re-download exactly what they need.
