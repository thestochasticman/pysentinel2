# pysentinel2

A **local Sentinel-2 datacube that fills itself on demand**. Every pixel
this machine ever downloads lands in one sparse, chunk-indexed store —
so nothing is ever downloaded twice: overlapping areas, extended date
ranges and repeat runs all reuse the same chunks. Part of the
[Borevitz Lab](https://borevitzlab.anu.edu.au/) ecosystem; the default
source is [Digital Earth Australia](https://explorer.dea.ga.gov.au/)'s
ARD collections (`ga_s2am_ard_3` / `ga_s2bm_ard_3`) via STAC.

Full documentation — architecture, grid geometry, storage, cleaning,
indices, robustness — is in [`docs/`](docs/README.md), with flowcharts
and figures generated from a real store.

![Every stored solar day for the example window](docs/images/cube_frames_rgb.png)
*Contents of the store for a 2 × 2 km example window. Clear, cloudy and
off-swath days are all stored raw and classified at read time.*

## How it works

```
{data_root}/sentinel2_cube/
├── index.db      # SQLite: populated chunks · seen scenes · past searches
└── cube.zarr/
    ├── 2024-01-03/   # one group per solar day
    │   ├── nbart_red # arrays on a fixed EPSG:6933 10 m global grid
    │   └── ...       # sparse: only written 256×256-px chunks exist on disk
    └── 2024-01-08/ ...
```

- Any bbox maps deterministically to a set of ~2.56 km grid chunks.
  `Cube.get_ds(bbox, start, end)` diffs the requested (day × chunk) cells against the
  index and downloads **only the missing cells**.
- STAC results are cached (full item JSON) in the index, so re-reads and
  re-fills of known regions work without re-searching. Cloud-cover
  filtering happens at read time from the index — relaxing the threshold
  later needs no re-search.
- Only raw bands (incl. fmask) are stored. `get_ds(..., clean=True)` applies
  cloud masking **on read** — there is no second "clean" copy on disk,
  roughly halving storage versus a raw+clean layout. See
  [Cleaning & masking](#cleaning--masking) for exactly what the mask does.
- Spectral indices — NDVI, CFI, NIRv, NDTI, CAI — are on-read
  derivatives too: `get_ds(..., indices=('NDVI', 'NIRv'))` computes them
  from cloud-masked reflectance and stores nothing.
- Writes are whole-chunk and the index is transactional (SQLite/WAL): a
  crash mid-fill just leaves cells unmarked, and the next run resumes.

## Usage

The core API is **query-agnostic** — just a bbox and dates, no setup:

```python
from datetime import date
from pysentinel2.cube import Cube

cube = Cube()
bbox = [148.36265, -33.52606, 148.38265, -33.50606]  # [W, S, E, N]

ds_raw = cube.get_ds(bbox, date(2024, 1, 1), date(2024, 12, 31))
ds     = cube.get_ds(bbox, date(2024, 1, 1), date(2024, 12, 31), clean=True)
ds     = cube.get_ds(bbox, date(2024, 1, 1), date(2024, 12, 31),
                     indices=('NDVI', 'CFI', 'NIRv', 'NDTI', 'CAI'))

cube.fill(bbox, date(2024, 1, 1), date(2024, 12, 31))  # → 0: already local
```

Pipelines that speak the shared `borevitz_lab.query.Query` (the
reproducibility layer — stubs, registry) use the adapters:

```python
ds = cube.get_ds_query(query)            # = cube.get_ds(query.bbox, query.start, query.end)
```

`download_sentinel2(query)` and `clean_sentinel2(query)` remain as thin
wrappers over `Cube.get_ds_query` for pipeline compatibility.

Package design (shared across the lab's packages — no inheritance,
composition only):

- **`Query`** (from `borevitz-lab`) — identity: what region, what dates.
- **`Sentinel2`** (`pysentinel2.sentinel2`) — config: STAC URL,
  collections, bands, CRS, cloud threshold, fmask codes.
- **`Paths`** (`pysentinel2.paths`) — derived locations of the store for
  a given `Config`.
- **`grid`** — the fixed global grid (pure, offline-testable math).
- **`Index`** (`pysentinel2.index`) — the SQLite ledger.
- **`Cube`** (`pysentinel2.cube`) — ties them together.

## Cleaning & masking

`clean=True` (and any `indices=` request, which implies it) runs the
window through `pysentinel2.cube.clean_dataset`. The design principle:
**invalid and contaminated are different things.**

| Pixel state | fmask | Meaning | Treatment |
|---|---|---|---|
| Invalid | 0 (nodata) | Outside the scene footprint / never sensed | → NaN; counts *against coverage*, not against cloudiness |
| Clear | 1 | Usable land observation | kept |
| Cloud | 2 | Contaminated | → NaN (dilated) |
| Shadow | 3 | Contaminated | → NaN (dilated) |
| Snow | 4 | Surface state; corrupts vegetation statistics | → NaN by default (`mask_snow=False` to keep); never counts toward the frame gate |
| Water | 5 | Legitimate signal (NDWI, dams, rivers) | kept by default (`mask_water=True` to drop) |

Contaminated pixels are dilated before masking, frames are gated on the
two fractions independently, and every read is annotated with the
statistics and thresholds that produced it. Nothing is persisted —
different thresholds on the same window are just different reads of the
same raw store. Full pipeline, tunables and figures:
[docs/cleaning.md](docs/cleaning.md).

## Performance

Live measurements against DEA — a ~2 × 2 km AOI, 11-band ARD at 10 m
(one *cell* = one 256 × 256-px chunk on one solar day):

| Scenario | Downloaded | Time |
|---|---|---|
| Cold fill — 3 weeks (3 clear scenes) | 12 cells | 5.7 s |
| Same request again | nothing | **0.0 s** |
| AOI shifted 1 km (inside cached chunks) | nothing | **0.0 s** |
| Date range extended +1 month | 32 cells — *new days only* | 17.2 s |
| Read cached window (512² px × 3 days × 11 bands) | — | 0.13 s |
| Read cached window, cloud-masked (`clean=True`) | — | 0.23 s |

Store footprint: **13.6 MB for 11 solar days** — raw + fmask only, since
the clean cube is a 0.1 s on-read transform rather than a second copy.

Absolute times vary with network and DEA load. The zero rows are the
significant ones: those requests are resolved by index lookups alone,
with no network access.

Multi-year fills are batched, not per-day: all 11 bands for up to 64
missing days come down in one bulk load per batch (see
[the fill algorithm](docs/architecture.md#the-fill-algorithm)), keeping
the I/O threads saturated across day boundaries — a two-month cold fill
measured 20-26 s where a per-day loop measured 33 s on a healthy DEA
and 270 s on a degraded one, and the gap widens with the length of the
range. An earlier fmask-first screening pass was removed after
measurement: it skipped 8.8% of days' reflectance while paying an extra
request round on every day.

## Install

### Conda (recommended)

```bash
conda install -c conda-forge -c thestochasticman pysentinel2
```

### From source

All lab repos share one conda environment, `borevitz_lab` — each repo's
`environment.yml` creates it if missing and adds its own packages if it
exists (never use `--prune`):

```bash
conda env update -n borevitz_lab -f environment.yml
conda activate borevitz_lab
pip install -e ../borevitz_lab   # shared core (not yet on PyPI)
pip install -e .
```

Without conda, you need the geospatial native stack (GDAL/PROJ)
available system-wide for `rasterio`/`rioxarray`, then the same two
`pip install -e` lines.

## Robustness notes

Hardening for DEA's public S3 + STAC quirks (cold-start 504s, stalled
reads, corrupt tiles) is built in — see
[docs/robustness.md](docs/robustness.md) and, for the underlying
investigations, [`diagnostics.md`](diagnostics.md).

## Test

```bash
# offline (pure math + synthetic store):
python pysentinel2/grid.py    # True
python pysentinel2/index.py   # True
python pysentinel2/paths.py   # True
python pysentinel2/cube.py    # True

# live (small real downloads from DEA, incl. dedup assertions):
python pysentinel2/download_sentinel2.py  # True
python pysentinel2/clean_sentinel2.py     # True
```
