# pysentinel2

A **local Sentinel-2 datacube that fills itself on demand**. Every pixel
this machine ever downloads lands in one sparse, chunk-indexed store —
so nothing is ever downloaded twice: overlapping areas, extended date
ranges and repeat runs all reuse the same chunks. Part of the
[Borevitz Lab](https://borevitzlab.anu.edu.au/) ecosystem; the default
source is [Digital Earth Australia](https://explorer.dea.ga.gov.au/)'s
ARD collections (`ga_s2am_ard_3` / `ga_s2bm_ard_3`) via STAC.

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
  `Cube.get(query)` diffs the requested (day × chunk) cells against the
  index and downloads **only the missing cells**.
- STAC results are cached (full item JSON) in the index, so re-reads and
  re-fills of known regions work without re-searching. Cloud-cover
  filtering happens at read time from the index — relaxing the threshold
  later needs no re-search.
- Only raw bands (incl. fmask) are stored. `get(clean=True)` applies
  cloud masking **on read** — there is no second "clean" copy on disk,
  roughly halving storage versus a raw+clean layout.
- Writes are whole-chunk and the index is transactional (SQLite/WAL): a
  crash mid-fill just leaves cells unmarked, and the next run resumes.

## Usage

```python
from datetime import date
from borevitz_lab.query import Query
from pysentinel2.cube import Cube

query = Query(
    bbox=[148.36265, -33.52606, 148.38265, -33.50606],  # [W, S, E, N]
    start=date(2024, 1, 1),
    end=date(2024, 12, 31),
    stub='my_aoi',
)

cube = Cube()
ds_raw = cube.get(query)              # raw window incl. fmask, dims (time, y, x)
ds     = cube.get(query, clean=True)  # cloud-masked view of the same pixels

cube.fill(query)                      # → 0: everything already local
```

`download_sentinel2(query)` and `clean_sentinel2(query)` remain as thin
wrappers over `Cube.get` for pipeline compatibility.

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

## Install

```bash
pip install -e ../borevitz_lab   # shared core (not yet on PyPI)
pip install -e .
```

Needs the geospatial native stack (GDAL/PROJ) available for
`rasterio`/`rioxarray`; under conda use the lab's environment.

## Robustness notes

Hardening for DEA's public S3 + STAC quirks is built in (see
`diagnostics.md`): STAC retries with backoff on cold-start 504s, GDAL
low-speed timeouts so stalled reads abort instead of hanging, and
`fail_on_error=False` so one corrupt tile costs a nodata gap rather
than a whole day's fill.

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
