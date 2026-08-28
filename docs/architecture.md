# Architecture

## Package composition

The package follows the lab-wide convention of **composition over
inheritance**: each concern lives in one module with a narrow interface,
and `Cube` assembles them. Every module carries its own offline test
suite (`python pysentinel2/<module>.py` → `True`).

```mermaid
flowchart TD
    subgraph borevitz_lab["borevitz-lab (shared core)"]
        CFG["Config<br/><i>where data lives</i>"]
        QRY["Query<br/><i>reproducibility layer:<br/>bbox, dates, stub, registry</i>"]
    end

    subgraph pysentinel2
        S2["Sentinel2<br/><i>STAC URL, collections, bands,<br/>fmask codes, cloud threshold</i>"]
        P["Paths<br/><i>derived store locations</i>"]
        G["grid<br/><i>pure EPSG:6933 chunk math</i>"]
        IX["Index<br/><i>SQLite ledger</i>"]
        CB["Cube<br/><i>fill · read · clean · derive</i>"]
        DRV["derive<br/><i>spectral indices</i>"]
        DL["download_sentinel2 / clean_sentinel2<br/><i>thin Query-compatible wrappers</i>"]
    end

    CFG --> P
    CFG --> CB
    S2 --> CB
    P --> CB
    G --> CB
    IX --> CB
    DRV --> CB
    QRY --> DL --> CB
```

| Component | Module | Responsibility |
|---|---|---|
| `Sentinel2` | `pysentinel2/sentinel2.py` | Immutable configuration: STAC endpoint, the two DEA ARD collections, the 11 stored bands, fmask class codes, per-scene cloud-cover threshold. |
| `Paths` | `pysentinel2/paths.py` | Derives the store locations (`{tmp_dir}/sentinel2_cube/{cube.zarr, index.db}`) from a `Config`. Inputs go in `Config`; derived locations in `Paths`. |
| `grid` | `pysentinel2/grid.py` | The fixed global grid. Pure functions, no I/O — see [The grid](grid.md). |
| `Index` | `pysentinel2/index.py` | The SQLite ledger of populated cells, seen scenes, and past searches — see [Storage & index](storage.md). |
| `Cube` | `pysentinel2/cube.py` | Orchestration: diff → fill → read → (optionally) clean → (optionally) derive. |
| `derive` | `pysentinel2/derive.py` | On-read spectral indices — see [Spectral indices](indices.md). |
| wrappers | `download_sentinel2.py`, `clean_sentinel2.py` | Compatibility entry points for pipelines that speak `borevitz_lab.query.Query`. |

The core API is deliberately query-agnostic — `Cube.get_ds(bbox,
start, end)` needs nothing but a region and a date range. The `Query`
adapters (`get_ds_query`, `fill_query`) exist so pipelines built on the
lab's reproducibility layer plug in without translation code.

## The `fill` algorithm

`fill()` is the write path; `get_ds()` is `fill()` followed by a read.
Accounting is pixel-exact per solar day (coverage rectangles); the unit of
*network work* is a multi-day **batch**: one bulk `odc.stac.load` of
every band (fmask included) per up to 64 missing days, so one dask
graph carries days × chunks × bands tasks and the I/O threads stay
saturated across day boundaries.

```mermaid
flowchart TD
    A["fill(bbox, start, end)"] --> B["window_for_bbox:<br/>snap bbox outward to whole chunks"]
    B --> C{"index.search_covered<br/>(bbox, start, end)?"}
    C -- "no" --> D["STAC search over the window<br/>(urllib3 retries, backoff)"]
    D --> E["upsert full item JSON per scene;<br/>record the search extent"]
    C -- "yes" --> F["index.scenes_for_range:<br/>solar days ≤ max_cloud_cover"]
    E --> F
    F --> G["plan: tight window − coverage rects<br/>per day → missing pixel rects"]
    G --> H["<b>bulk load</b> — all 11 bands,<br/>one odc.stac.load per ≤64-day batch"]
    H --> I{"integrity gate per day:<br/>reflectance mostly nodata where<br/>fmask has valid ground?"}
    I -- "fails (bad HTTP reads)" --> J["skip: day left unwritten and<br/>unmarked — next fill retries"]
    I -- "passes" --> K["write whole chunks into the<br/>day's global Zarr arrays"]
    K --> L["index.mark_chunks per day, transactional,<br/><i>after</i> the pixels are on disk"]
    L --> M["return number of cells downloaded<br/>(0 = fully served from cache)"]
    J --> M
```

Batch sizes are capped by window area (~256 MB of int16 in flight) so
small AOIs get maximal batching while large AOIs degrade gracefully
toward per-day loads.

An earlier design fetched the fmask band first and *screened* chunks —
downloading reflectance only where a chunk's own fmask showed enough
clear ground. Measured on a real store the screen skipped reflectance
on 8.8% of days while paying an extra round of COG opens on 100% of
them; in a phase dominated by request latency rather than bytes, that
is a net loss, so the screen was removed. Cloud handling now lives in
the STAC-level `eo:cloud_cover ≤ 80` prefilter and in read-time
[cleaning](cleaning.md); the fmask band is stored for every day either
way, so nothing downstream changed.

Two properties underpin the correctness of the scheme:

- **Downloads are pinned to the grid.** `odc.stac.load` receives a
  `GeoBox` constructed from the chunk-aligned window
  (`grid.geobox_for_window`), so every downloaded pixel lands at its
  single canonical position in the global array. Writes therefore
  cover whole chunks only, and a ledger row is accurate from the
  moment its transaction commits.
- **Search caching is separate from pixel caching.** The `searches`
  table distinguishes a region and range that was queried and returned
  no scenes from one that was never queried; without it, empty regions
  would be re-searched on every request. Scene records store the full
  STAC item JSON, so re-fills and cloud-threshold changes do not
  contact the STAC API again.

## The read path

Reads never touch the network beyond what `fill()` requires:

1. `scenes_for_range` selects solar days from the index, applying the
   per-scene `eo:cloud_cover ≤ max_cloud_cover` filter **at read time**.
   Relaxing the threshold later requires no re-search — newly eligible
   days simply become visible (and are filled on the next request).
2. Each day's bands are sliced from the day's global Zarr arrays at the
   window's pixel coordinates — an index lookup plus a few chunk reads.
3. `clean=True` applies the [cleaning pipeline](cleaning.md);
   `indices=(...)` additionally computes [spectral indices](indices.md).
   Both are pure array transforms on the returned window.

On a warm store, reading a 512 × 512-pixel × 11-band window across
several days takes on the order of 0.1–0.3 s (see the
[performance table](../README.md#performance) in the main README).

## Concurrency model

Downloads run under Dask's in-process threaded scheduler
(`scheduler='threads'`) rather than a distributed cluster. The choice
is deliberate: the workload is network-bound (S3 range reads), threads
share the GDAL/CURL configuration set in the main process, and the
distributed client exhibited a startup race in which workers opened
assets before receiving the unsigned-S3 configuration (see
[Robustness](robustness.md)).
