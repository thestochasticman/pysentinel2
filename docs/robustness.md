# Robustness

DEA's public STAC endpoint and S3 bucket are reliable in aggregate but
exhibit reproducible transient failure modes. The hardening below is
built into `pysentinel2/cube.py`; the underlying investigations, with
shell reproductions, are recorded in
[`diagnostics.md`](../diagnostics.md).

## Failure modes and mitigations

| Symptom | Root cause | Mitigation (built in) |
|---|---|---|
| `APIError: 504 Gateway Time-out` on the first STAC request | DEA's STAC load balancer takes ~30 s to serve a cold-cache request, then warms up | `StacApiIO` with `urllib3.Retry`: 5 attempts, exponential backoff (1 s → 8 s), on HTTP 408/429/502/503/504 |
| Download hangs indefinitely mid-fill | DEA S3 intermittently half-closes connections (sockets stuck in `CLOSE_WAIT`); a blocked read inside `rasterio` never times out by default | GDAL low-speed timeout: abort any transfer making < 1 B/s for 60 s (`GDAL_HTTP_LOW_SPEED_*`), then up to 5 GDAL-level retries with delay |
| `RasterioIOError('Unsupported Authorization Type')` from a worker | On `dask` distributed clusters, the unsigned-S3 GDAL configuration broadcast can race the first task, so a worker opens an asset with default auth settings | In-process **threaded** scheduler only — workers share the main process's GDAL/CURL configuration by construction; no broadcast exists to race |
| One corrupt/unreadable tile aborts a whole day's fill | `odc.stac.load` defaults to fail-fast | `fail_on_error=False`: a bad asset costs a nodata gap in that band, not the fill |

The cumulative retry backoff (≈ 15 s) covers DEA's observed ~30 s
warm-up window without hammering the endpoint. Because search results
(full STAC item JSON) are cached in the [index](storage.md), this cost
is paid at most once per new region/date-range; every later request
over covered extents performs no STAC traffic at all.

## GDAL/CURL environment

Set at import time (via `os.environ.setdefault`, so user overrides win)
and mirrored into `odc.stac.configure_rio`:

```python
GDAL_HTTP_CONNECTTIMEOUT = 20    # s to establish a connection
GDAL_HTTP_LOW_SPEED_TIME = 60    # abort after 60 s below…
GDAL_HTTP_LOW_SPEED_LIMIT = 1    # …1 byte/s (i.e. no progress)
GDAL_HTTP_MAX_RETRY = 5
GDAL_HTTP_RETRY_DELAY = 1        # s, GDAL doubles internally
CPL_VSIL_CURL_USE_HEAD = NO      # DEA S3 serves ranged GETs; skip HEAD probes
```

The low-speed settings carry most of the value: they distinguish a
stalled transfer (no progress; abort and retry) from a slow but
progressing one (left alone), so large fills over slow links are not
terminated by an absolute timeout.

## Consequences of failure

The [crash-safety ordering](storage.md#crash-safety-semantics) bounds
the cost of any of these failing terminally mid-fill: cells not yet
marked in the ledger are re-fetched on the next run, and already-marked
cells are never touched again. The mitigations above therefore affect
latency, not correctness.

## Verifying against a live endpoint

```bash
python pysentinel2/download_sentinel2.py   # live: download + dedup assertions
python pysentinel2/clean_sentinel2.py      # live: clean read end-to-end
```

If the cold-cache 504 reappears, confirm with the `curl` loop in
[`diagnostics.md`](../diagnostics.md) before suspecting a regression in
this package.
