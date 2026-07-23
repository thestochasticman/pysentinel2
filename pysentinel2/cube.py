"""One machine-wide Sentinel-2 datacube that fills itself on demand.

Instead of one isolated Zarr per (bbox, time) query, every pixel this
machine ever downloads lands in a single sparse store on the fixed
EPSG:6933/10 m grid (:mod:`pysentinel2.grid`):

    {config.tmp_dir}/sentinel2_cube/
    ├── index.db      # what's populated / what's been searched (pysentinel2.index)
    └── cube.zarr/
        └── 2024-01-03/   # one group per solar day, one array per band
            ├── nbart_red # global-grid array; only written chunks exist on disk
            └── ...

``Cube.get(bbox, start, end)`` diffs the requested (day x chunk) cells
against the index, downloads only the missing cells, then reads the
window. Nothing is ever fetched twice — overlapping bboxes, extended
date ranges and repeat runs all reuse the same chunks. Only the raw
bands (incl. fmask) are stored; ``get(..., clean=True)`` applies cloud
masking on read, so no second "clean" copy exists on disk.

The core API is query-agnostic (bbox + dates — the data layer);
``get_query`` / ``fill_query`` adapt a :class:`borevitz_lab.query.Query`
(the reproducibility layer) onto it.
"""

import os
import dask
import numpy as np
import xarray as xr
import zarr
from attrs import frozen, field
from datetime import date, datetime, timedelta, timezone
from os import makedirs
from xarray import Dataset

from borevitz_lab.config import Config, config as default_config
from pysentinel2 import grid
from pysentinel2.index import Index
from pysentinel2.paths import Paths
from pysentinel2.sentinel2 import Sentinel2, defaultsentinel2

# GDAL/CURL hardening for reads from DEA's public S3 over /vsicurl. Without a
# timeout, a stalled connection (DEA's S3 intermittently half-closes sockets —
# they show up in CLOSE_WAIT) leaves the rasterio.warp.reproject read blocked
# forever, hanging the whole download. The low-speed timeout aborts only a
# read that makes *no progress*, then MAX_RETRY rounds of backoff recover the
# transient stall. See diagnostics.md §1/§2.
_GDAL_HTTP_CONFIG = {
    'GDAL_HTTP_CONNECTTIMEOUT': '20',
    'GDAL_HTTP_LOW_SPEED_TIME': '60',
    'GDAL_HTTP_LOW_SPEED_LIMIT': '1',
    'GDAL_HTTP_MAX_RETRY': '5',
    'GDAL_HTTP_RETRY_DELAY': '1',
    'CPL_VSIL_CURL_USE_HEAD': 'NO',
}
for _k, _v in _GDAL_HTTP_CONFIG.items():
    os.environ.setdefault(_k, _v)

_rio_configured = False


def _configure_rio():
    global _rio_configured
    if not _rio_configured:
        import odc.stac
        odc.stac.configure_rio(cloud_defaults=True, aws={"aws_unsigned": True}, **_GDAL_HTTP_CONFIG)
        _rio_configured = True


def _solar_day(item) -> str:
    """Solar day of a STAC item: UTC datetime shifted by its centre longitude."""
    lon = (item.bbox[0] + item.bbox[2]) / 2
    return str((item.datetime + timedelta(hours=lon / 15.0)).date())


def clean_dataset(ds: Dataset, sentinel2: Sentinel2 = defaultsentinel2,
                  max_nan_fraction: float = 0.5) -> Dataset:
    """Cloud-mask a raw cube window and drop too-cloudy frames.

    Same semantics as the old persisted "clean" zarr — fmask-based clear-sky
    mask (cloud + shadow pixels → NaN, fmask band dropped), then scenes whose
    NaN fraction exceeds ``max_nan_fraction`` are removed — but computed on
    read instead of stored, so the clean copy costs no disk.
    """
    fmask = ds[sentinel2.cloud_mask_band]
    clear_mask = (fmask != sentinel2.fmask_cloud) & (fmask != sentinel2.fmask_shadow)
    ds = ds.drop_vars(sentinel2.cloud_mask_band).where(clear_mask)

    nan_frac = ds.to_array().isnull().mean(dim=['variable', 'x', 'y'])
    ds = ds.sel(time=nan_frac < max_nan_fraction)

    ds = ds.rio.write_crs(sentinel2.crs, inplace=False)
    return ds.assign_attrs(max_nan_fraction=max_nan_fraction)


@frozen
class Cube:
    """The machine-wide Sentinel-2 store: one grid, one index, zero re-downloads.

    Composed from :class:`borevitz_lab.config.Config` (where the store
    lives) and :class:`pysentinel2.sentinel2.Sentinel2` (what to fetch and
    from where). No inheritance.

    Example:
        ```python
        from datetime import date
        from pysentinel2.cube import Cube

        cube = Cube()
        ds  = cube.get(bbox, date(2024, 1, 1), date(2024, 6, 30))  # fills gaps, returns raw window
        dsc = cube.get(bbox, date(2024, 1, 1), date(2024, 6, 30), clean=True)
        dsq = cube.get_query(query)        # same, for pipelines that speak Query
        ```
    """

    config: Config = default_config
    sentinel2: Sentinel2 = defaultsentinel2
    paths: Paths = field(init=False)

    paths.default(lambda s: Paths(s.config))

    def __attrs_post_init__(s):
        makedirs(s.paths.root, exist_ok=True)

    def _index(s) -> Index:
        return Index(s.paths.index_db)

    # -- fill -------------------------------------------------------------

    def fill(s, bbox: list[float], start: date, end: date, threads: int = 8) -> int:
        """Ensure every (day x chunk) cell covering ``bbox`` x ``[start, end]``
        is populated.

        Query-agnostic: takes the region and range directly, no
        :class:`borevitz_lab.query.Query` (and none of its registry/dir side
        effects) required. Returns the number of cells actually downloaded —
        0 means the request was already fully covered and no network was
        touched (beyond a STAC search if this exact region/range was never
        searched before).
        """
        window = grid.window_for_bbox(bbox)
        wanted = grid.chunks_in_window(window)
        bbox6933 = (grid.X0 + window[2] * grid.RES, grid.Y_TOP - window[1] * grid.RES,
                    grid.X0 + window[3] * grid.RES, grid.Y_TOP - window[0] * grid.RES)
        ix = s._index()
        try:
            if not ix.search_covered(bbox6933, start, end):
                s._search_stac(ix, window, bbox6933, start, end)

            by_day = ix.scenes_for_range(start, end, s.sentinel2.max_cloud_cover)
            downloaded = 0
            for day, item_dicts in by_day.items():
                missing = set(wanted) - ix.chunks_done(day, wanted)
                if not missing:
                    continue
                day_window = grid.window_of_chunks(sorted(missing))
                s._download_day(day, item_dicts, day_window, threads)
                ix.mark_chunks(day, grid.chunks_in_window(day_window))
                downloaded += len(missing)
            return downloaded
        finally:
            ix.close()

    def _search_stac(s, ix: Index, window, bbox6933, start: date, end: date) -> None:
        """STAC-search the window (no cloud filter — that's applied at read
        from the index, so a laxer threshold later needs no re-search)."""
        import pystac_client
        from urllib3 import Retry
        from pystac_client.stac_api_io import StacApiIO

        # DEA STAC's first request after a cold cache often hits a 504; the
        # next one succeeds. See diagnostics.md.
        stac_io = StacApiIO(max_retries=Retry(
            total=5,
            backoff_factor=1.0,
            status_forcelist=[408, 429, 502, 503, 504],
            allowed_methods=['GET', 'POST'],
        ))
        catalog = pystac_client.Client.open(s.sentinel2.stac_url, stac_io=stac_io)
        result = catalog.search(
            bbox=grid.window_to_4326(window),
            collections=s.sentinel2.collections,
            datetime=f'{start}/{end}',
        )
        items = list(result.items())
        ix.upsert_scenes([
            (item.id, _solar_day(item),
             item.properties.get('eo:cloud_cover'), item.to_dict())
            for item in items
        ])
        ix.record_search(bbox6933, start, end)

    def _download_day(s, day: str, item_dicts: list[dict], window, threads: int) -> None:
        """Fetch one solar day's pixels for ``window`` and write them into
        the day's global-grid arrays. Chunk-aligned, so writes are whole
        chunks and the index rows marked afterwards are truthful."""
        import odc.stac
        import pystac
        _configure_rio()

        items = [pystac.Item.from_dict(d) for d in item_dicts]
        # In-process threaded scheduler — no distributed cluster (see the
        # deadlock notes in diagnostics.md).
        with dask.config.set(scheduler='threads', num_workers=threads):
            ds: Dataset = odc.stac.load(
                items,
                bands=s.sentinel2.bands,
                geobox=grid.geobox_for_window(window),
                groupby=s.sentinel2.groupby,
                chunks={'time': 1, 'x': grid.CHUNK, 'y': grid.CHUNK},
                # One corrupt DEA tile costs a nodata gap, not the whole day.
                fail_on_error=False,
            ).compute()

        if ds.time.size == 0:
            data_slice = None
        else:
            # Items were grouped by solar day before the call, so expect one
            # timestep; collapse defensively if odc still yields several.
            data_slice = ds.isel(time=0) if ds.time.size == 1 else ds.max(dim='time', keep_attrs=True)

        root = zarr.open_group(s.paths.store, mode='a')
        try:
            day_group = root[day]
        except KeyError:
            day_group = root.create_group(day)
        row0, row1, col0, col1 = window
        for band in s.sentinel2.bands:
            da = data_slice[band] if data_slice is not None else None
            dtype = da.dtype if da is not None else np.dtype('int16')
            nodata = (da.attrs.get('nodata') if da is not None else None)
            if nodata is None:
                nodata = 0 if dtype.kind == 'u' else -999
            try:
                arr = day_group[band]
            except KeyError:
                arr = day_group.create_array(
                    band,
                    shape=(grid.HEIGHT_PX, grid.WIDTH_PX),
                    chunks=(grid.CHUNK, grid.CHUNK),
                    dtype=dtype,
                    fill_value=nodata,
                )
                arr.attrs['nodata'] = int(nodata)
            if da is not None:
                arr[row0:row1, col0:col1] = da.values
        day_group.attrs['crs'] = grid.CRS

    # -- read -------------------------------------------------------------

    def get(s, bbox: list[float], start: date, end: date, clean: bool = False,
            max_nan_fraction: float = 0.5, threads: int = 8) -> Dataset:
        """Return the Sentinel-2 window for ``bbox`` x ``[start, end]``,
        downloading only what's missing first.

        Query-agnostic — the data layer of the package. Pipelines that
        speak :class:`borevitz_lab.query.Query` use :meth:`get_query`.

        Args:
            bbox: ``[west, south, east, north]`` in EPSG:4326.
            start: Inclusive start date.
            end: Inclusive end date.
            clean: Apply :func:`clean_dataset` (cloud mask + frame filter)
                to the window before returning it.
            max_nan_fraction: Frame-filter threshold used when ``clean=True``.
            threads: I/O concurrency for any downloads triggered.

        Returns:
            xarray.Dataset with dims ``(time, y, x)`` on the fixed grid,
            time = solar days (cloud-filtered per the ``Sentinel2`` config).
        """
        s.fill(bbox, start, end, threads=threads)
        window = grid.window_for_bbox(bbox)
        ix = s._index()
        try:
            by_day = ix.scenes_for_range(start, end, s.sentinel2.max_cloud_cover)
        finally:
            ix.close()

        ds = s._read_window(window, sorted(by_day))
        if clean:
            ds = clean_dataset(ds, s.sentinel2, max_nan_fraction)
        return ds

    # -- Query adapters (the reproducibility layer speaks Query) ----------

    def fill_query(s, query, threads: int = 8) -> int:
        """:meth:`fill` for a :class:`borevitz_lab.query.Query`."""
        return s.fill(query.bbox, query.start, query.end, threads=threads)

    def get_query(s, query, clean: bool = False, max_nan_fraction: float = 0.5,
                  threads: int = 8) -> Dataset:
        """:meth:`get` for a :class:`borevitz_lab.query.Query`."""
        return s.get(query.bbox, query.start, query.end, clean=clean,
                     max_nan_fraction=max_nan_fraction, threads=threads)

    def _read_window(s, window, days: list[str]) -> Dataset:
        row0, row1, col0, col1 = window
        y, x = grid.coords_for_window(window)
        root = zarr.open_group(s.paths.store, mode='r') if days else None

        bands: dict[str, list[np.ndarray]] = {b: [] for b in s.sentinel2.bands}
        nodata: dict[str, int] = {}
        kept_days = []
        for day in days:
            try:
                day_group = root[day]
            except KeyError:
                continue  # a searched day with no scenes written
            kept_days.append(day)
            for band in s.sentinel2.bands:
                arr = day_group[band]
                nodata[band] = arr.attrs.get('nodata', arr.fill_value)
                bands[band].append(arr[row0:row1, col0:col1])

        time = np.array([np.datetime64(d) for d in kept_days], dtype='datetime64[ns]')
        data_vars = {
            band: xr.DataArray(
                np.stack(stack) if stack else np.empty((0, len(y), len(x))),
                dims=('time', 'y', 'x'),
                attrs={'nodata': nodata.get(band)} if band in nodata else {},
            )
            for band, stack in bands.items()
        }
        ds = Dataset(data_vars, coords={'time': time, 'y': y, 'x': x})
        import rioxarray  # noqa: F401 — registers the .rio accessor
        return ds.rio.write_crs(grid.CRS, inplace=False).assign_attrs(
            read_at=datetime.now(timezone.utc).isoformat()
        )


# -- offline tests (synthetic writes, no network) --------------------------

def _tmp_cube() -> Cube:
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix='pysentinel2_cube_test_')
    return Cube(config=Config(out_dir=tmpdir, tmp_dir=tmpdir))


_TEST_BBOX = [148.36265, -33.52606, 148.38265, -33.50606]


def _write_synthetic_day(cube: Cube, day: str, window, value: int) -> None:
    """Write a constant-valued synthetic day directly, bypassing the network."""
    root = zarr.open_group(cube.paths.store, mode='a')
    try:
        day_group = root[day]
    except KeyError:
        day_group = root.create_group(day)
    row0, row1, col0, col1 = window
    for band in cube.sentinel2.bands:
        dtype = np.dtype('uint8') if band == cube.sentinel2.cloud_mask_band else np.dtype('int16')
        fill = 0 if dtype.kind == 'u' else -999
        try:
            arr = day_group[band]
        except KeyError:
            arr = day_group.create_array(
                band, shape=(grid.HEIGHT_PX, grid.WIDTH_PX),
                chunks=(grid.CHUNK, grid.CHUNK), dtype=dtype, fill_value=fill,
            )
            arr.attrs['nodata'] = fill
        band_value = 1 if band == cube.sentinel2.cloud_mask_band else value  # fmask 1 = clear
        arr[row0:row1, col0:col1] = np.full((row1 - row0, col1 - col0), band_value, dtype=dtype)


_TEST_START, _TEST_END = date(2024, 1, 1), date(2024, 1, 31)


def _prime_synthetic(cube: Cube, day: str, item_id: str, value: int):
    """Mark a fully-populated synthetic day in the index and store."""
    window = grid.window_for_bbox(_TEST_BBOX)
    ix = cube._index()
    ix.upsert_scenes([(item_id, day, 1.0, {'id': item_id})])
    ix.record_search((-1e9, -1e9, 1e9, 1e9), _TEST_START, _TEST_END)
    ix.mark_chunks(day, grid.chunks_in_window(window))
    ix.close()
    _write_synthetic_day(cube, day, window, value=value)


def test_synthetic_write_read_roundtrip():
    cube = _tmp_cube()
    _prime_synthetic(cube, '2024-01-03', 'synth_a', value=1234)
    ds = cube.get(_TEST_BBOX, _TEST_START, _TEST_END)
    return (
        ds.time.size == 1
        and int(ds['nbart_red'].isel(time=0)[0, 0]) == 1234
        and ds.rio.crs is not None
    )


def test_clean_masks_and_drops_fmask():
    cube = _tmp_cube()
    _prime_synthetic(cube, '2024-01-08', 'synth_b', value=42)
    ds = cube.get(_TEST_BBOX, _TEST_START, _TEST_END, clean=True)
    return (
        cube.sentinel2.cloud_mask_band not in ds.data_vars
        and ds.time.size == 1  # fully clear frame survives the filter
        and float(ds['nbart_red'].isel(time=0)[0, 0]) == 42.0
    )


def test_fill_skips_populated_cells():
    """With search covered and all chunks marked, fill() must return 0
    without touching the network (no STAC client is even constructed)."""
    cube = _tmp_cube()
    _prime_synthetic(cube, '2024-01-13', 'synth_c', value=7)
    return cube.fill(_TEST_BBOX, _TEST_START, _TEST_END) == 0


def test_query_adapters_match_agnostic_calls():
    """get_query/fill_query are pure delegations to the bbox+dates core."""
    from borevitz_lab.query import Query
    cube = _tmp_cube()
    _prime_synthetic(cube, '2024-01-18', 'synth_d', value=99)
    q = Query(bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
              stub='cube_adapter', config=cube.config)
    ds_agnostic = cube.get(_TEST_BBOX, _TEST_START, _TEST_END)
    ds_query = cube.get_query(q)
    return (
        cube.fill_query(q) == 0
        and ds_query.time.size == ds_agnostic.time.size
        and float(ds_query['nbart_red'].isel(time=0)[0, 0])
            == float(ds_agnostic['nbart_red'].isel(time=0)[0, 0])
    )


def test():
    return all([
        test_synthetic_write_read_roundtrip(),
        test_clean_masks_and_drops_fmask(),
        test_fill_skips_populated_cells(),
        test_query_adapters_match_agnostic_calls(),
    ])


if __name__ == '__main__':
    print(test())
