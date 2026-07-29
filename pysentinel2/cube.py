"""One machine-wide Sentinel-2 datacube that fills itself on demand.

Every pixel this machine ever downloads lands in a single sparse store
on the fixed EPSG:6933/10 m grid (:mod:`pysentinel2.grid`):

    {config.tmp_dir}/sentinel2_cube/
    ├── index.db      # what's populated / what's been searched (pysentinel2.index)
    └── cube.zarr/
        └── 2024-01-03/   # one group per solar day, one array per band
            ├── nbart_red # global-grid array; only written chunks exist on disk
            └── ...

``Cube.get_ds(bbox, start, end)`` diffs the requested (day x chunk) cells
against the index, downloads only the missing cells, then reads the
window. Nothing is ever fetched twice — overlapping bboxes, extended
date ranges and repeat runs all reuse the same chunks. Only the raw
bands (incl. fmask) are stored; ``get_ds(..., clean=True)`` applies cloud
masking on read, so no second "clean" copy exists on disk.

The core API is query-agnostic (bbox + dates — the data layer);
``get_ds_query`` / ``fill_query`` adapt a :class:`borevitz_lab.query.Query`
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
                  max_cloud_fraction: float = 0.5,
                  min_valid_fraction: float = 0.2,
                  mask_snow: bool = True,
                  mask_water: bool = False,
                  buffer_px: int = 5) -> Dataset:
    """Cloud-mask a raw cube window and drop unusable frames.

    Three pixel roles are kept distinct, because they mean different
    things for frame quality:

    - *invalid* (fmask nodata — outside the scene footprint or never
      sensed) counts against ``valid_fraction`` only;
    - *contaminated* (cloud and shadow) drives both the per-pixel mask
      and the ``cloud_fraction`` frame gate;
    - *snow and water* are correctly classified surface states, not
      sensing failures: they are masked per pixel (snow by default,
      water on request) but never count toward the frame gate — a
      cloud-free, fully snow-covered frame survives as an observation
      (analysis/ finding 2: the previous conflation dropped 16 of 37
      usable alpine winter frames).

    Frames are dropped when ``cloud_fraction`` (cloud+shadow share of
    the valid pixels) exceeds ``max_cloud_fraction`` or
    ``valid_fraction`` (valid share of the window) falls below
    ``min_valid_fraction``.

    Cloud and shadow are dilated by ``buffer_px`` (circular structuring
    element) before masking, excluding the bright halo and penumbra that
    fmask misses; the measured edge bias supports the default of 5 px
    (analysis/ finding 3). Snow and water are masked without dilation.

    Per-frame statistics attach as ``time`` coordinates
    (``cloud_fraction``, ``valid_fraction``, ``snow_fraction``) so
    downstream consumers can see why frames survived and filter on snow
    explicitly. Computed on read, never stored.
    """
    fmask = ds[sentinel2.cloud_mask_band].values
    time_dims = ds[sentinel2.cloud_mask_band].dims

    valid = fmask != sentinel2.fmask_nodata
    gate = np.isin(fmask, [sentinel2.fmask_cloud, sentinel2.fmask_shadow])

    if buffer_px:
        from scipy.ndimage import binary_dilation
        yy, xx = np.ogrid[-buffer_px:buffer_px + 1, -buffer_px:buffer_px + 1]
        disk = (yy ** 2 + xx ** 2) <= buffer_px ** 2
        gate = np.stack([binary_dilation(frame, structure=disk) for frame in gate])

    gate &= valid                      # nodata is invalid, not "cloudy"
    snow = (fmask == sentinel2.fmask_snow) & valid

    masked_classes = [sentinel2.fmask_cloud, sentinel2.fmask_shadow]
    mask = gate
    if mask_snow:
        mask = mask | snow
        masked_classes.append(sentinel2.fmask_snow)
    if mask_water:
        mask = mask | ((fmask == sentinel2.fmask_water) & valid)
        masked_classes.append(sentinel2.fmask_water)
    clear = valid & ~mask

    n_valid = valid.sum(axis=(1, 2))
    valid_frac = n_valid / (valid.shape[1] * valid.shape[2])
    cloud_frac = np.where(n_valid > 0, gate.sum(axis=(1, 2)) / np.maximum(n_valid, 1), 1.0)
    snow_frac = np.where(n_valid > 0, snow.sum(axis=(1, 2)) / np.maximum(n_valid, 1), 0.0)
    keep = (valid_frac >= min_valid_fraction) & (cloud_frac <= max_cloud_fraction)

    nodatas = {name: ds[name].attrs.get('nodata') for name in ds.data_vars
               if name != sentinel2.cloud_mask_band}
    clear_da = xr.DataArray(clear, dims=time_dims)
    ds = ds.drop_vars(sentinel2.cloud_mask_band).where(clear_da)
    for name, nodata in nodatas.items():
        if nodata is not None:
            ds[name] = ds[name].where(ds[name] != nodata)

    ds = ds.assign_coords(
        valid_fraction=('time', np.round(valid_frac, 4)),
        cloud_fraction=('time', np.round(cloud_frac, 4)),
        snow_fraction=('time', np.round(snow_frac, 4)),
    )
    ds = ds.isel(time=np.flatnonzero(keep))

    ds = ds.rio.write_crs(sentinel2.crs, inplace=False)
    return ds.assign_attrs(
        max_cloud_fraction=max_cloud_fraction,
        min_valid_fraction=min_valid_fraction,
        cloud_buffer_px=buffer_px,
        frame_gate_classes=[sentinel2.fmask_cloud, sentinel2.fmask_shadow],
        masked_fmask_classes=sorted(masked_classes),
    )


def _chunks_passing_screen(fm: np.ndarray, window, sentinel2: Sentinel2,
                           threshold: float | None = None) -> list[tuple[int, int]]:
    """Chunk ids in ``window`` whose fmask justifies a reflectance download.

    A chunk passes when it has valid pixels and its cloud+shadow share
    of them is at most ``sentinel2.screen_cloud_fraction``. The decision
    is a pure function of the chunk's own fmask, so it is deterministic:
    a screened-out cell would be screened out identically by any later
    query, and the index can mark it done.
    """
    thr = sentinel2.screen_cloud_fraction if threshold is None else threshold
    row0, _, col0, _ = window
    passing = []
    for cy, cx in grid.chunks_in_window(window):
        r, c = cy * grid.CHUNK - row0, cx * grid.CHUNK - col0
        sub = fm[r:r + grid.CHUNK, c:c + grid.CHUNK]
        valid = sub != sentinel2.fmask_nodata
        n_valid = int(valid.sum())
        if n_valid == 0:
            continue
        bad = np.isin(sub, [sentinel2.fmask_cloud, sentinel2.fmask_shadow]) & valid
        if bad.sum() / n_valid <= thr:
            passing.append((cy, cx))
    return passing


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
        ds  = cube.get_ds(bbox, date(2024, 1, 1), date(2024, 6, 30))  # fills gaps, returns raw window
        dsc = cube.get_ds(bbox, date(2024, 1, 1), date(2024, 6, 30), clean=True)
        dsq = cube.get_ds_query(query)        # same, for pipelines that speak Query
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

    def _load_day(s, items, bands, window, threads):
        """One solar day's pixels for ``bands`` over ``window``, collapsed
        to a single timestep (None if the day yields no data)."""
        import odc.stac
        # In-process threaded scheduler — no distributed cluster (see the
        # deadlock notes in diagnostics.md).
        with dask.config.set(scheduler='threads', num_workers=threads):
            ds: Dataset = odc.stac.load(
                items,
                bands=bands,
                geobox=grid.geobox_for_window(window),
                groupby=s.sentinel2.groupby,
                chunks={'time': 1, 'x': grid.CHUNK, 'y': grid.CHUNK},
                # One corrupt DEA tile costs a nodata gap, not the whole day.
                fail_on_error=False,
            ).compute()
        if ds.time.size == 0:
            return None
        # Items were grouped by solar day before the call, so expect one
        # timestep; collapse defensively if odc still yields several.
        return ds.isel(time=0) if ds.time.size == 1 else ds.max(dim='time', keep_attrs=True)

    @staticmethod
    def _write_band(day_group, band, da, window, default_dtype, default_nodata):
        dtype = da.dtype if da is not None else default_dtype
        nodata = da.attrs.get('nodata') if da is not None else None
        if nodata is None:
            nodata = 0 if dtype.kind == 'u' else default_nodata
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
            row0, row1, col0, col1 = window
            arr[row0:row1, col0:col1] = da.values

    def _download_day(s, day: str, item_dicts: list[dict], window, threads: int) -> None:
        """Fetch one solar day for ``window`` and write it into the day's
        global-grid arrays, fmask first.

        Two phases, both chunk-aligned:

        1. Download the fmask band for the whole window and write it —
           the cheap band (it compresses ~50x) is always stored.
        2. Screen each chunk on its own fmask statistics
           (:func:`_chunks_passing_screen`); reflectance bands are
           downloaded only for the window of passing chunks. Screened
           chunks keep their real fmask and read as nodata reflectance,
           and ``clean_dataset`` drops such frames via the same fmask.

        Because the screen is deterministic per chunk, every chunk in
        ``window`` — passing or screened — is safely marked done by the
        caller.
        """
        import pystac
        _configure_rio()
        items = [pystac.Item.from_dict(d) for d in item_dicts]

        root = zarr.open_group(s.paths.store, mode='a')
        try:
            day_group = root[day]
        except KeyError:
            day_group = root.create_group(day)
        day_group.attrs['crs'] = grid.CRS

        fmask_band = s.sentinel2.cloud_mask_band
        fm_slice = s._load_day(items, [fmask_band], window, threads)
        fm = fm_slice[fmask_band] if fm_slice is not None else None
        s._write_band(day_group, fmask_band, fm, window, np.dtype('uint8'), 0)
        if fm is None:
            return

        passing = _chunks_passing_screen(fm.values, window, s.sentinel2)
        if not passing:
            return
        sub = grid.window_of_chunks(passing)
        bands = [b for b in s.sentinel2.bands if b != fmask_band]
        data = s._load_day(items, bands, sub, threads)
        for band in bands:
            da = data[band] if data is not None else None
            s._write_band(day_group, band, da, sub, np.dtype('int16'), -999)

    # -- read -------------------------------------------------------------

    def get_ds(s, bbox: list[float], start: date, end: date, clean: bool = False,
            indices: tuple[str, ...] = (), threads: int = 8,
            **clean_kwargs) -> Dataset:
        """Return the Sentinel-2 window for ``bbox`` x ``[start, end]``,
        downloading only what's missing first.

        Query-agnostic — the data layer of the package. Pipelines that
        speak :class:`borevitz_lab.query.Query` use :meth:`get_ds_query`.

        Args:
            bbox: ``[west, south, east, north]`` in EPSG:4326.
            start: Inclusive start date.
            end: Inclusive end date.
            clean: Apply :func:`clean_dataset` (cloud/shadow masking with
                dilation + frame filtering) to the window before returning it.
            indices: Spectral indices to compute on read (never stored) —
                any of ``'NDVI'``, ``'CFI'``, ``'NIRv'``, ``'NDTI'``,
                ``'CAI'`` (:mod:`pysentinel2.derive`). Requesting indices
                implies ``clean=True`` so formulas see cloud-masked
                reflectance.
            threads: I/O concurrency for any downloads triggered.
            **clean_kwargs: Forwarded to :func:`clean_dataset` —
                ``max_cloud_fraction``, ``min_valid_fraction``,
                ``mask_snow``, ``mask_water``, ``buffer_px``.

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
        if clean or indices:
            ds = clean_dataset(ds, s.sentinel2, **clean_kwargs)
        if indices:
            from pysentinel2.derive import add_indices
            ds = add_indices(ds, indices)
        return ds

    # -- Query adapters (the reproducibility layer speaks Query) ----------

    def fill_query(s, query, threads: int = 8) -> int:
        """:meth:`fill` for a :class:`borevitz_lab.query.Query`."""
        return s.fill(query.bbox, query.start, query.end, threads=threads)

    def get_ds_query(s, query, clean: bool = False,
                  indices: tuple[str, ...] = (), threads: int = 8,
                  **clean_kwargs) -> Dataset:
        """:meth:`get_ds` for a :class:`borevitz_lab.query.Query`."""
        return s.get_ds(query.bbox, query.start, query.end, clean=clean,
                     indices=indices, threads=threads, **clean_kwargs)

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
                try:
                    arr = day_group[band]
                except KeyError:
                    # A screened day: fmask was stored, reflectance was
                    # never fetched — read as nodata.
                    is_mask = band == s.sentinel2.cloud_mask_band
                    fill = 0 if is_mask else -999
                    nodata.setdefault(band, fill)
                    bands[band].append(np.full(
                        (row1 - row0, col1 - col0), fill,
                        dtype='uint8' if is_mask else 'int16'))
                    continue
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
    ds = cube.get_ds(_TEST_BBOX, _TEST_START, _TEST_END)
    return (
        ds.time.size == 1
        and int(ds['nbart_red'].isel(time=0)[0, 0]) == 1234
        and ds.rio.crs is not None
    )


def test_clean_masks_and_drops_fmask():
    cube = _tmp_cube()
    _prime_synthetic(cube, '2024-01-08', 'synth_b', value=42)
    ds = cube.get_ds(_TEST_BBOX, _TEST_START, _TEST_END, clean=True)
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


def _synthetic_window(fmask: np.ndarray, band_value: int = 4000) -> Dataset:
    """In-memory raw window with a given (time, y, x) fmask array."""
    shape = fmask.shape
    data_vars = {
        band: xr.DataArray(np.full(shape, band_value, dtype='int16'),
                           dims=('time', 'y', 'x'), attrs={'nodata': -999})
        for band in defaultsentinel2.bands if band != defaultsentinel2.cloud_mask_band
    }
    data_vars[defaultsentinel2.cloud_mask_band] = xr.DataArray(
        fmask.astype('uint8'), dims=('time', 'y', 'x'), attrs={'nodata': 0})
    times = np.array([np.datetime64('2024-01-01') + np.timedelta64(5 * i, 'D')
                      for i in range(shape[0])])
    return Dataset(data_vars, coords={
        'time': times,
        'y': np.arange(shape[1], dtype='float64'),
        'x': np.arange(shape[2], dtype='float64'),
    })


def test_clean_drops_cloudy_frame_keeps_clear():
    fmask = np.ones((2, 40, 40))
    fmask[1, :, :24] = 2                      # frame 1: 60% cloud
    ds = clean_dataset(_synthetic_window(fmask), buffer_px=0)
    return ds.time.size == 1 and float(ds.cloud_fraction[0]) == 0.0


def test_nodata_margin_does_not_drop_frame():
    """A clear frame with a 60% off-swath margin must survive — the old
    whole-window NaN rule would have dropped it."""
    fmask = np.ones((1, 40, 40))
    fmask[0, :, :24] = 0                      # 60% outside the footprint
    ds = clean_dataset(_synthetic_window(fmask), buffer_px=0)
    return (
        ds.time.size == 1
        and abs(float(ds.valid_fraction[0]) - 0.4) < 1e-6
        and float(ds.cloud_fraction[0]) == 0.0
    )


def test_sliver_frame_dropped():
    fmask = np.zeros((1, 40, 40))
    fmask[0, :, :4] = 1                       # only 10% of the window sensed
    ds = clean_dataset(_synthetic_window(fmask), buffer_px=0)
    return ds.time.size == 0


def test_cloud_buffer_dilates():
    """Pixels adjacent to a cloud must be masked when buffer_px > 0."""
    fmask = np.ones((1, 40, 40))
    fmask[0, 20, 20] = 2                      # single cloud pixel
    ds0 = clean_dataset(_synthetic_window(fmask), buffer_px=0)
    ds3 = clean_dataset(_synthetic_window(fmask), buffer_px=3)
    neighbour0 = float(ds0['nbart_red'][0, 20, 22])
    neighbour3 = ds3['nbart_red'][0, 20, 22]
    return neighbour0 == 4000.0 and bool(np.isnan(neighbour3))


def test_snow_masked_water_kept_by_default():
    fmask = np.ones((1, 40, 40))
    fmask[0, 0, 0] = defaultsentinel2.fmask_snow
    fmask[0, 0, 1] = defaultsentinel2.fmask_water
    ds = clean_dataset(_synthetic_window(fmask), buffer_px=0)
    snow_px = ds['nbart_red'][0, 0, 0]
    water_px = ds['nbart_red'][0, 0, 1]
    ds_w = clean_dataset(_synthetic_window(fmask), buffer_px=0, mask_water=True)
    return (
        bool(np.isnan(snow_px)) and float(water_px) == 4000.0
        and bool(np.isnan(ds_w['nbart_red'][0, 0, 1]))
    )


def test_snow_excluded_from_frame_gate():
    """A cloud-free frame that is 75% snow must survive the gate with its
    snow pixels masked and its snow_fraction reported (analysis/ finding 2)."""
    fmask = np.ones((1, 40, 40))
    fmask[0, :30, :] = defaultsentinel2.fmask_snow
    ds = clean_dataset(_synthetic_window(fmask), buffer_px=0)
    return (
        ds.time.size == 1
        and float(ds.cloud_fraction[0]) == 0.0
        and abs(float(ds.snow_fraction[0]) - 0.75) < 1e-6
        and bool(np.isnan(ds['nbart_red'][0, 0, 0]))       # snow pixel masked
        and float(ds['nbart_red'][0, 35, 0]) == 4000.0     # clear pixel kept
    )


def test_screen_passes_clear_rejects_cloudy_chunks():
    """Per-chunk download screen: clear and off-swath chunks decided
    correctly, threshold honoured."""
    from pysentinel2.cube import _chunks_passing_screen
    win = (0, grid.CHUNK, 0, 3 * grid.CHUNK)                 # 1 x 3 chunks
    fm = np.ones((grid.CHUNK, 3 * grid.CHUNK), dtype='uint8')
    fm[:, grid.CHUNK:2 * grid.CHUNK] = defaultsentinel2.fmask_cloud   # chunk 1 fully cloudy
    fm[:, 2 * grid.CHUNK:] = defaultsentinel2.fmask_nodata            # chunk 2 off-swath
    passing = _chunks_passing_screen(fm, win, defaultsentinel2)
    all_pass = _chunks_passing_screen(fm, win, defaultsentinel2, threshold=1.0)
    return passing == [(0, 0)] and all_pass == [(0, 0), (0, 1)]


def test_query_adapters_match_agnostic_calls():
    """get_ds_query/fill_query are pure delegations to the bbox+dates core."""
    from borevitz_lab.query import Query
    cube = _tmp_cube()
    _prime_synthetic(cube, '2024-01-18', 'synth_d', value=99)
    q = Query(bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
              stub='cube_adapter', config=cube.config)
    ds_agnostic = cube.get_ds(_TEST_BBOX, _TEST_START, _TEST_END)
    ds_query = cube.get_ds_query(q)
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
        test_clean_drops_cloudy_frame_keeps_clear(),
        test_nodata_margin_does_not_drop_frame(),
        test_sliver_frame_dropped(),
        test_cloud_buffer_dilates(),
        test_snow_masked_water_kept_by_default(),
        test_snow_excluded_from_frame_gate(),
        test_screen_passes_clear_rejects_cloudy_chunks(),
        test_query_adapters_match_agnostic_calls(),
    ])


if __name__ == '__main__':
    print(test())
