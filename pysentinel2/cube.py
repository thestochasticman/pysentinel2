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
_rio_lock = __import__('threading').Lock()


def _configure_rio():
    global _rio_configured
    with _rio_lock:
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
    import rioxarray  # noqa: F401 — registers the .rio accessor used below
    fmask = ds[sentinel2.cloud_mask_band].values
    time_dims = ds[sentinel2.cloud_mask_band].dims

    valid = fmask != sentinel2.fmask_nodata
    gate = np.isin(fmask, [sentinel2.fmask_cloud, sentinel2.fmask_shadow])

    if buffer_px:
        # cv2.dilate over uint8 frames: same circular structuring element
        # as the previous scipy binary_dilation, ~15x faster on multi-year
        # windows (SIMD + no per-frame python overhead in the kernel).
        import cv2
        yy, xx = np.ogrid[-buffer_px:buffer_px + 1, -buffer_px:buffer_px + 1]
        disk = ((yy ** 2 + xx ** 2) <= buffer_px ** 2).astype(np.uint8)
        gate = np.stack([
            cv2.dilate(frame.astype(np.uint8), disk).astype(bool)
            for frame in gate
        ])

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

    # Masking preserves the storage dtype: masked pixels are set to the
    # band's nodata sentinel in the int16 array itself, after dropped
    # frames are discarded. The previous implementation used ``.where``
    # which promotes every band to float — a multi-year window ballooned
    # from ~3.5 GB int16 to ~7-14 GB float and pushed 8 GB machines into
    # swap, slowing every downstream consumer by 5-20x. Consumers convert
    # to float per batch via :func:`to_float`.
    #
    # The input dataset is CONSUMED: each raw band is deleted as soon as
    # its cleaned copy exists, so raw + cleaned never fully coexist —
    # that transient doubling (~9 GB for 8 years) was itself a swap trip.
    keep_idx = np.flatnonzero(keep)
    not_clear = ~clear[keep_idx]
    band_names = [n for n in ds.data_vars if n != sentinel2.cloud_mask_band]
    time_coord = ds.coords['time'].values[keep_idx]
    y_coord, x_coord = ds.coords['y'], ds.coords['x']

    new_vars = {}
    for name in band_names:
        da = ds[name]
        nodata = da.attrs.get('nodata')
        if nodata is None:
            nodata = -999
        vals = da.values[keep_idx]           # one band's kept frames
        vals[not_clear] = nodata
        new_vars[name] = xr.DataArray(
            vals, dims=('time', 'y', 'x'), attrs={**da.attrs, 'nodata': nodata})
        del ds[name]                          # free the raw band now

    ds = xr.Dataset(new_vars, coords={'time': time_coord, 'y': y_coord, 'x': x_coord})
    ds = ds.assign_coords(
        valid_fraction=('time', np.round(valid_frac[keep_idx], 4)),
        cloud_fraction=('time', np.round(cloud_frac[keep_idx], 4)),
        snow_fraction=('time', np.round(snow_frac[keep_idx], 4)),
    )

    ds = ds.rio.write_crs(sentinel2.crs, inplace=False)
    return ds.assign_attrs(
        max_cloud_fraction=max_cloud_fraction,
        min_valid_fraction=min_valid_fraction,
        cloud_buffer_px=buffer_px,
        frame_gate_classes=[sentinel2.fmask_cloud, sentinel2.fmask_shadow],
        masked_fmask_classes=sorted(masked_classes),
    )


def to_float(ds: Dataset, dtype: str = 'float32') -> Dataset:
    """Cast a cleaned window's bands to float, nodata sentinel -> NaN.

    :func:`clean_dataset` keeps bands in their compact storage dtype
    (int16, masked pixels at the band's ``nodata`` attr) so a multi-year
    window fits in RAM. Consumers that need float-with-NaN semantics
    call this at point of use — ideally on a time slice or batch, not
    the whole window, to keep the float footprint bounded.
    """
    out = ds.copy(deep=False)
    for name in ds.data_vars:
        da = ds[name]
        nodata = da.attrs.get('nodata')
        vals = da.values.astype(dtype)
        if nodata is not None:
            vals[da.values == nodata] = np.nan
        out[name] = da.copy(data=vals)
        out[name].attrs.pop('nodata', None)
    return out


def _plan_day_windows(window, covered,
                      min_fill: float = 0.9, max_rects: int = 8):
    """Load windows for one day: the missing region as its own rects, or
    their bbox when that is barely bigger.

    Returns ``[]`` when the day is fully covered. When the missing rects
    fill at least ``min_fill`` of their bounding box (or fragment past
    ``max_rects``), one bbox load wins. Otherwise (a ring or L around
    covered data — a district requested around an already-downloaded
    farm) each missing rect loads separately, so covered pixels are
    never re-fetched. The threshold is deliberately high: loads batch
    across days, so an extra rect costs a handful of load calls per
    ~60-day batch while the byte saving applies to every day.
    """
    missing = grid.rect_subtract(window, covered)
    if not missing:
        return []
    bbox = grid.rects_bbox(missing)
    bbox_area = (bbox[1] - bbox[0]) * (bbox[3] - bbox[2])
    missing_area = sum((r1 - r0) * (c1 - c0) for r0, r1, c0, c1 in missing)
    if missing_area / bbox_area >= min_fill or len(missing) > max_rects:
        return [bbox]
    return missing


def _window_wants_reflectance(fm: np.ndarray, sentinel2: Sentinel2,
                              threshold: float | None = None) -> bool:
    """Would the removed download screen have fetched reflectance here?

    True when the window has valid pixels and their cloud+shadow share
    is at most ``sentinel2.screen_cloud_fraction``. Used only by
    :meth:`Cube.repair` to recognise fmask-only days written by old
    screening-era fills as legitimately reflectance-free rather than
    download-failure debris.
    """
    thr = sentinel2.screen_cloud_fraction if threshold is None else threshold
    valid = fm != sentinel2.fmask_nodata
    n_valid = int(valid.sum())
    if n_valid == 0:
        return False
    bad = np.isin(fm, [sentinel2.fmask_cloud, sentinel2.fmask_shadow]) & valid
    return bad.sum() / n_valid <= thr


_FAILED_NODATA_FRACTION = 0.5


def _reflectance_looks_failed(data: Dataset, bands,
                              sentinel2: Sentinel2) -> bool:
    """True when a day's freshly loaded reflectance is mostly nodata on
    pixels its own fmask calls valid — the signature of failed HTTP reads
    (``fail_on_error=False``) rather than a real acquisition gap. A single
    corrupt tile stays below the threshold and is accepted as a small gap,
    matching the long-standing behaviour. ``data`` is the day's
    window-local dataset including the fmask band."""
    fm = data[sentinel2.cloud_mask_band].values
    valid = fm != sentinel2.fmask_nodata
    n_valid = int(valid.sum())
    if n_valid == 0:
        return False
    for band in bands:
        da = data[band]
        nodata = da.attrs.get('nodata', -999)
        vals = da.values
        bad = ((vals == nodata) | np.isnan(vals)) & valid
        if bad.sum() / n_valid > _FAILED_NODATA_FRACTION:
            return True
    return False


def _day_group(store, day: str):
    """The zarr group for ``day``, created on first use.

    Opened by path with ``use_consolidated=False``: child-group lookups
    through a parent Group object honour any ``consolidated_metadata``
    embedded in that group's own metadata document, and a stale snapshot
    (e.g. from a one-off ``zarr.consolidate_metadata`` call) then hides
    arrays that exist on disk — lookups KeyError, the create fallback
    hits ContainsArrayError, and fills crash. Path-based opens with the
    flag off are immune regardless of what metadata a store carries.
    """
    group = zarr.open_group(store, path=day, mode='a', use_consolidated=False)
    group.attrs['crs'] = grid.CRS
    return group


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

    def fill(s, bbox: list[float], start: date, end: date, threads: int = 16,
             batch_days: int = 64) -> int:
        """Ensure every pixel of ``bbox``'s tight window is populated for
        every candidate solar day in ``[start, end]``.

        Query-agnostic: takes the region and range directly, no
        :class:`borevitz_lab.query.Query` (and none of its registry/dir side
        effects) required. Coverage accounting is pixel-exact: each day's
        missing region is the tight window minus its recorded coverage
        rects, so small farms pay no chunk padding (the previous 256 px
        chunk accounting padded a 2.5 km bbox by 3.6x). Returns the number
        of (day x pixel) cells downloaded — 0 means fully covered, no
        network touched (beyond a STAC search if this exact region/range
        was never searched before).

        The unit of network work is a multi-day *batch*, not a day: one
        bulk ``odc.stac.load`` of every band (fmask included) per up to
        ``batch_days`` days, so the dask graph carries days x chunks x
        bands tasks and ``threads`` I/O workers stay saturated across day
        boundaries. (An earlier design fetched fmask first and screened
        chunks before reflectance; measured on this store the screen
        skipped 8.8% of days' reflectance while paying an extra COG-open
        round on every day — a net loss in a latency-bound phase.)

        Coverage is marked per completed day, strictly after its pixels
        are on disk — an interrupted fill never records an unwritten
        region (and never re-downloads a written one).
        """
        window = grid.tight_window_for_bbox(bbox)
        bbox6933 = (grid.X0 + window[2] * grid.RES, grid.Y_TOP - window[1] * grid.RES,
                    grid.X0 + window[3] * grid.RES, grid.Y_TOP - window[0] * grid.RES)
        ix = s._index()
        try:
            if not ix.search_covered(bbox6933, start, end):
                s._search_stac(ix, window, bbox6933, start, end)

            by_day = ix.scenes_for_range(start, end, s.sentinel2.max_cloud_cover)
            # Work grouped by load window: a day whose missing region is a
            # ring/L around covered data contributes one entry per missing
            # rect (never re-fetching the covered middle); a mostly-missing
            # day contributes its single bbox. Days sharing a window load
            # together, so the ring case costs at most 4 group loads.
            groups: dict[tuple, dict[str, list]] = {}
            for day, item_dicts in by_day.items():
                for win in _plan_day_windows(window, ix.covered_rects(day)):
                    groups.setdefault(win, {})[day] = item_dicts
            if not groups:
                return 0

            _configure_rio()
            zarr.open_group(s.paths.store, mode='a', use_consolidated=False)
            fmask_band = s.sentinel2.cloud_mask_band
            all_bands = list(s.sentinel2.bands)
            refl_bands = [b for b in all_bands if b != fmask_band]

            # Single pass: every band (fmask included) for a batch of days
            # in one bulk load per window group. (An earlier fmask-first
            # screen was removed after measurement: 8.8% of days' bytes
            # saved for an extra request round on every day.)
            downloaded = 0
            failed = []
            for win, day_items in groups.items():
                # A batch is held in RAM while it's split and written, so
                # cap its day count by window area: ~256 MB of int16 in
                # flight — the lab's field machines have 8 GB RAM, and
                # macOS answers memory pressure with SIGKILL.
                n_px = (win[1] - win[0]) * (win[3] - win[2])
                win_area = n_px
                day_batch = max(1, min(batch_days,
                                       int(256e6 / max(len(all_bands) * 2 * n_px, 1))))
                days_sorted = sorted(day_items)
                for i in range(0, len(days_sorted), day_batch):
                    batch = days_sorted[i:i + day_batch]
                    data_by_day = s._load_days(
                        {d: day_items[d] for d in batch}, all_bands, win, threads)
                    for day in batch:
                        data = data_by_day.get(day)
                        day_group = _day_group(s.paths.store, day)
                        if data is None:
                            # A searched day whose items yielded no pixels:
                            # store the empty fmask so reads see nodata, mark.
                            s._write_band(day_group, fmask_band, None, win,
                                          np.dtype('uint8'), 0)
                            ix.mark_rect(day, win)
                            downloaded += win_area
                            continue
                        # Integrity gate: fail_on_error=False turns failed
                        # HTTP reads into silent nodata. A band mostly nodata
                        # where the day's own fmask has valid ground is a
                        # failed download, not a data gap — writing and
                        # marking it would poison the store permanently.
                        # Leave it unmarked (and unwritten) for retry.
                        if _reflectance_looks_failed(data, refl_bands, s.sentinel2):
                            failed.append((day, win))
                            continue
                        s._write_band(day_group, fmask_band, data[fmask_band],
                                      win, np.dtype('uint8'), 0)
                        for band in refl_bands:
                            s._write_band(day_group, band, data[band], win,
                                          np.dtype('int16'), -999)
                        ix.mark_rect(day, win)
                        downloaded += win_area
            if failed:
                print(f'pysentinel2: {len(failed)} (day, window) load(s) failed '
                      f'download integrity and were left unmarked for retry: '
                      f'{failed[:3]}{"..." if len(failed) > 3 else ""}')
            return downloaded
        finally:
            ix.close()

    def repair(s, bbox: list[float], start: date, end: date) -> int:
        """Unmark stored days whose reflectance is download-failure debris.

        Scans every day with recorded coverage intersecting ``bbox`` over
        ``[start, end]`` and applies the same integrity test as the fill
        gate: a band that is mostly nodata where the day's own fmask has
        valid ground is a failed download, not a data gap. Matching days
        are removed from the ledger so the next :meth:`fill` re-fetches
        them. Returns the number of days unmarked.

        Exists because fills before the integrity gate (or interrupted by
        SIGKILL during a service outage) could mark such days as complete
        — after which nothing would ever re-download them.
        """
        window = grid.tight_window_for_bbox(bbox)
        row0, row1, col0, col1 = window
        fmask_band = s.sentinel2.cloud_mask_band
        bands = [b for b in s.sentinel2.bands if b != fmask_band]

        root = zarr.open_group(s.paths.store, mode='r', use_consolidated=False)
        ix = s._index()
        repaired = 0
        try:
            for day in sorted(ix.scenes_for_range(start, end, s.sentinel2.max_cloud_cover)):
                # only days whose coverage intersects the window
                if not any(max(r0, row0) < min(r1, row1) and max(c0, col0) < min(c1, col1)
                           for r0, r1, c0, c1 in ix.covered_rects(day)):
                    continue
                try:
                    day_group = zarr.open_group(s.paths.store, path=day,
                                                mode='r', use_consolidated=False)
                    fm = day_group[fmask_band][row0:row1, col0:col1]
                except (KeyError, FileNotFoundError, zarr.errors.GroupNotFoundError):
                    continue
                valid = fm != s.sentinel2.fmask_nodata
                n_valid = int(valid.sum())
                if n_valid == 0:
                    continue
                for band in bands:
                    try:
                        arr = day_group[band]
                    except KeyError:
                        # fmask-only day: screened out, legitimately no
                        # reflectance — but only if the screen would still
                        # reject it; a clear day with no reflectance array
                        # at all is failure debris.
                        if _window_wants_reflectance(fm, s.sentinel2):
                            break
                        continue
                    nodata = arr.attrs.get('nodata', arr.fill_value)
                    vals = arr[row0:row1, col0:col1]
                    bad = ((vals == nodata) | np.isnan(vals.astype('float64'))) & valid
                    if bad.sum() / n_valid > _FAILED_NODATA_FRACTION:
                        break
                else:
                    continue
                ix.unmark_day(day)
                repaired += 1
            if repaired:
                print(f'pysentinel2: repair unmarked {repaired} day(s) for re-download')
            return repaired
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

    def _load_days(s, items_by_day: dict[str, list[dict]], bands, window,
                   threads) -> dict[str, Dataset]:
        """Pixels for ``bands`` over ``window`` for every day at once.

        One ``odc.stac.load`` for the whole set of days: the dask graph has
        (days x chunks x bands) tasks, so ``threads`` I/O workers stay busy
        across day boundaries instead of paying each day's round-trip
        latency serially. Returns ``{day: single-timestep Dataset}``; days
        that yield no data are absent from the result.
        """
        import pystac
        import odc.stac
        items = [pystac.Item.from_dict(d)
                 for dicts in items_by_day.values() for d in dicts]
        if not items:
            return {}
        # Scheduler goes to compute() rather than dask.config: config
        # mutation is process-global and callers may be concurrent.
        ds: Dataset = odc.stac.load(
            items,
            bands=bands,
            geobox=grid.geobox_for_window(window),
            groupby=s.sentinel2.groupby,
            chunks={'time': 1, 'x': grid.CHUNK, 'y': grid.CHUNK},
            # One corrupt DEA tile costs a nodata gap, not the whole day.
            fail_on_error=False,
        ).compute(scheduler='threads', num_workers=threads)

        out: dict[str, Dataset] = {}
        for t in range(ds.time.size):
            day = str(ds.time.values[t])[:10]
            slice_ = ds.isel(time=t)
            if day in out:
                # groupby='solar_day' should yield one timestep per day;
                # collapse defensively if odc still splits one.
                out[day] = xr.concat([out[day], slice_], dim='__dup__').max(
                    dim='__dup__', keep_attrs=True)
            else:
                out[day] = slice_
        return out

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


    # -- read -------------------------------------------------------------

    def get_ds(s, bbox: list[float], start: date, end: date, clean: bool = False,
            indices: tuple[str, ...] = (), threads: int = 16,
            bands: tuple[str, ...] | None = None,
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
            bands: Optional subset of stored bands to read. ``None`` reads
                everything. The fill always covers all bands; this trims
                only the *read*, which matters for RAM — an 8-year window
                of all 11 bands is several GB once cleaning promotes to
                float, while e.g. the RGB bands are a fraction of that.
                The fmask band is added automatically when ``clean`` or
                ``indices`` need it.
            **clean_kwargs: Forwarded to :func:`clean_dataset` —
                ``max_cloud_fraction``, ``min_valid_fraction``,
                ``mask_snow``, ``mask_water``, ``buffer_px``.

        Returns:
            xarray.Dataset with dims ``(time, y, x)`` on the fixed grid,
            time = solar days (cloud-filtered per the ``Sentinel2`` config).
        """
        s.fill(bbox, start, end, threads=threads)
        ix = s._index()
        try:
            by_day = ix.scenes_for_range(start, end, s.sentinel2.max_cloud_cover)
        finally:
            ix.close()

        if bands is not None and (clean or indices):
            bands = tuple(dict.fromkeys((*bands, s.sentinel2.cloud_mask_band)))

        # Read and deliver exactly the requested bbox. Storage, fill and
        # dedup stay chunk-aligned (``window``), but the chunk-snapped
        # window carries up to 3.6x the requested pixels as padding, and
        # every downstream cost — this read's RAM, cleaning, SAM,
        # unmixing, videos — pays for it (e.g. SAM 737 s padded vs
        # ~210 s tight). Zarr slices across chunk boundaries natively.
        ds = s._read_window(grid.tight_window_for_bbox(bbox),
                            sorted(by_day), bands=bands)

        if clean or indices:
            ds = clean_dataset(ds, s.sentinel2, **clean_kwargs)
        if indices:
            from pysentinel2.derive import add_indices
            ds = add_indices(ds, indices)
        return ds

    # -- Query adapters (the reproducibility layer speaks Query) ----------

    def fill_query(s, query, threads: int = 16) -> int:
        """:meth:`fill` for a :class:`borevitz_lab.query.Query`."""
        return s.fill(query.bbox, query.start, query.end, threads=threads)

    def get_ds_query(s, query, clean: bool = False,
                  indices: tuple[str, ...] = (), threads: int = 16,
                  **clean_kwargs) -> Dataset:
        """:meth:`get_ds` for a :class:`borevitz_lab.query.Query`."""
        return s.get_ds(query.bbox, query.start, query.end, clean=clean,
                     indices=indices, threads=threads, **clean_kwargs)

    def _read_window(s, window, days: list[str],
                     bands: tuple[str, ...] | None = None) -> Dataset:
        row0, row1, col0, col1 = window
        y, x = grid.coords_for_window(window)
        root = zarr.open_group(s.paths.store, mode='r', use_consolidated=False) if days else None

        band_names = tuple(bands) if bands is not None else s.sentinel2.bands

        def _read_day(day):
            """All requested bands for one day — run on a worker thread.

            Returns ``None`` for a searched day with no scenes written,
            else ``{band: (array, nodata)}``.
            """
            try:
                day_group = zarr.open_group(s.paths.store, path=day,
                                            mode='r', use_consolidated=False)
            except (FileNotFoundError, zarr.errors.GroupNotFoundError):
                return None
            out = {}
            for band in band_names:
                try:
                    arr = day_group[band]
                except KeyError:
                    # A screened day: fmask was stored, reflectance was
                    # never fetched — read as nodata.
                    is_mask = band == s.sentinel2.cloud_mask_band
                    fill = 0 if is_mask else -999
                    out[band] = (np.full(
                        (row1 - row0, col1 - col0), fill,
                        dtype='uint8' if is_mask else 'int16'), fill)
                    continue
                out[band] = (arr[row0:row1, col0:col1],
                             arr.attrs.get('nodata', arr.fill_value))
            return out

        # Chunk decompression releases the GIL, and a multi-year window is
        # thousands of independent little reads — serially they run at a
        # fraction of disk bandwidth (measured ~5 s per year of a 4-chunk
        # window; ~4x faster with the pool).
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=8) as pool:
            per_day = list(pool.map(_read_day, days))

        bands: dict[str, list[np.ndarray]] = {b: [] for b in band_names}
        nodata: dict[str, int] = {}
        kept_days = []
        for day, day_data in zip(days, per_day):
            if day_data is None:
                continue
            kept_days.append(day)
            for band in band_names:
                arr, nd = day_data[band]
                nodata[band] = nd
                bands[band].append(arr)

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
    root = zarr.open_group(cube.paths.store, mode='a', use_consolidated=False)
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
    window = grid.window_for_bbox(_TEST_BBOX)   # generous cover, superset of tight
    ix = cube._index()
    ix.upsert_scenes([(item_id, day, 1.0, {'id': item_id})])
    ix.record_search((-1e9, -1e9, 1e9, 1e9), _TEST_START, _TEST_END)
    ix.mark_rect(day, window)
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
    ds0 = to_float(clean_dataset(_synthetic_window(fmask), buffer_px=0))
    ds3 = to_float(clean_dataset(_synthetic_window(fmask), buffer_px=3))
    neighbour0 = float(ds0['nbart_red'][0, 20, 22])
    neighbour3 = ds3['nbart_red'][0, 20, 22]
    return neighbour0 == 4000.0 and bool(np.isnan(neighbour3))


def test_snow_masked_water_kept_by_default():
    fmask = np.ones((1, 40, 40))
    fmask[0, 0, 0] = defaultsentinel2.fmask_snow
    fmask[0, 0, 1] = defaultsentinel2.fmask_water
    ds = to_float(clean_dataset(_synthetic_window(fmask), buffer_px=0))
    snow_px = ds['nbart_red'][0, 0, 0]
    water_px = ds['nbart_red'][0, 0, 1]
    ds_w = to_float(clean_dataset(_synthetic_window(fmask), buffer_px=0, mask_water=True))
    return (
        bool(np.isnan(snow_px)) and float(water_px) == 4000.0
        and bool(np.isnan(ds_w['nbart_red'][0, 0, 1]))
    )


def test_snow_excluded_from_frame_gate():
    """A cloud-free frame that is 75% snow must survive the gate with its
    snow pixels masked and its snow_fraction reported (analysis/ finding 2)."""
    fmask = np.ones((1, 40, 40))
    fmask[0, :30, :] = defaultsentinel2.fmask_snow
    ds = to_float(clean_dataset(_synthetic_window(fmask), buffer_px=0))
    return (
        ds.time.size == 1
        and float(ds.cloud_fraction[0]) == 0.0
        and abs(float(ds.snow_fraction[0]) - 0.75) < 1e-6
        and bool(np.isnan(ds['nbart_red'][0, 0, 0]))       # snow pixel masked
        and float(ds['nbart_red'][0, 35, 0]) == 4000.0     # clear pixel kept
    )


def test_plan_day_windows():
    """Ring around a covered island -> the 4 frame rects; small covered
    corner -> single bbox; fully covered -> no work; fully missing ->
    the window itself; heavy fragmentation -> bbox fallback."""
    from pysentinel2.cube import _plan_day_windows
    w = (0, 100, 0, 100)
    ring = _plan_day_windows(w, [(30, 60, 40, 80)])         # 12% island
    corner = _plan_day_windows(w, [(0, 10, 0, 10)])          # 1% covered
    covered = _plan_day_windows(w, [(0, 100, 0, 100)])
    empty = _plan_day_windows(w, [])
    frag = _plan_day_windows(w, [(i, i + 2, 10, 90) for i in range(10, 90, 10)])
    ring_area = sum((r1 - r0) * (c1 - c0) for r0, r1, c0, c1 in ring)
    return (len(ring) == 4 and ring_area == 100 * 100 - 30 * 40
            and len(corner) == 1                              # bbox fallback
            and covered == [] and empty == [(0, 100, 0, 100)]
            and len(frag) == 1)                               # fragmentation cap


def test_screen_legacy_recognition():
    """Repair's fmask-only-day recognition: clear window wants
    reflectance, cloudy window does not, off-swath window does not."""
    from pysentinel2.cube import _window_wants_reflectance
    clear = np.ones((64, 64), dtype='uint8')
    cloudy = np.full((64, 64), defaultsentinel2.fmask_cloud, dtype='uint8')
    nodata = np.full((64, 64), defaultsentinel2.fmask_nodata, dtype='uint8')
    return (_window_wants_reflectance(clear, defaultsentinel2)
            and not _window_wants_reflectance(cloudy, defaultsentinel2)
            and _window_wants_reflectance(cloudy, defaultsentinel2, threshold=1.0)
            and not _window_wants_reflectance(nodata, defaultsentinel2))


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
        test_screen_legacy_recognition(),
        test_plan_day_windows(),
        test_query_adapters_match_agnostic_calls(),
    ])


if __name__ == '__main__':
    print(test())
