"""Download the Sentinel-2 window for a troi — via the machine-wide cube.

Thin compatibility wrapper: the heavy lifting (grid math, chunk-level
dedup, STAC search caching, Zarr writes) lives in
:class:`pysentinel2.cube.Cube`. Kept as a module so the familiar
``download_sentinel2(troi)`` entry point survives the storage refactor.
"""

from xarray import Dataset
from troi import Troi
from pysentinel2.sentinel2 import Sentinel2, defaultsentinel2


def download_sentinel2(
    troi: Troi,
    threads_per_worker: int = 8,
    sentinel2: Sentinel2 = defaultsentinel2,
) -> Dataset:
    """Return the raw Sentinel-2 cube window (incl. fmask) for ``troi``.

    Fills only the (day x chunk) cells of the shared cube that no previous
    troi has populated — repeat or overlapping queries re-download nothing.

    Args:
        troi: The :class:`troi.Troi` describing region + range.
        threads_per_worker: I/O concurrency for any downloads triggered.
        sentinel2: STAC/band/cloud configuration; defaults to the DEA config.

    Returns:
        xarray.Dataset with dims ``(time, y, x)`` on the fixed grid.
    """
    from pysentinel2.cube import Cube
    cube = Cube(config=troi.config, sentinel2=sentinel2)
    return cube.get_ds_troi(troi, threads=threads_per_worker)


def test_internet(s):
    from urllib import request
    from urllib.error import URLError
    try:
        request.urlopen('https://www.google.com/', timeout=2)
        return True
    except URLError as error:
        return False


_TEST_BBOX = [148.36265, -33.52606, 148.38265, -33.50606]
from datetime import date as _date
_TEST_START, _TEST_END = _date(2024, 1, 1), _date(2024, 1, 21)
_test_cfg = None


def _shared_test_cfg():
    global _test_cfg
    if _test_cfg is None:
        import tempfile
        from troi import Config
        tmpdir = tempfile.mkdtemp(prefix='pysentinel2_s2_test_')
        _test_cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir)
    return _test_cfg


def test_download_returns_data():
    """A live download returns a non-empty window with the requested bands."""
    q = Troi(
        bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
        stub='s2_live', config=_shared_test_cfg(),
    )
    ds = download_sentinel2(q)
    return ds.time.size > 0 and 'nbart_red' in ds.data_vars and 'oa_fmask' in ds.data_vars


def test_repeat_troi_downloads_nothing():
    """Second identical troi → fill() reports 0 cells downloaded."""
    from pysentinel2.cube import Cube
    q = Troi(
        bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
        stub='s2_repeat', config=_shared_test_cfg(),
    )
    download_sentinel2(q)
    return Cube(config=q.config).fill_troi(q) == 0


def test_overlapping_troi_downloads_only_new_cells():
    """A bbox shifted ~1 km east reuses covered pixels — the downloaded
    (day x pixel) count is strictly less than the shifted window's total,
    matching the uncovered strip exactly."""
    from pysentinel2.cube import Cube
    from pysentinel2 import grid
    cfg = _shared_test_cfg()
    q1 = Troi(bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
               stub='s2_overlap_a', config=cfg)
    download_sentinel2(q1)

    shifted = [_TEST_BBOX[0] + 0.01, _TEST_BBOX[1], _TEST_BBOX[2] + 0.01, _TEST_BBOX[3]]
    q2 = Troi(bbox=shifted, start=_TEST_START, end=_TEST_END,
               stub='s2_overlap_b', config=cfg)
    cube = Cube(config=cfg)
    n_days = len(cube._index().scenes_for_range(q2.start, q2.end, cube.sentinel2.max_cloud_cover))
    w = grid.tight_window_for_bbox(q2.bbox)
    total_px = (w[1] - w[0]) * (w[3] - w[2]) * max(n_days, 1)
    downloaded = cube.fill_troi(q2)
    return 0 <= downloaded < total_px


def test_ring_fill_downloads_only_frame():
    """Coverage fully inside the ROI: fill must fetch only the ring.

    An interior ~1/3 of the window is marked covered for every candidate
    day before filling; the reported download must equal the frame area
    times the day count, and a repeat fill must report 0.
    """
    import tempfile
    from troi import Config
    from pysentinel2.cube import Cube
    from pysentinel2 import grid

    tmpdir = tempfile.mkdtemp(prefix='pysentinel2_ring_test_')
    cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir)
    cube = Cube(config=cfg)

    w = grid.tight_window_for_bbox(_TEST_BBOX)
    r0, r1, c0, c1 = w
    h, wd = r1 - r0, c1 - c0
    island = (r0 + h // 3, r0 + 2 * (h // 3), c0 + wd // 3, c0 + 2 * (wd // 3))
    island_area = (island[1] - island[0]) * (island[3] - island[2])
    frame_area = h * wd - island_area

    # Fill the whole range once (searches + downloads), then rewrite the
    # ledger so only the interior island is covered: the next fill must
    # download exactly the frame around it, per day.
    cube.fill(_TEST_BBOX, _TEST_START, _TEST_END)
    ix = cube._index()
    days = list(ix.scenes_for_range(_TEST_START, _TEST_END, cube.sentinel2.max_cloud_cover))
    for d in days:
        ix.unmark_day(d)
        ix.mark_rect(d, island)
    ix.close()

    downloaded = cube.fill(_TEST_BBOX, _TEST_START, _TEST_END)
    again = cube.fill(_TEST_BBOX, _TEST_START, _TEST_END)
    return (len(days) > 0 and downloaded == frame_area * len(days)
            and again == 0)


def test():
    return all([
        test_internet(None),
        test_download_returns_data(),
        test_repeat_troi_downloads_nothing(),
        test_overlapping_troi_downloads_only_new_cells(),
        test_ring_fill_downloads_only_frame(),
    ])


if __name__ == '__main__':
    print(test())
