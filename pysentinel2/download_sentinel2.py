"""Download the Sentinel-2 window for a query — via the machine-wide cube.

Thin compatibility wrapper: the heavy lifting (grid math, chunk-level
dedup, STAC search caching, Zarr writes) lives in
:class:`pysentinel2.cube.Cube`. Kept as a module so the familiar
``download_sentinel2(query)`` entry point survives the storage refactor.
"""

from xarray import Dataset
from borevitz_lab.query import Query
from pysentinel2.sentinel2 import Sentinel2, defaultsentinel2


def download_sentinel2(
    query: Query,
    threads_per_worker: int = 8,
    sentinel2: Sentinel2 = defaultsentinel2,
) -> Dataset:
    """Return the raw Sentinel-2 cube window (incl. fmask) for ``query``.

    Fills only the (day x chunk) cells of the shared cube that no previous
    query has populated — repeat or overlapping queries re-download nothing.

    Args:
        query: The :class:`borevitz_lab.query.Query` describing region + range.
        threads_per_worker: I/O concurrency for any downloads triggered.
        sentinel2: STAC/band/cloud configuration; defaults to the DEA config.

    Returns:
        xarray.Dataset with dims ``(time, y, x)`` on the fixed grid.
    """
    from pysentinel2.cube import Cube
    cube = Cube(config=query.config, sentinel2=sentinel2)
    return cube.get_ds_query(query, threads=threads_per_worker)


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
        from borevitz_lab.config import Config
        tmpdir = tempfile.mkdtemp(prefix='pysentinel2_s2_test_')
        _test_cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir)
    return _test_cfg


def test_download_returns_data():
    """A live download returns a non-empty window with the requested bands."""
    q = Query(
        bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
        stub='s2_live', config=_shared_test_cfg(),
    )
    ds = download_sentinel2(q)
    return ds.time.size > 0 and 'nbart_red' in ds.data_vars and 'oa_fmask' in ds.data_vars


def test_repeat_query_downloads_nothing():
    """Second identical query → fill() reports 0 cells downloaded."""
    from pysentinel2.cube import Cube
    q = Query(
        bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
        stub='s2_repeat', config=_shared_test_cfg(),
    )
    download_sentinel2(q)
    return Cube(config=q.config).fill_query(q) == 0


def test_overlapping_query_downloads_only_new_cells():
    """A bbox shifted ~1 km east reuses the shared chunks — the number of
    newly downloaded cells is strictly less than its total cell count."""
    from pysentinel2.cube import Cube
    from pysentinel2 import grid
    cfg = _shared_test_cfg()
    q1 = Query(bbox=_TEST_BBOX, start=_TEST_START, end=_TEST_END,
               stub='s2_overlap_a', config=cfg)
    download_sentinel2(q1)

    shifted = [_TEST_BBOX[0] + 0.01, _TEST_BBOX[1], _TEST_BBOX[2] + 0.01, _TEST_BBOX[3]]
    q2 = Query(bbox=shifted, start=_TEST_START, end=_TEST_END,
               stub='s2_overlap_b', config=cfg)
    cube = Cube(config=cfg)
    n_days = len(cube._index().scenes_for_range(q2.start, q2.end, cube.sentinel2.max_cloud_cover))
    total_cells = len(grid.chunks_in_window(grid.window_for_bbox(q2.bbox))) * max(n_days, 1)
    downloaded = cube.fill_query(q2)
    return 0 <= downloaded < total_cells


def test():
    return all([
        test_internet(None),
        test_download_returns_data(),
        test_repeat_query_downloads_nothing(),
        test_overlapping_query_downloads_only_new_cells(),
    ])


if __name__ == '__main__':
    print(test())
