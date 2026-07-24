"""Cloud-masked Sentinel-2 window for a query — computed on read.

Thin compatibility wrapper over :func:`pysentinel2.cube.Cube.get_ds` with
``clean=True``. No "clean" copy is ever persisted: masking a window is
cheap, so the clean cube is a view of the raw store rather than a
second store.
"""

from xarray import Dataset
from borevitz_lab.query import Query
from pysentinel2.sentinel2 import Sentinel2, defaultsentinel2


def clean_sentinel2(
    query: Query,
    ds_sentinel2: Dataset | None = None,
    sentinel2: Sentinel2 = defaultsentinel2,
    **clean_kwargs,
) -> Dataset:
    """Produce a cloud-masked, frame-filtered Sentinel-2 window.

    Args:
        query: The :class:`borevitz_lab.query.Query`. Pixels come from the
            shared cube (downloading only what's missing).
        ds_sentinel2: Optional in-memory raw dataset (must still include
            the fmask band); if given, it is masked directly and the cube
            is not touched.
        sentinel2: Config supplying the fmask band and class codes.
            Defaults to the bundled DEA config.
        **clean_kwargs: Forwarded to
            :func:`pysentinel2.cube.clean_dataset` —
            ``max_cloud_fraction``, ``min_valid_fraction``, ``mask_snow``,
            ``mask_water``, ``buffer_px``.

    Returns:
        xarray.Dataset: The cleaned window, fmask band removed, only the
        retained timesteps, with per-frame ``cloud_fraction`` /
        ``valid_fraction`` coordinates.
    """
    from pysentinel2.cube import Cube, clean_dataset
    if ds_sentinel2 is not None:
        return clean_dataset(ds_sentinel2, sentinel2, **clean_kwargs)
    cube = Cube(config=query.config, sentinel2=sentinel2)
    return cube.get_ds_query(query, clean=True, **clean_kwargs)


def test_clean_drops_fmask_band():
    """The fmask band must not appear in the cleaned dataset (live)."""
    import tempfile
    from datetime import date
    from borevitz_lab.config import Config

    tmpdir = tempfile.mkdtemp(prefix='pysentinel2_clean_test_')
    cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir)
    q = Query(
        bbox=[148.36265, -33.52606, 148.38265, -33.50606],
        start=date(2024, 1, 1), end=date(2024, 1, 21),
        stub='clean_no_fmask', config=cfg,
    )
    ds = clean_sentinel2(q, max_cloud_fraction=0.7)
    return defaultsentinel2.cloud_mask_band not in ds.data_vars and ds.time.size > 0


def test():
    from pysentinel2.download_sentinel2 import test_internet
    return all([
        test_internet(None),
        test_clean_drops_fmask_band(),
    ])


if __name__ == '__main__':
    print(test())
