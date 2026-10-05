"""Derived on-disk locations of the machine-wide Sentinel-2 cube.

The cube is keyed by :class:`troi.Config` (one store per
data root, shared by every request on this machine). Rule of thumb
across the lab's packages: user-settable inputs → Config, derived
locations → Paths. No inheritance — composition only.
"""
from attrs import frozen, field
from troi import Config, config as default_config


@frozen
class Paths:
    """Where the pysentinel2 cube lives for a given Config.

    Attributes:
        config: The :class:`troi.Config` supplying the data root.
        root: Cube directory (``{config.tmp_dir}/sentinel2_cube``). Cross-node
            claims live under ``{root}/claims`` (see :mod:`troi.ledger`).
        store: The sparse Zarr store holding every downloaded pixel.
        index: Marker trees of coverage rects, seen scenes and past
            searches (:mod:`pysentinel2.index`).

    Example:
        ```python
        from pysentinel2.paths import Paths

        paths = Paths()
        paths.store  # '~/Downloads/Troi-Tmp/sentinel2_cube/cube.zarr'
        paths.index  # '~/Downloads/Troi-Tmp/sentinel2_cube/index'
        ```
    """

    config: Config = default_config

    root: str = field(init=False)
    store: str = field(init=False)
    index: str = field(init=False)

    root.default(lambda s: f'{s.config.tmp_dir}/sentinel2_cube')
    store.default(lambda s: f'{s.root}/cube.zarr')
    index.default(lambda s: f'{s.root}/index')


def test_paths_derive_from_config():
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix='pysentinel2_paths_test_')
    cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir)
    paths = Paths(cfg)
    return (
        paths.root == f'{tmpdir}/sentinel2_cube'
        and paths.store == f'{tmpdir}/sentinel2_cube/cube.zarr'
        and paths.index == f'{tmpdir}/sentinel2_cube/index'
    )


def test():
    return test_paths_derive_from_config()


if __name__ == '__main__':
    print(test())
