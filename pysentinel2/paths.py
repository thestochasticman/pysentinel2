"""Derived on-disk locations of the machine-wide Sentinel-2 cube.

The cube is keyed by :class:`borevitz_lab.config.Config` (one store per
data root, shared by every query on this machine), not per-Query — that's
the whole point: all queries read and fill the same store. Rule of thumb
across the lab's packages: user-settable inputs → Config, derived
locations → Paths. No inheritance — composition only.
"""
from attrs import frozen, field
from borevitz_lab.config import Config, config as default_config


@frozen
class Paths:
    """Where the pysentinel2 cube lives for a given Config.

    Attributes:
        config: The :class:`borevitz_lab.config.Config` supplying the data root.
        root: Cube directory (``{config.tmp_dir}/sentinel2_cube``).
        store: The sparse Zarr store holding every downloaded pixel.
        index_db: SQLite index of populated chunks / seen scenes / searches.

    Example:
        ```python
        from pysentinel2.paths import Paths

        paths = Paths()
        paths.store     # '~/Downloads/BorevitzLab-Tmp/sentinel2_cube/cube.zarr'
        paths.index_db  # '~/Downloads/BorevitzLab-Tmp/sentinel2_cube/index.db'
        ```
    """

    config: Config = default_config

    root: str = field(init=False)
    store: str = field(init=False)
    index_db: str = field(init=False)

    root.default(lambda s: f'{s.config.tmp_dir}/sentinel2_cube')
    store.default(lambda s: f'{s.root}/cube.zarr')
    index_db.default(lambda s: f'{s.root}/index.db')


def test_paths_derive_from_config():
    import tempfile
    tmpdir = tempfile.mkdtemp(prefix='pysentinel2_paths_test_')
    cfg = Config(out_dir=tmpdir, tmp_dir=tmpdir)
    paths = Paths(cfg)
    return (
        paths.root == f'{tmpdir}/sentinel2_cube'
        and paths.store == f'{tmpdir}/sentinel2_cube/cube.zarr'
        and paths.index_db == f'{tmpdir}/sentinel2_cube/index.db'
    )


def test():
    return test_paths_derive_from_config()


if __name__ == '__main__':
    print(test())
