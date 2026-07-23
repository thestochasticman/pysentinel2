# Light-weight exports only: pysentinel2.cube (and the download/clean
# wrappers) pull in the heavy geospatial stack (odc.stac, rioxarray, zarr),
# so those stay behind explicit submodule imports.
from pysentinel2.paths import Paths
from pysentinel2.sentinel2 import Sentinel2, defaultsentinel2

__all__ = [
    'Paths',
    'Sentinel2',
    'defaultsentinel2',
]
