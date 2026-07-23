"""Spectral indices, computed on read — never stored.

Vegetation and tillage indices over the cube's Sentinel-2 reflectance,
in the same on-read-derivative spirit as cloud masking: request them
per call, pay array math, store nothing. Formulas are written against
DEA's ARD band names (``nbart_*``); DN values are scaled by 1/10000 to
reflectance, and both DN 0 and the band's nodata value are treated as
missing.
"""
import numpy as np
from xarray import Dataset


def _band(ds: Dataset, name: str) -> np.ndarray:
    """Band as float32 reflectance ``(time, y, x)``; 0 and nodata → NaN."""
    da = ds[name]
    b = da.values.astype(np.float32)
    b[b == 0] = np.nan
    nodata = da.attrs.get('nodata')
    if nodata is not None:
        b[b == float(nodata)] = np.nan
    b /= 10000.0
    return b


def _normalised_diff(a, b):
    with np.errstate(invalid='ignore', divide='ignore'):
        nd = (a - b) / (a + b)
    nd[~np.isfinite(nd)] = np.nan
    return nd


def ndvi(ds: Dataset) -> np.ndarray:
    """Normalised Difference Vegetation Index ``(NIR - Red) / (NIR + Red)``."""
    return _normalised_diff(_band(ds, 'nbart_nir_1'), _band(ds, 'nbart_red'))


def cfi(ds: Dataset) -> np.ndarray:
    """Crop Foliage Index ``NDVI * (Red + 2*Green - Blue)``."""
    red = _band(ds, 'nbart_red')
    green = _band(ds, 'nbart_green')
    blue = _band(ds, 'nbart_blue')
    return ndvi(ds) * (red + green + green - blue)


def nirv(ds: Dataset) -> np.ndarray:
    """Near-Infrared Reflectance of Vegetation ``NDVI * NIR``."""
    return ndvi(ds) * _band(ds, 'nbart_nir_1')


def ndti(ds: Dataset) -> np.ndarray:
    """Normalised Difference Tillage Index ``(SWIR2 - SWIR3) / (SWIR2 + SWIR3)``."""
    return _normalised_diff(_band(ds, 'nbart_swir_2'), _band(ds, 'nbart_swir_3'))


def cai(ds: Dataset) -> np.ndarray:
    """Cellulose Absorption Index ``0.5 * (SWIR2 + SWIR3) - NIR``."""
    return 0.5 * (_band(ds, 'nbart_swir_2') + _band(ds, 'nbart_swir_3')) - _band(ds, 'nbart_nir_1')


INDICES = {'NDVI': ndvi, 'CFI': cfi, 'NIRv': nirv, 'NDTI': ndti, 'CAI': cai}


def add_indices(ds: Dataset, indices) -> Dataset:
    """Return ``ds`` with one float32 ``(time, y, x)`` variable per index."""
    unknown = set(indices) - set(INDICES)
    if unknown:
        raise ValueError(f'Unknown index/indices: {sorted(unknown)} — pick from {sorted(INDICES)}')
    for name in indices:
        ds[name] = (('time', 'y', 'x'), INDICES[name](ds))
    return ds


def _synthetic(nir=6000, red=2000, green=1500, blue=1000, swir2=3000, swir3=1000):
    from xarray import DataArray
    shape = (1, 2, 2)
    return Dataset({
        band: DataArray(np.full(shape, dn, dtype='int16'), dims=('time', 'y', 'x'),
                        attrs={'nodata': -999})
        for band, dn in {
            'nbart_nir_1': nir, 'nbart_red': red, 'nbart_green': green,
            'nbart_blue': blue, 'nbart_swir_2': swir2, 'nbart_swir_3': swir3,
        }.items()
    })


def test_ndvi_known_value():
    # refl nir=.6 red=.2 -> (0.6-0.2)/(0.6+0.2) = 0.5
    return abs(float(ndvi(_synthetic())[0, 0, 0]) - 0.5) < 1e-6


def test_ndti_known_value():
    # (0.3-0.1)/(0.3+0.1) = 0.5
    return abs(float(ndti(_synthetic())[0, 0, 0]) - 0.5) < 1e-6


def test_cai_known_value():
    # 0.5*(0.3+0.1) - 0.6 = -0.4
    return abs(float(cai(_synthetic())[0, 0, 0]) - (-0.4)) < 1e-6


def test_cfi_known_value():
    # 0.5 * (0.2 + 2*0.15 - 0.1) = 0.2
    return abs(float(cfi(_synthetic())[0, 0, 0]) - 0.2) < 1e-6


def test_nodata_and_zero_are_masked():
    ds = _synthetic()
    ds['nbart_red'].values[0, 0, 0] = 0
    ds['nbart_red'].values[0, 0, 1] = -999
    out = ndvi(ds)
    return np.isnan(out[0, 0, 0]) and np.isnan(out[0, 0, 1]) and np.isfinite(out[0, 1, 1])


def test_add_indices_rejects_unknown():
    try:
        add_indices(_synthetic(), ('NDWI',))
    except ValueError:
        return True
    return False


def test():
    return all([
        test_ndvi_known_value(),
        test_ndti_known_value(),
        test_cai_known_value(),
        test_cfi_known_value(),
        test_nodata_and_zero_are_masked(),
        test_add_indices_rejects_unknown(),
    ])


if __name__ == '__main__':
    print(test())
