from attrs import frozen, Factory as F
from typing_extensions import Self

@frozen
class Sentinel2:
    stub: str = 'ds2'
    stac_url: str = 'https://explorer.dea.ga.gov.au/stac'
    collections: tuple[str, ...] = (
        'ga_s2am_ard_3',
        'ga_s2bm_ard_3',
    )
    bands: tuple = (
        'oa_fmask',
        'nbart_blue',
        'nbart_green',
        'nbart_red',
        'nbart_red_edge_1',
        'nbart_red_edge_2',
        'nbart_red_edge_3',
        'nbart_nir_1',
        'nbart_nir_2',
        'nbart_swir_2',
        'nbart_swir_3',
    )

    cloud_mask_band: str = 'oa_fmask'
    # Applied per-scene when selecting days from the cube's index (NOT at
    # STAC search time — all scenes are recorded, so relaxing the threshold
    # later needs no re-search, just fills the newly-eligible days).
    # Granule-level eo:cloud_cover is a weak proxy for window-level
    # cloudiness (analysis/: 30 excluded 18% of usable frames), so this
    # stays loose; read-time cleaning is the real cloud gate.
    max_cloud_cover: float = 80.0
    # Legacy of the removed fmask-first download screen (it skipped 8.8%
    # of days' reflectance for an extra request round on every day).
    # Still used by Cube.repair to recognise fmask-only days written by
    # old fills as legitimately screened rather than failure debris.
    screen_cloud_fraction: float = 0.9
    fmask_nodata: int = 0
    fmask_cloud: int = 2
    fmask_shadow: int = 3
    fmask_snow: int = 4
    fmask_water: int = 5
    crs: str = 'EPSG:6933'
    resolution: int = 10
    groupby: str = 'solar_day'

    def __post_init__(s: Self):
        object.__setattr__(s, 'bands', tuple(sorted(list(s.bands))))
        object.__setattr__(s, 'collections', tuple(sorted(list(s.collections))))

defaultsentinel2 = Sentinel2()