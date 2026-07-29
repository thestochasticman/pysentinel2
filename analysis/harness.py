"""Parameter-selection harness: difficult-region survey for cleaning/download settings.

Downloads small windows over regions chosen to stress different failure
modes of Sentinel-2 cleaning, with the scene-level cloud filter disabled
(max_cloud_cover=100) so the analysis sees every acquisition:

- alpine_snow      Perisher Valley, NSW — winter snow cover
- tas_west_cloud   near Queenstown, TAS — persistent frontal cloud
- wet_tropics      Tully, QLD — monsoonal convective cloud
- arid_control     near Alice Springs, NT — near-permanent clear sky
- coastal_mixed    Jervis Bay, NSW — open water + swath edges
- cropping_control Grenfell, NSW — the documented example window

Outputs, under analysis/data/:
- frames.csv  one row per (region, solar day): scene-level eo:cloud_cover
  vs window-level fmask statistics (valid/cloud/shadow/snow/water fractions)
- rings.csv   near-cloud reflectance profiles: for partly cloudy frames,
  median band reflectance of clear pixels at distance d (in pixels) from
  the nearest fmask cloud/shadow, normalised by the far-field (>12 px)
  value of the same frame — the empirical basis for choosing buffer_px

Run inside the borevitz_lab environment:  python analysis/harness.py
Fills are incremental; re-runs are served from the machine-wide cube.
"""
import os
import csv
import numpy as np
from datetime import date

from scipy.ndimage import binary_dilation
from pysentinel2.cube import Cube
from pysentinel2.sentinel2 import Sentinel2

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')
os.makedirs(DATA, exist_ok=True)

REGIONS = {
    'alpine_snow':      ([148.395, -36.415, 148.415, -36.395], date(2023, 6, 1), date(2023, 8, 31)),
    'tas_west_cloud':   ([145.545, -42.090, 145.565, -42.070], date(2023, 6, 1), date(2023, 8, 31)),
    'wet_tropics':      ([145.920, -17.940, 145.940, -17.920], date(2024, 1, 1), date(2024, 2, 29)),
    'arid_control':     ([133.860, -23.720, 133.880, -23.700], date(2024, 1, 1), date(2024, 2, 29)),
    'coastal_mixed':    ([150.740, -35.060, 150.760, -35.040], date(2024, 1, 1), date(2024, 2, 29)),
    'cropping_control': ([148.36265, -33.52606, 148.38265, -33.50606], date(2023, 12, 15), date(2024, 2, 5)),
}

RING_MAX = 12          # px; distances analysed are 1..RING_MAX
s2_all = Sentinel2(max_cloud_cover=100.0)   # disable the scene-level filter
cube = Cube(sentinel2=s2_all)
BANDS = [b for b in s2_all.bands if b != s2_all.cloud_mask_band]


def frame_stats(fm: np.ndarray) -> dict:
    n = fm.size
    valid = fm != s2_all.fmask_nodata
    nv = int(valid.sum())
    frac = lambda cls: float((fm == cls).sum()) / nv if nv else 0.0
    contaminated = np.isin(fm, [s2_all.fmask_cloud, s2_all.fmask_shadow]) & valid
    return {
        'valid_frac': nv / n,
        'cloud_frac': float(contaminated.sum()) / nv if nv else 1.0,
        'cloud_only_frac': frac(s2_all.fmask_cloud),
        'shadow_frac': frac(s2_all.fmask_shadow),
        'snow_frac': frac(s2_all.fmask_snow),
        'water_frac': frac(s2_all.fmask_water),
    }


def ring_profile(frame, fm) -> list[dict]:
    """Median clear-pixel reflectance vs distance to nearest cloud/shadow,
    normalised by the same frame's far-field (> RING_MAX px) median."""
    bad = np.isin(fm, [s2_all.fmask_cloud, s2_all.fmask_shadow])
    clear = fm == 1          # fmask class 1 = clear land
    prev = bad
    disks = {}
    for d in range(1, RING_MAX + 1):
        yy, xx = np.ogrid[-d:d + 1, -d:d + 1]
        disks[d] = (yy ** 2 + xx ** 2) <= d ** 2
    far = clear & ~binary_dilation(bad, structure=disks[RING_MAX])
    if far.sum() < 2000:
        return []
    base = {b: float(np.median(frame[b].values[far])) for b in BANDS}
    if any(v <= 0 for v in base.values()):
        return []
    out = []
    for d in range(1, RING_MAX + 1):
        cur = binary_dilation(bad, structure=disks[d])
        ring = clear & cur & ~prev
        prev = cur
        if ring.sum() < 300:
            continue
        row = {'distance_px': d, 'n_pixels': int(ring.sum())}
        for b in BANDS:
            row[b] = float(np.median(frame[b].values[ring])) / base[b]
        out.append(row)
    return out


def min_scene_cloud(ix, day: str) -> float | None:
    rows = ix.db.execute(
        'SELECT cloud_cover FROM scenes WHERE solar_day = ?', (day,)).fetchall()
    vals = [r[0] for r in rows if r[0] is not None]
    return min(vals) if vals else None


frames_rows, rings_rows = [], []
for name, (bbox, start, end) in REGIONS.items():
    print(f'=== {name}: filling...', flush=True)
    n = cube.fill(bbox, start, end)
    print(f'    downloaded {n} cells', flush=True)
    ds = cube.get_ds(bbox, start, end)
    ix = cube._index()
    try:
        for i in range(ds.time.size):
            day = str(ds.time.values[i])[:10]
            fm = ds[s2_all.cloud_mask_band].isel(time=i).values
            st = frame_stats(fm)
            st.update(region=name, solar_day=day, scene_cloud=min_scene_cloud(ix, day))
            frames_rows.append(st)
            if 0.03 <= st['cloud_frac'] <= 0.7 and st['valid_frac'] > 0.5:
                for row in ring_profile(ds.isel(time=i), fm):
                    row.update(region=name, solar_day=day)
                    rings_rows.append(row)
    finally:
        ix.close()
    print(f'    {ds.time.size} frames analysed', flush=True)

with open(f'{DATA}/frames.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['region', 'solar_day', 'scene_cloud', 'valid_frac',
                                      'cloud_frac', 'cloud_only_frac', 'shadow_frac',
                                      'snow_frac', 'water_frac'])
    w.writeheader()
    w.writerows(frames_rows)

with open(f'{DATA}/rings.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['region', 'solar_day', 'distance_px', 'n_pixels'] + BANDS)
    w.writeheader()
    w.writerows(rings_rows)

print(f'wrote {len(frames_rows)} frame rows, {len(rings_rows)} ring rows')
