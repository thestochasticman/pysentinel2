"""Generate the documentation figures from the real machine-wide cube.

Pixel data comes from the cube at the default data root; anything missing
for the example window is downloaded once from DEA on the first run, and
every later run is served entirely from cache. Outputs land in docs/images/.

Run with the troi conda environment:

    python docs/make_figures.py
"""
import os
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle
from matplotlib.lines import Line2D
from datetime import date
from scipy.ndimage import binary_dilation

from pysentinel2.cube import Cube, clean_dataset
from pysentinel2 import grid

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'images')
os.makedirs(OUT, exist_ok=True)

# ---- palette (validated reference set) ------------------------------------
SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#898781'
GRID = '#e1e0d9'
BASE = '#c3c2b7'
BLUE = '#2a78d6'
ORANGE = '#eb6834'
AQUA = '#1baf7a'
YELLOW = '#eda100'
MAGENTA = '#e87ba4'
GREEN = '#008300'
VIOLET = '#4a3aa7'
BLUE_LIGHT = '#cde2fb'

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 9,
    'text.color': INK,
    'axes.edgecolor': BASE,
    'axes.labelcolor': INK2,
    'xtick.color': MUTED,
    'ytick.color': MUTED,
    'figure.facecolor': SURFACE,
    'axes.facecolor': SURFACE,
    'savefig.facecolor': SURFACE,
    'savefig.dpi': 150,
})

BBOX = [148.36265, -33.52606, 148.38265, -33.50606]
START, END = date(2023, 12, 15), date(2024, 2, 5)

cube = Cube()
s2 = cube.sentinel2
ds = cube.get_ds(BBOX, START, END)          # raw, fully cached
print('raw frames:', ds.time.size, flush=True)

# per-frame statistics over ALL frames (gates disabled)
stats = clean_dataset(ds, s2, max_cloud_fraction=1.0, min_valid_fraction=0.0)
cloud_frac = stats.cloud_fraction.values
valid_frac = stats.valid_fraction.values
days = [str(t)[:10] for t in ds.time.values]
print(list(zip(days, np.round(cloud_frac, 3), np.round(valid_frac, 3))), flush=True)


def rgb(frame_ds):
    """(y, x, 3) true-colour array, reflectance stretched 0–0.25."""
    chans = []
    for b in ('nbart_red', 'nbart_green', 'nbart_blue'):
        v = frame_ds[b].values.astype('float32')
        v[v == -999] = np.nan
        chans.append(v / 10000.0)
    img = np.dstack(chans)
    return np.clip(img / 0.25, 0, 1)


def style_map_ax(ax):
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_color(GRID)
        sp.set_linewidth(0.8)


# ---- 1. every sensed frame in the cube window -----------------------------
sensed = [i for i in range(ds.time.size) if valid_frac[i] > 0]
n_off = ds.time.size - len(sensed)
ncols = 4
nrows = int(np.ceil((len(sensed) + 1) / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(9.6, 2.55 * nrows))
for k, ax in enumerate(axes.flat):
    if k >= len(sensed):
        ax.axis('off')
        if k == len(sensed):
            ax.text(0.5, 0.5, f'+ {n_off} stored days\nentirely off-swath\n(valid = 0)',
                    transform=ax.transAxes, fontsize=9, color=INK2,
                    ha='center', va='center')
        continue
    i = sensed[k]
    img = rgb(ds.isel(time=i))
    img_show = np.where(np.isnan(img), 0.93, img)
    ax.imshow(img_show)
    style_map_ax(ax)
    ax.set_title(f'{days[i]}', fontsize=9, color=INK, pad=4)
    ax.text(0.02, 0.03, f'cloud {cloud_frac[i]:.0%} · valid {valid_frac[i]:.0%}',
            transform=ax.transAxes, fontsize=7.2, color=INK2,
            bbox=dict(facecolor='white', alpha=0.85, edgecolor='none', pad=1.6))
fig.suptitle('Every sensed solar day stored for the example window (true colour, 512 × 512 px @ 10 m)',
             fontsize=11, color=INK, y=1.0)
fig.tight_layout(rect=(0, 0, 1, 0.985))
fig.savefig(f'{OUT}/cube_frames_rgb.png', bbox_inches='tight')
plt.close(fig)
print('fig 1 done', flush=True)

# ---- 2. cleaning pipeline on the cloudiest usable frame -------------------
cand = np.where(valid_frac > 0.5, cloud_frac, -1)
ci = int(np.argmax(cand))
frame = ds.isel(time=ci)
fm = frame[s2.cloud_mask_band].values

valid = fm != s2.fmask_nodata
bad_raw = np.isin(fm, [s2.fmask_cloud, s2.fmask_shadow])
r = 5
yy, xx = np.ogrid[-r:r + 1, -r:r + 1]
disk = (yy ** 2 + xx ** 2) <= r ** 2
bad_dil = binary_dilation(bad_raw, structure=disk) & valid
bad_dil |= (fm == s2.fmask_snow) & valid          # snow masked, undilated

class_colors = {0: GRID, 1: AQUA, 2: YELLOW, 3: VIOLET, 4: MAGENTA, 5: BLUE}
class_names = {0: 'nodata (0)', 1: 'clear (1)', 2: 'cloud (2)',
               3: 'shadow (3)', 4: 'snow (4)', 5: 'water (5)'}
cmap = ListedColormap([class_colors[k] for k in range(6)])
norm = BoundaryNorm(np.arange(-0.5, 6), cmap.N)

img = rgb(frame)
img_show = np.where(np.isnan(img), 0.93, img)
img_clean = img.copy()
img_clean[bad_dil | ~valid] = np.nan
img_clean_show = np.where(np.isnan(img_clean), 0.93, img_clean)

fig, axes = plt.subplots(1, 4, figsize=(12.6, 3.6))
axes[0].imshow(img_show)
axes[0].set_title('a  raw true colour', loc='left', fontsize=9.5, color=INK)
axes[1].imshow(fm, cmap=cmap, norm=norm, interpolation='nearest')
axes[1].set_title('b  fmask classification', loc='left', fontsize=9.5, color=INK)
axes[2].imshow(img_show)
overlay = np.zeros((*bad_dil.shape, 4))
halo = bad_dil & ~bad_raw
overlay[bad_raw] = matplotlib.colors.to_rgba(YELLOW, 0.75)
overlay[halo] = matplotlib.colors.to_rgba(ORANGE, 0.75)
axes[2].imshow(overlay, interpolation='nearest')
axes[2].set_title('c  contaminated mask, dilated 5 px', loc='left', fontsize=9.5, color=INK)
axes[3].imshow(img_clean_show)
axes[3].set_title('d  cleaned observation', loc='left', fontsize=9.5, color=INK)
for ax in axes:
    style_map_ax(ax)

legend1 = [Patch(facecolor=class_colors[k], label=class_names[k]) for k in range(6)]
legend2 = [Patch(facecolor=YELLOW, label='fmask cloud/shadow'),
           Patch(facecolor=ORANGE, label='+ 5 px dilation halo'),
           Patch(facecolor='#ededec', label='masked → NaN')]
fig.legend(handles=legend1 + legend2, loc='lower center', ncol=5, frameon=False,
           fontsize=8, bbox_to_anchor=(0.5, -0.06))
fig.suptitle(f'On-read cleaning of the {days[ci]} frame '
             f'(cloud fraction {cloud_frac[ci]:.0%} of valid pixels)',
             fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0.02, 1, 0.97))
fig.savefig(f'{OUT}/cleaning_pipeline.png', bbox_inches='tight')
plt.close(fig)
print('fig 2 done', flush=True)

# ---- 3. frame gates -------------------------------------------------------
MAXC, MINV = 0.5, 0.2
keep = (valid_frac >= MINV) & (cloud_frac <= MAXC)
fig, ax = plt.subplots(figsize=(6.6, 4.4))
ax.axvspan(MINV, 1.02, ymin=0, ymax=MAXC / 1.05, color=BLUE_LIGHT, alpha=0.45, zorder=0)
ax.axhline(MAXC, color=BASE, lw=1, ls=(0, (4, 3)))
ax.axvline(MINV, color=BASE, lw=1, ls=(0, (4, 3)))
ax.scatter(valid_frac[keep], cloud_frac[keep], s=52, color=BLUE, zorder=3,
           label=f'kept ({keep.sum()})')
ax.scatter(valid_frac[~keep], cloud_frac[~keep], s=52, color=MUTED, zorder=3,
           marker='X', label=f'dropped ({(~keep).sum()})')
n_offswath = int(((valid_frac == 0)).sum())
ax.annotate(f'{n_offswath} off-swath days\n(valid = 0 → dropped)', (0.0, 1.0),
            textcoords='offset points', xytext=(14, -26), fontsize=7.8,
            color=INK2, ha='left', va='top')
n_stack = int((keep & (cloud_frac < 0.02)).sum())
ax.annotate(f'{n_stack} clear frames', (1.0, 0.0),
            textcoords='offset points', xytext=(-10, 8), fontsize=7.8,
            color=INK2, ha='right')
for i in range(ds.time.size):
    if valid_frac[i] > 0 and cloud_frac[i] > 0.1:
        ax.annotate(days[i], (valid_frac[i], cloud_frac[i]),
                    textcoords='offset points', xytext=(-10, -3),
                    fontsize=7.4, color=INK2, ha='right')
ax.text(0.98, MAXC + 0.02, 'max_cloud_fraction = 0.5', fontsize=8, color=INK2,
        va='bottom', ha='right')
ax.text(MINV + 0.012, 0.06, 'min_valid_fraction = 0.2', fontsize=8, color=INK2,
        rotation=90, va='bottom', ha='left')
ax.set_xlabel('valid fraction  (valid ÷ all window pixels)')
ax.set_ylabel('cloud fraction  (contaminated ÷ valid pixels)')
ax.set_xlim(0, 1.02)
ax.set_ylim(-0.03, 1.05)
ax.set_title('Frame gating: the two statistics are independent', loc='left',
             fontsize=11, color=INK, pad=10)
ax.legend(frameon=False, loc='upper center', fontsize=8.5, bbox_to_anchor=(0.55, 0.98))
ax.grid(color=GRID, lw=0.6)
ax.set_axisbelow(True)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(f'{OUT}/frame_gates.png', bbox_inches='tight')
plt.close(fig)
print('fig 3 done', flush=True)

# ---- 4. spectral index maps on the clearest frame -------------------------
best = int(np.argmin(np.where(valid_frac > 0.95, cloud_frac, 2)))
ds_idx = cube.get_ds(BBOX, START, END, indices=('NDVI', 'NIRv', 'NDTI', 'CAI'))
di = [str(t)[:10] for t in ds_idx.time.values].index(days[best])
frame_idx = ds_idx.isel(time=di)

green_ramp = LinearSegmentedColormap.from_list('g', ['#f1f7ee', '#bcd9b0', '#6fb35f', GREEN, '#014b01'])
blue_ramp = LinearSegmentedColormap.from_list('b', ['#eef3fb', BLUE_LIGHT, '#6da7ec', BLUE, '#0d366b'])
orange_ramp = LinearSegmentedColormap.from_list('o', ['#fdf1ec', '#f6c4ab', '#f0926a', ORANGE, '#8f2f0d'])

panels = [
    ('NDVI', green_ramp, (0, 0.9), 'vegetation vigour'),
    ('NIRv', green_ramp, (0, 0.35), 'photosynthetic capacity'),
    ('NDTI', blue_ramp, (0, 0.4), 'residue / tillage'),
    ('CAI', orange_ramp, (-0.4, 0.1), 'cellulose absorption'),
]
fig, axes = plt.subplots(1, 4, figsize=(12.6, 3.5))
for ax, (name, cm, (v0, v1), sub) in zip(axes, panels):
    v = frame_idx[name].values
    cm = cm.copy()
    cm.set_bad('#ededec')
    im = ax.imshow(v, cmap=cm, vmin=v0, vmax=v1, interpolation='nearest')
    style_map_ax(ax)
    ax.set_title(name, loc='left', fontsize=10, color=INK)
    ax.text(1.0, 1.035, sub, transform=ax.transAxes, fontsize=7.6, color=MUTED,
            ha='right', va='bottom')
    cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=7, color=BASE)
    cb.outline.set_edgecolor(GRID)
fig.suptitle(f'On-read spectral indices, {days[best]} (cloud-masked reflectance; nothing stored)',
             fontsize=11, color=INK, y=1.06)
fig.tight_layout()
fig.savefig(f'{OUT}/indices_maps.png', bbox_inches='tight')
plt.close(fig)
print('fig 4 done', flush=True)

# ---- 5. NDVI time series over the season ----------------------------------
ndvi = ds_idx['NDVI'].values
t = [str(v)[:10] for v in ds_idx.time.values]
med = np.nanmedian(ndvi, axis=(1, 2))
q25 = np.nanpercentile(ndvi, 25, axis=(1, 2))
q75 = np.nanpercentile(ndvi, 75, axis=(1, 2))
xs = ds_idx.time.values.astype('datetime64[D]').astype('int64')

fig, ax = plt.subplots(figsize=(7.6, 3.6))
ax.fill_between(xs, q25, q75, color=BLUE_LIGHT, alpha=0.8, lw=0, label='interquartile range')
ax.plot(xs, med, color=BLUE, lw=2, marker='o', ms=5, label='spatial median')
ax.set_ylabel('NDVI')
ax.set_ylim(0, max(0.5, np.nanmax(q75) * 1.15))
ticklab = [t[i][5:] for i in range(len(t))]
ax.set_xticks(xs)
ax.set_xticklabels(ticklab, rotation=0, fontsize=7.4)
ax.set_xlabel('solar day (2023-12 → 2024-02)')
ax.set_title('Window-median NDVI across the cleaned frames', loc='left',
             fontsize=11, color=INK, pad=10)
ax.grid(color=GRID, lw=0.6)
ax.set_axisbelow(True)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
ax.legend(frameon=False, fontsize=8.5, loc='upper right')
fig.tight_layout()
fig.savefig(f'{OUT}/ndvi_timeseries.png', bbox_inches='tight')
plt.close(fig)
print('fig 5 done', flush=True)

# ---- 6. grid & dedup schematic --------------------------------------------
bbox_a = BBOX
bbox_b = [BBOX[0] + 0.022, BBOX[1] - 0.012, BBOX[2] + 0.022, BBOX[3] - 0.012]
wa, wb = grid.window_for_bbox(bbox_a), grid.window_for_bbox(bbox_b)
ca, cb = set(grid.chunks_in_window(wa)), set(grid.chunks_in_window(wb))
shared = ca & cb
print('chunks A only / B only / shared:', len(ca - cb), len(cb - ca), len(shared), flush=True)

allc = ca | cb
cys = [c[0] for c in allc]
cxs = [c[1] for c in allc]
XREF = grid.X0 + min(cxs) * grid.CHUNK_M          # local origin, km axes
YREF = grid.Y_TOP - (max(cys) + 1) * grid.CHUNK_M


def to_km(x, y):
    return (x - XREF) / 1000.0, (y - YREF) / 1000.0


def bbox_km(b):
    x0, y0, x1, y1 = grid.bbox_to_6933(b)
    (x0, y0), (x1, y1) = to_km(x0, y0), to_km(x1, y1)
    return x0, y0, x1, y1


fig, ax = plt.subplots(figsize=(7.6, 5.2))
ck = grid.CHUNK_M / 1000.0
pad = 0.5
for cy in range(min(cys), max(cys) + 2):
    _, y = to_km(0, grid.Y_TOP - cy * grid.CHUNK_M)
    ax.axhline(y, color=GRID, lw=0.8, zorder=0)
for cx in range(min(cxs), max(cxs) + 2):
    x, _ = to_km(grid.X0 + cx * grid.CHUNK_M, 0)
    ax.axvline(x, color=GRID, lw=0.8, zorder=0)

for (cy, cx) in allc:
    x, y = to_km(grid.X0 + cx * grid.CHUNK_M, grid.Y_TOP - (cy + 1) * grid.CHUNK_M)
    if (cy, cx) in shared:
        fc = '#bcd9f5'
    elif (cy, cx) in ca:
        fc = '#e3eefb'
    else:
        fc = '#fbe6dd'
    ax.add_patch(Rectangle((x, y), ck, ck, facecolor=fc, edgecolor='none', zorder=1))

for b, col in ((bbox_a, BLUE), (bbox_b, ORANGE)):
    x0, y0, x1, y1 = bbox_km(b)
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, facecolor='none',
                           edgecolor=col, lw=2, zorder=3))
xa = bbox_km(bbox_a)
xb = bbox_km(bbox_b)
ax.text(xa[0], xa[3] + 0.09, 'troi A', color=BLUE, fontsize=9.5, weight='bold')
ax.text(xb[2], xb[1] - 0.24, 'troi B  (≈2 km E, ≈1.3 km S)', color=ORANGE,
        fontsize=9.5, weight='bold', ha='right')

handles = [
    Patch(facecolor='#e3eefb', label=f'chunks only A ({len(ca - cb)})'),
    Patch(facecolor='#fbe6dd', label=f'chunks only B ({len(cb - ca)})'),
    Patch(facecolor='#bcd9f5', label=f'shared — downloaded once ({len(shared)})'),
    Line2D([], [], color=BLUE, lw=2, label='troi A bbox'),
    Line2D([], [], color=ORANGE, lw=2, label='troi B bbox'),
]
ax.legend(handles=handles, frameon=False, fontsize=8.4, loc='center left',
          bbox_to_anchor=(1.01, 0.5))
ax.set_xlim(-pad, (max(cxs) - min(cxs) + 1) * ck + pad)
ax.set_ylim(-pad, (max(cys) - min(cys) + 1) * ck + pad)
ax.set_aspect('equal')
ax.set_xlabel('easting within the window (km) — gridlines are 2.56 km chunk edges')
ax.set_ylabel('northing within the window (km)')
ax.tick_params(labelsize=7.5)
ax.set_title('Two overlapping queries resolve to overlapping chunk sets',
             loc='left', fontsize=11, color=INK, pad=10)
for sp in ('top', 'right'):
    ax.spines[sp].set_visible(False)
fig.tight_layout()
fig.savefig(f'{OUT}/grid_dedup.png', bbox_inches='tight')
plt.close(fig)
print('fig 6 done', flush=True)
print('ALL DONE')
