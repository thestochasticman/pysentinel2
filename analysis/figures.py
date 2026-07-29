"""Figures for the parameter-selection analysis (reads analysis/data/*.csv).

Run after harness.py:  python analysis/figures.py
"""
import os
import csv
import numpy as np
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA, FIGS = os.path.join(HERE, 'data'), os.path.join(HERE, 'figures')
os.makedirs(FIGS, exist_ok=True)

SURFACE, INK, INK2, MUTED = '#fcfcfb', '#0b0b0b', '#52514e', '#898781'
GRID, BASE = '#e1e0d9', '#c3c2b7'
BLUE, ORANGE, AQUA = '#2a78d6', '#eb6834', '#1baf7a'
BLUE_LIGHT = '#cde2fb'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 9, 'text.color': INK,
    'axes.edgecolor': BASE, 'axes.labelcolor': INK2,
    'xtick.color': MUTED, 'ytick.color': MUTED,
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
    'savefig.facecolor': SURFACE, 'savefig.dpi': 150,
})

REGION_LABELS = {
    'alpine_snow': 'Alpine (Perisher, winter)',
    'tas_west_cloud': 'West Tasmania (winter)',
    'wet_tropics': 'Wet tropics (Tully, monsoon)',
    'arid_control': 'Arid (Alice Springs)',
    'coastal_mixed': 'Coastal (Jervis Bay)',
    'cropping_control': 'Cropping (Grenfell)',
}
ORDER = list(REGION_LABELS)

with open(f'{DATA}/frames.csv') as f:
    frames = [
        {k: (v if k in ('region', 'solar_day') else (None if v == '' else float(v)))
         for k, v in row.items()}
        for row in csv.DictReader(f)
    ]
with open(f'{DATA}/rings.csv') as f:
    rings = [
        {k: (v if k in ('region', 'solar_day') else float(v)) for k, v in row.items()}
        for row in csv.DictReader(f)
    ]


def style(ax):
    ax.grid(color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)


# ---- 1. scene-level vs window-level cloud ---------------------------------
fig, axes = plt.subplots(2, 3, figsize=(11.4, 7.2), sharex=True, sharey=True)
for ax, region in zip(axes.flat, ORDER):
    rows = [r for r in frames if r['region'] == region and r['scene_cloud'] is not None]
    sc = np.array([r['scene_cloud'] for r in rows])
    wc = np.array([100 * r['cloud_frac'] for r in rows])
    vf = np.array([r['valid_frac'] for r in rows])
    usable = vf >= 0.2
    ax.axvspan(30, 102, color='#f3f2ef', zorder=0)
    ax.axvline(30, color=BASE, lw=1, ls=(0, (4, 3)))
    ax.axhline(50, color=BASE, lw=1, ls=(0, (4, 3)))
    ax.scatter(sc[usable], wc[usable], s=26, color=BLUE, zorder=3)
    ax.scatter(sc[~usable], wc[~usable], s=26, color=MUTED, marker='X', zorder=3)
    lost = int(((sc > 30) & (wc <= 50) & usable).sum())
    total_usable = int((usable & (wc <= 50)).sum())
    ax.set_title(REGION_LABELS[region], loc='left', fontsize=9.5, color=INK)
    ax.text(0.97, 0.03, f'usable frames excluded\nby scene filter: {lost}/{total_usable}',
            transform=ax.transAxes, fontsize=7.6, color=INK2, ha='right', va='bottom')
    ax.set_xlim(-3, 102)
    ax.set_ylim(-3, 102)
    style(ax)
for ax in axes[1]:
    ax.set_xlabel('scene eo:cloud_cover (%)')
for ax in axes[:, 0]:
    ax.set_ylabel('window cloud fraction (%)')
fig.suptitle('Scene-level eo:cloud_cover versus window-level cloud fraction\n'
             'blue = frame with ≥20 % valid pixels; grey ✕ = off-swath frame; '
             'shaded = excluded by max_cloud_cover = 30; dashed = the two thresholds',
             fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(f'{FIGS}/scene_vs_window.png', bbox_inches='tight')
plt.close(fig)

# ---- 2. frame yield vs max_cloud_fraction ---------------------------------
sweep = np.linspace(0, 1, 51)
fig, axes = plt.subplots(2, 3, figsize=(11.4, 6.6), sharex=True, sharey=True)
for ax, region in zip(axes.flat, ORDER):
    rows = [r for r in frames if r['region'] == region]
    wc = np.array([r['cloud_frac'] for r in rows])
    vf = np.array([r['valid_frac'] for r in rows])
    n_usable = max(int((vf >= 0.2).sum()), 1)
    kept = [int(((vf >= 0.2) & (wc <= t)).sum()) for t in sweep]
    ax.plot(sweep, kept, color=BLUE, lw=2)
    ax.axvline(0.5, color=BASE, lw=1, ls=(0, (4, 3)))
    k05 = int(((vf >= 0.2) & (wc <= 0.5)).sum())
    ax.set_title(f'{REGION_LABELS[region]} — {k05}/{len(rows)} kept at 0.5',
                 loc='left', fontsize=9.3, color=INK)
    style(ax)
for ax in axes[1]:
    ax.set_xlabel('max_cloud_fraction threshold')
for ax in axes[:, 0]:
    ax.set_ylabel('frames kept')
fig.suptitle('Frame yield as a function of the cloud-fraction gate '
             '(min_valid_fraction = 0.2 applied throughout)',
             fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(f'{FIGS}/yield_curves.png', bbox_inches='tight')
plt.close(fig)

# ---- 3. near-cloud reflectance bias ---------------------------------------
SHOW = [('nbart_blue', BLUE, 'blue (B02)'),
        ('nbart_nir_1', AQUA, 'NIR (B08)'),
        ('nbart_swir_2', ORANGE, 'SWIR (B11)')]
fig, ax = plt.subplots(figsize=(7.6, 4.4))
dists = sorted({int(r['distance_px']) for r in rings})
for band, col, label in SHOW:
    med, lo, hi = [], [], []
    for d in dists:
        v = np.array([r[band] for r in rings if int(r['distance_px']) == d])
        med.append(np.median(v))
        lo.append(np.percentile(v, 25))
        hi.append(np.percentile(v, 75))
    ax.fill_between(dists, lo, hi, color=col, alpha=0.12, lw=0)
    ax.plot(dists, med, color=col, lw=2, marker='o', ms=4, label=label)
ax.axhline(1.0, color=BASE, lw=1)
ax.axvline(3, color=BASE, lw=1, ls=(0, (4, 3)))
ax.text(3.15, ax.get_ylim()[0] + 0.01, 'default buffer_px = 3', fontsize=8, color=INK2)
ax.set_xlabel('distance of clear pixel from nearest fmask cloud/shadow (px, 10 m each)')
ax.set_ylabel('median reflectance ÷ same-frame far-field (> 12 px)')
ax.set_title('Near-cloud reflectance bias in clear-classified pixels',
             loc='left', fontsize=11, color=INK, pad=10)
ax.legend(frameon=False, fontsize=8.5)
style(ax)
fig.tight_layout()
fig.savefig(f'{FIGS}/cloud_edge_bias.png', bbox_inches='tight')
plt.close(fig)

# ---- 4. per-region fmask composition --------------------------------------
comps = [('cloud_only_frac', '#eda100', 'cloud'), ('shadow_frac', '#4a3aa7', 'shadow'),
         ('snow_frac', '#e87ba4', 'snow'), ('water_frac', BLUE, 'water')]
fig, ax = plt.subplots(figsize=(8.6, 4.2))
x = np.arange(len(ORDER))
bottom = np.zeros(len(ORDER))
for key, col, label in comps:
    vals = []
    for region in ORDER:
        rows = [r for r in frames if r['region'] == region and r['valid_frac'] >= 0.2]
        vals.append(100 * float(np.mean([r[key] for r in rows])) if rows else 0.0)
    ax.bar(x, vals, 0.62, bottom=bottom, color=col, label=label,
           edgecolor=SURFACE, linewidth=2)
    bottom += np.array(vals)
ax.set_xticks(x)
ax.set_xticklabels([REGION_LABELS[r].split(' (')[0] for r in ORDER], fontsize=8.2)
ax.set_ylabel('mean share of valid pixels (%)')
ax.set_title('Mean fmask composition of valid pixels by region (frames with ≥20 % valid)',
             loc='left', fontsize=11, color=INK, pad=10)
ax.legend(frameon=False, fontsize=8.5, ncol=4, loc='upper right')
style(ax)
fig.tight_layout()
fig.savefig(f'{FIGS}/region_composition.png', bbox_inches='tight')
plt.close(fig)

print('figures written to', FIGS)
