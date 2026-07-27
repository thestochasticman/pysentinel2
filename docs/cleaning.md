# Cleaning & masking

`get_ds(..., clean=True)` — and any `indices=` request, which implies
it — runs the window through `pysentinel2.cube.clean_dataset`. Nothing
here is persisted: the clean cube is a view of the raw store, so
different thresholds on the same window are simply different reads, and
storage cost is roughly half that of a raw-plus-clean layout.

The masking source is DEA's `oa_fmask` band — a per-pixel classification
produced by the Fmask algorithm (Zhu & Woodcock, 2012; Frantz et al.,
2018, as implemented in DEA's ARD pipeline).

## Invalid versus contaminated pixels

A pixel can be unusable for two distinct reasons, and conflating them
biases frame selection:

| State | fmask | Physical meaning | Treatment |
|---|---|---|---|
| Invalid | 0 (nodata) | Outside the swath / never sensed | → NaN; counts against **coverage**, not cloudiness |
| Clear | 1 | Usable land observation | kept |
| Cloud | 2 | Contaminated radiometry | → NaN (dilated) |
| Cloud shadow | 3 | Contaminated radiometry | → NaN (dilated) |
| Snow | 4 | Corrupts reflectance statistics like cloud | → NaN by default (`mask_snow=False` to keep) |
| Water | 5 | Legitimate surface signal (dams, rivers, NDWI) | kept by default (`mask_water=True` to drop) |

A frame whose window is 60 % off-swath but cloud-free over the
remaining 40 % is a valid observation of 40 % of the window. A rule
that scores frames by total NaN fraction — as an earlier implementation
of this package did — discards such frames and systematically biases
the record against swath-edge regions.

## Pipeline

```mermaid
flowchart TD
    A["raw window (time, y, x)<br/>11 bands incl. oa_fmask"] --> B["1 · classify per pixel:<br/>valid = fmask ≠ 0<br/>contaminated = cloud ∪ shadow (∪ snow ∪ water)"]
    B --> C["2 · dilate the contaminated mask<br/>by buffer_px (circular element, default 3 px)"]
    C --> D["3 · per-frame gates:<br/>cloud_fraction = contaminated ÷ valid ≤ 0.5<br/>valid_fraction = valid ÷ window ≥ 0.2"]
    D -- "frame fails either gate" --> X["frame dropped"]
    D -- "frame passes" --> E["mask contaminated + invalid → NaN;<br/>drop the fmask band; mask band nodata (−999)"]
    E --> F["4 · annotate: cloud_fraction & valid_fraction<br/>as time coordinates; thresholds as dataset attrs"]
```

**Step 2 — dilation.** Fmask draws tight cloud boundaries, and the
bright halo and penumbra immediately outside a detected cloud are a
well-documented source of residual contamination in nominally clear
pixels (Zhu & Woodcock, 2012, recommend buffering for this reason). The
contaminated mask is therefore dilated with a circular structuring
element of radius `buffer_px` (default 3 px ≈ 30 m) before masking.
Dilation does not convert invalid pixels to contaminated: the mask is
intersected with the valid mask afterwards, so swath margins remain
classified as invalid rather than cloudy.

**Step 3 — two independent gates.** For each frame:

$$
\mathrm{cloud\_fraction} = \frac{\#\,\text{contaminated}}{\#\,\text{valid}}
\qquad
\mathrm{valid\_fraction} = \frac{\#\,\text{valid}}{\#\,\text{window pixels}}
$$

A frame is retained when both
$\mathrm{cloud\_fraction} \le$ `max_cloud_fraction` (default 0.5) and
$\mathrm{valid\_fraction} \ge$ `min_valid_fraction` (default 0.2) hold.
Because the cloudiness denominator is valid pixels, partial-swath
frames are not penalised for their margin; because coverage is gated
separately, a clear but near-empty sliver of swath is still rejected as
an unusable observation.

## The pipeline on a real frame

The 2023-12-23 overpass of the example window is 72 % cloud-contaminated
over its valid pixels:

![On-read cleaning of a cloudy frame](images/cleaning_pipeline.png)

Panel b is the raw fmask classification; panel c overlays the
contaminated mask after dilation, with the thin ring marking the
buffer; panel d shows the masked observation as `clean=True` would
return it, had the frame passed the gates — at the default
`max_cloud_fraction = 0.5` this frame is dropped.

## Frame gating on the example window

Across the 12 stored solar days, the two statistics separate the frames
into three groups: fully valid clear or partly cloudy frames (kept),
one heavily clouded frame (dropped by the cloud gate), and six days on
which the swath missed the window (dropped by the valid gate, rather
than being misrecorded as fully cloudy days):

![Frame gating scatter](images/frame_gates.png)

## Auditability

Every survival decision is reconstructible from the returned dataset
alone:

```python
ds = cube.get_ds(bbox, start, end, clean=True)

ds.cloud_fraction          # (time,) — contamination of each surviving frame
ds.valid_fraction          # (time,) — coverage of each surviving frame
ds.attrs['max_cloud_fraction']    # the gates this read was made with
ds.attrs['min_valid_fraction']
ds.attrs['cloud_buffer_px']
ds.attrs['masked_fmask_classes']  # e.g. [2, 3, 4]
```

All parameters are per-read:

```python
ds = cube.get_ds(bbox, start, end, clean=True,
                 max_cloud_fraction=0.3,  # stricter contamination gate
                 min_valid_fraction=0.5,  # require ≥ half the window sensed
                 mask_water=True,         # pure-vegetation statistics
                 buffer_px=5)             # wider halo exclusion
```

## References

- Zhu, Z. & Woodcock, C. E. (2012). Object-based cloud and cloud shadow
  detection in Landsat imagery. *Remote Sensing of Environment*, 118,
  83–94.
- Frantz, D. et al. (2018). Improvement of the Fmask algorithm for
  Sentinel-2 images. *Remote Sensing of Environment*, 215, 471–481.
- Digital Earth Australia, *Sentinel-2 ARD* product documentation
  (`oa_fmask` observation attribute).
