"""The fixed global grid every Sentinel-2 pixel in the cube lives on.

One grid, defined once: EPSG:6933 (equal-area), 10 m pixels, 256x256-pixel
chunks (~2.56 km squares). Any EPSG:4326 bbox maps deterministically to a
pixel window and a set of chunk ids on this grid, which is what makes the
cube dedup-able: two overlapping queries resolve to overlapping chunk sets,
and a chunk is only ever downloaded once.

All functions here are pure — no I/O, no store access — so the grid math
is unit-testable offline.
"""
from pyproj import Transformer

CRS = 'EPSG:6933'
RES = 10                      # metres per pixel
CHUNK = 256                   # pixels per chunk edge
CHUNK_M = CHUNK * RES         # metres per chunk edge (2560)

# Grid origin: top-left corner, aligned to whole chunks and covering the
# full EPSG:6933 valid area (x within ±17,367,530 m, y within ±7,314,541 m).
X0 = -6785 * CHUNK_M          # -17,369,600
Y_TOP = 2862 * CHUNK_M        # 7,326,720
WIDTH_PX = 2 * 6785 * CHUNK   # 3,473,920
HEIGHT_PX = 2 * 2862 * CHUNK  # 1,465,344

_to_6933 = Transformer.from_crs('EPSG:4326', CRS, always_xy=True)
_to_4326 = Transformer.from_crs(CRS, 'EPSG:4326', always_xy=True)


def bbox_to_6933(bbox: list[float]) -> tuple[float, float, float, float]:
    """Project an EPSG:4326 ``[west, south, east, north]`` bbox to EPSG:6933.

    EPSG:6933 is cylindrical (x depends only on lon, y only on lat, both
    monotonic), so transforming the two corners is exact.
    """
    west, south, east, north = bbox
    x0, y0 = _to_6933.transform(west, south)
    x1, y1 = _to_6933.transform(east, north)
    return (x0, y0, x1, y1)


def window_for_bbox(bbox: list[float]) -> tuple[int, int, int, int]:
    """Pixel window ``(row0, row1, col0, col1)`` covering ``bbox``, snapped
    outward to whole chunks.

    Rows increase downward (row 0 at ``Y_TOP``); the half-open window
    ``[row0, row1) x [col0, col1)`` is always a whole number of chunks so
    every write touches complete chunks only — the unit of dedup.
    """
    x0, y0, x1, y1 = bbox_to_6933(bbox)
    col0 = int((x0 - X0) // CHUNK_M) * CHUNK
    col1 = -int(-(x1 - X0) // CHUNK_M) * CHUNK   # ceil-div
    row0 = int((Y_TOP - y1) // CHUNK_M) * CHUNK
    row1 = -int(-(Y_TOP - y0) // CHUNK_M) * CHUNK
    return (row0, row1, col0, col1)


def tight_window_for_bbox(bbox: list[float]) -> tuple[int, int, int, int]:
    """Pixel window ``(row0, row1, col0, col1)`` covering ``bbox`` exactly —
    snapped outward to whole pixels only, not to chunks.

    :func:`window_for_bbox` (chunk-snapped) is the unit of storage and
    dedup accounting; this is the unit of *delivery*. A chunk-snapped
    window hands consumers up to (2*CHUNK-1)^2 extra pixels of padding —
    3.6x the requested area for a ~2.5 km bbox — and every downstream
    compute (cleaning, unmixing, segmentation) pays for it.
    """
    x0, y0, x1, y1 = bbox_to_6933(bbox)
    col0 = int((x0 - X0) // RES)
    col1 = -int(-(x1 - X0) // RES)
    row0 = int((Y_TOP - y1) // RES)
    row1 = -int(-(Y_TOP - y0) // RES)
    return (row0, row1, col0, col1)


def rect_subtract(rect: tuple[int, int, int, int],
                  covers: list[tuple[int, int, int, int]],
                  ) -> list[tuple[int, int, int, int]]:
    """Disjoint rectangles of ``rect`` not covered by any rect in ``covers``.

    All rects are half-open pixel windows ``(row0, row1, col0, col1)``.
    Each subtraction splits a remainder into at most four pieces (above,
    below, left, right of the intersection). The basis of pixel-exact
    coverage accounting: a fill's missing work for a day is
    ``rect_subtract(tight_window, covered_rects(day))``.
    """
    remaining = [rect] if rect[0] < rect[1] and rect[2] < rect[3] else []
    for cr0, cr1, cc0, cc1 in covers:
        next_remaining = []
        for r0, r1, c0, c1 in remaining:
            ir0, ir1 = max(r0, cr0), min(r1, cr1)
            ic0, ic1 = max(c0, cc0), min(c1, cc1)
            if ir0 >= ir1 or ic0 >= ic1:          # no overlap
                next_remaining.append((r0, r1, c0, c1))
                continue
            if r0 < ir0:                           # strip above
                next_remaining.append((r0, ir0, c0, c1))
            if ir1 < r1:                           # strip below
                next_remaining.append((ir1, r1, c0, c1))
            if c0 < ic0:                           # strip left of intersection
                next_remaining.append((ir0, ir1, c0, ic0))
            if ic1 < c1:                           # strip right
                next_remaining.append((ir0, ir1, ic1, c1))
        remaining = next_remaining
    return remaining


def rects_bbox(rects: list[tuple[int, int, int, int]]) -> tuple[int, int, int, int]:
    """Smallest pixel window containing every rect in ``rects``."""
    return (
        min(r[0] for r in rects), max(r[1] for r in rects),
        min(r[2] for r in rects), max(r[3] for r in rects),
    )


def chunks_in_window(window: tuple[int, int, int, int]) -> list[tuple[int, int]]:
    """All chunk ids ``(cy, cx)`` inside a chunk-aligned pixel window."""
    row0, row1, col0, col1 = window
    return [
        (cy, cx)
        for cy in range(row0 // CHUNK, row1 // CHUNK)
        for cx in range(col0 // CHUNK, col1 // CHUNK)
    ]


def window_of_chunks(chunks: list[tuple[int, int]]) -> tuple[int, int, int, int]:
    """Smallest chunk-aligned pixel window containing every chunk in ``chunks``."""
    cys = [cy for cy, _ in chunks]
    cxs = [cx for _, cx in chunks]
    return (
        min(cys) * CHUNK, (max(cys) + 1) * CHUNK,
        min(cxs) * CHUNK, (max(cxs) + 1) * CHUNK,
    )


def window_to_4326(window: tuple[int, int, int, int]) -> list[float]:
    """EPSG:4326 ``[west, south, east, north]`` bbox of a pixel window."""
    row0, row1, col0, col1 = window
    x0, x1 = X0 + col0 * RES, X0 + col1 * RES
    y1, y0 = Y_TOP - row0 * RES, Y_TOP - row1 * RES
    west, south = _to_4326.transform(x0, y0)
    east, north = _to_4326.transform(x1, y1)
    return [west, south, east, north]


def coords_for_window(window: tuple[int, int, int, int]):
    """Pixel-centre coordinate arrays ``(y, x)`` for a window (y descending)."""
    import numpy as np
    row0, row1, col0, col1 = window
    x = X0 + (np.arange(col0, col1) + 0.5) * RES
    y = Y_TOP - (np.arange(row0, row1) + 0.5) * RES
    return y, x


def geobox_for_window(window: tuple[int, int, int, int]):
    """An ``odc.geo.geobox.GeoBox`` for a window — pinning odc.stac.load to
    exactly this grid so downloaded pixels land chunk-aligned."""
    from affine import Affine
    from odc.geo.geobox import GeoBox
    row0, row1, col0, col1 = window
    transform = Affine(RES, 0, X0 + col0 * RES, 0, -RES, Y_TOP - row0 * RES)
    return GeoBox((row1 - row0, col1 - col0), transform, CRS)


def test_window_is_chunk_aligned():
    bbox = [148.36265, -33.52606, 148.38265, -33.50606]
    row0, row1, col0, col1 = window_for_bbox(bbox)
    return all(v % CHUNK == 0 for v in (row0, row1, col0, col1)) and row1 > row0 and col1 > col0


def test_window_contains_bbox():
    bbox = [148.36265, -33.52606, 148.38265, -33.50606]
    x0, y0, x1, y1 = bbox_to_6933(bbox)
    row0, row1, col0, col1 = window_for_bbox(bbox)
    return (
        X0 + col0 * RES <= x0 and X0 + col1 * RES >= x1
        and Y_TOP - row0 * RES >= y1 and Y_TOP - row1 * RES <= y0
    )


def test_overlapping_bboxes_share_chunks():
    a = window_for_bbox([148.36265, -33.52606, 148.38265, -33.50606])
    b = window_for_bbox([148.37265, -33.52606, 148.39265, -33.50606])  # shifted ~1km east
    shared = set(chunks_in_window(a)) & set(chunks_in_window(b))
    return len(shared) > 0


def test_roundtrip_4326():
    bbox = [148.36265, -33.52606, 148.38265, -33.50606]
    w = window_for_bbox(bbox)
    back = window_to_4326(w)
    # The window bbox must contain the original bbox
    return back[0] <= bbox[0] and back[1] <= bbox[1] and back[2] >= bbox[2] and back[3] >= bbox[3]


def test_coords_match_window():
    import numpy as np
    w = window_for_bbox([148.36265, -33.52606, 148.38265, -33.50606])
    y, x = coords_for_window(w)
    return len(y) == w[1] - w[0] and len(x) == w[3] - w[2] and y[0] > y[-1] and x[0] < x[-1]


def test_rect_subtract_against_mask():
    """rect_subtract must agree exactly with brute-force boolean masking
    across randomised cover sets: disjoint output, full missing coverage."""
    import numpy as np
    rng = np.random.default_rng(7)
    for _ in range(50):
        H = W = 40
        rect = (0, H, 0, W)
        covers = []
        for _ in range(rng.integers(0, 6)):
            r0, c0 = rng.integers(0, H - 1), rng.integers(0, W - 1)
            covers.append((int(r0), int(rng.integers(r0 + 1, H + 5)),
                           int(c0), int(rng.integers(c0 + 1, W + 5))))
        missing = rect_subtract(rect, covers)
        # reference mask
        mask = np.ones((H, W), dtype=bool)
        for r0, r1, c0, c1 in covers:
            mask[max(r0, 0):r1, max(c0, 0):c1] = False
        out = np.zeros((H, W), dtype=int)
        for r0, r1, c0, c1 in missing:
            out[r0:r1, c0:c1] += 1
        if not ((out <= 1).all()                      # disjoint
                and ((out == 1) == mask).all()):      # exactly the missing set
            return False
    return True


def test_rects_bbox():
    return rects_bbox([(2, 5, 1, 4), (0, 3, 3, 9)]) == (0, 5, 1, 9)


def test():
    return all([
        test_window_is_chunk_aligned(),
        test_window_contains_bbox(),
        test_overlapping_bboxes_share_chunks(),
        test_roundtrip_4326(),
        test_coords_match_window(),
        test_rect_subtract_against_mask(),
        test_rects_bbox(),
    ])


if __name__ == '__main__':
    print(test())
