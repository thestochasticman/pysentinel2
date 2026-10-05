"""File index of what the cube holds and what it has already looked for.

Three marker trees under ``{root}/index/``, no pixels, no database:

- ``scenes/<YYYY>/<YYYY-MM-DD>/<item_id>.json`` -- every STAC item ever
  seen (solar day, cloud cover, full item JSON), so day selection and
  re-fills work offline without re-hitting STAC.
- ``coverage/<YYYY-MM-DD>/<uuid>.json`` -- which pixel rectangles of the
  grid are populated per solar day, one rect per file. A region covered
  here is on disk and complete; anything outside is never-fetched. This
  is the dedup ledger at pixel exactness: fills only ever write
  axis-aligned windows, so a day's coverage is a short list of rects,
  and missing work is ``grid.rect_subtract(window, covered)``.
- ``searches/<uuid>.json`` -- which (bbox, date-range) regions STAC has
  been queried for, so "no scenes" is distinguishable from "never
  asked".

Every marker is committed by write-to-temp + atomic rename
(:class:`troi.ledger.Markers`), so a reader sees it whole or not at
all, and a crash mid-fill leaves cells unmarked -- never a half-written
cell recorded as done. Files rather than SQLite because the cube is
filled by many PBS jobs on many Gadi nodes against one store on
Lustre, where file locks are node-local (see ``troi/docs/ledger.md``).
"""
import json
import os
import uuid
from datetime import date, datetime, timezone

from troi.ledger import Markers


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class Index:
    """The cube's ledger of scenes, coverage rects and searches, as files."""

    def __init__(self, path: str):
        self.path = path
        self.m = Markers(path)

    # -- searches ---------------------------------------------------------

    def search_covered(self, bbox6933: tuple, start: date, end: date) -> bool:
        """True iff a previous STAC search fully contains this bbox + range."""
        x0, y0, x1, y1 = bbox6933
        for leaf in self.m.list(('searches',)):
            s = self.m.read(('searches', leaf))
            if (s and s['x0'] <= x0 and s['y0'] <= y0 and s['x1'] >= x1 and s['y1'] >= y1
                    and s['start'] <= str(start) and s['end'] >= str(end)):
                return True
        return False

    def record_search(self, bbox6933: tuple, start: date, end: date) -> None:
        x0, y0, x1, y1 = bbox6933
        self.m.write(('searches', uuid.uuid4().hex),
                     {'x0': x0, 'y0': y0, 'x1': x1, 'y1': y1,
                      'start': str(start), 'end': str(end), 'searched_at': _now()})

    # -- scenes -----------------------------------------------------------

    def upsert_scenes(self, scenes: list[tuple[str, str, float, dict]]) -> None:
        """Insert/refresh ``(item_id, solar_day, cloud_cover, item_dict)`` markers."""
        for item_id, day, cc, item in scenes:
            self.m.write(('scenes', day[:4], day, item_id.replace('/', '_')),
                         {'item_id': item_id, 'solar_day': day, 'cloud_cover': cc, 'item': item})

    def scenes_for_range(self, start: date, end: date, max_cloud_cover: float) -> dict[str, list[dict]]:
        """Item dicts per solar day within ``[start, end]``, cloud-filtered."""
        by_day: dict[str, list[dict]] = {}
        for year in range(start.year, end.year + 1):
            ydir = os.path.join(self.path, 'scenes', str(year))
            try:
                days = sorted(os.listdir(ydir))
            except FileNotFoundError:
                continue
            for day in days:
                if not (str(start) <= day <= str(end)):
                    continue
                for leaf in self.m.list(('scenes', str(year), day)):
                    s = self.m.read(('scenes', str(year), day, leaf))
                    if s and (s['cloud_cover'] is None or s['cloud_cover'] <= max_cloud_cover):
                        by_day.setdefault(day, []).append(s['item'])
        return by_day

    # -- coverage ---------------------------------------------------------

    def covered_rects(self, solar_day: str) -> list[tuple[int, int, int, int]]:
        """Pixel rects ``(row0, row1, col0, col1)`` populated for the day,
        in the order they were written."""
        rects = []
        for leaf in self.m.list(('coverage', solar_day)):
            r = self.m.read(('coverage', solar_day, leaf))
            if r:
                rects.append((r['written_at'], (r['row0'], r['row1'], r['col0'], r['col1'])))
        return [rect for _, rect in sorted(rects)]

    def mark_rect(self, solar_day: str, rect: tuple[int, int, int, int]) -> None:
        """Record a written pixel window. Call strictly after the write."""
        row0, row1, col0, col1 = rect
        self.m.write(('coverage', solar_day, uuid.uuid4().hex),
                     {'row0': int(row0), 'row1': int(row1), 'col0': int(col0), 'col1': int(col1)})

    def unmark_day(self, solar_day: str) -> None:
        """Drop a day's coverage so the next fill re-downloads it.

        Used by :meth:`pysentinel2.cube.Cube.repair` when stored pixels
        turn out to be download failures rather than real data gaps.
        """
        for leaf in self.m.list(('coverage', solar_day)):
            self.m.remove(('coverage', solar_day, leaf))

    def close(self) -> None:
        pass


def _tmp_index():
    import tempfile
    return Index(f'{tempfile.mkdtemp(prefix="pysentinel2_index_test_")}/index')


def test_search_coverage():
    ix = _tmp_index()
    ix.record_search((0.0, 0.0, 100.0, 100.0), date(2024, 1, 1), date(2024, 6, 30))
    inside = ix.search_covered((10.0, 10.0, 90.0, 90.0), date(2024, 2, 1), date(2024, 3, 1))
    outside_space = ix.search_covered((10.0, 10.0, 200.0, 90.0), date(2024, 2, 1), date(2024, 3, 1))
    outside_time = ix.search_covered((10.0, 10.0, 90.0, 90.0), date(2024, 2, 1), date(2024, 7, 1))
    return inside and not outside_space and not outside_time


def test_coverage_ledger():
    ix = _tmp_index()
    ix.mark_rect('2024-01-03', (100, 200, 50, 300))
    ix.mark_rect('2024-01-03', (200, 250, 50, 300))
    same_day = ix.covered_rects('2024-01-03')
    other_day = ix.covered_rects('2024-01-08')
    ix.unmark_day('2024-01-03')
    return (sorted(same_day) == [(100, 200, 50, 300), (200, 250, 50, 300)]
            and other_day == [] and ix.covered_rects('2024-01-03') == [])


def test_markers_are_whole_files():
    """No temp files survive a write; a reopened Index sees everything."""
    ix = _tmp_index()
    ix.mark_rect('2024-01-03', (0, 10, 0, 10))
    ix.record_search((0, 0, 1, 1), date(2024, 1, 1), date(2024, 1, 2))
    ix.upsert_scenes([('a/b', '2024-01-03', 1.0, {'id': 'a/b'})])
    names = [f for _, _, fs in os.walk(ix.path) for f in fs]
    again = Index(ix.path)
    return (all(n.endswith('.json') for n in names) and len(names) == 3
            and again.covered_rects('2024-01-03') == [(0, 10, 0, 10)]
            and again.scenes_for_range(date(2024, 1, 1), date(2024, 1, 31), 50.0)['2024-01-03'][0]['id'] == 'a/b')


def test_scenes_cloud_filter():
    ix = _tmp_index()
    ix.upsert_scenes([
        ('item_a', '2024-01-03', 5.0, {'id': 'item_a'}),
        ('item_b', '2024-01-03', 80.0, {'id': 'item_b'}),
        ('item_c', '2024-01-08', 10.0, {'id': 'item_c'}),
        ('item_d', '2023-12-30', 1.0, {'id': 'item_d'}),
    ])
    by_day = ix.scenes_for_range(date(2024, 1, 1), date(2024, 1, 31), max_cloud_cover=30.0)
    across = ix.scenes_for_range(date(2023, 12, 1), date(2024, 1, 5), max_cloud_cover=30.0)
    return (
        [i['id'] for i in by_day.get('2024-01-03', [])] == ['item_a']
        and [i['id'] for i in by_day.get('2024-01-08', [])] == ['item_c']
        and sorted(across) == ['2023-12-30', '2024-01-03']
    )


def test():
    return all([
        test_search_coverage(),
        test_coverage_ledger(),
        test_markers_are_whole_files(),
        test_scenes_cloud_filter(),
    ])


if __name__ == '__main__':
    print(test())
