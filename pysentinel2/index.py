"""SQLite index of what the cube holds and what it has already looked for.

Three tables, no pixels:

- ``scenes`` — every STAC item ever seen (solar day, cloud cover, full item
  JSON), so day selection and re-fills work offline without re-hitting STAC.
- ``coverage`` — which pixel rectangles of the grid are populated per solar
  day. A region covered here is on disk and complete; anything outside is
  never-fetched. This is the dedup ledger at pixel exactness: fills only
  ever write axis-aligned windows, so a day's coverage is a short list of
  rects, and missing work is ``grid.rect_subtract(window, covered)``.
  (Earlier versions accounted in fixed 256 px chunks — the ``chunks``
  table — which forced downloads to pad small farms by up to ~3.6x; old
  ledgers migrate automatically on open.)
- ``searches`` — which (bbox, date-range) regions STAC has been queried
  for, so "no scenes" is distinguishable from "never asked".

SQLite (WAL mode) makes completeness transactional: a crash mid-fill
leaves cells unmarked — never a half-written cell recorded as done.
"""
import json
import sqlite3
from datetime import date, datetime, timezone

_SCHEMA = """
CREATE TABLE IF NOT EXISTS scenes (
    item_id     TEXT PRIMARY KEY,
    solar_day   TEXT NOT NULL,
    cloud_cover REAL,
    item_json   TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS scenes_by_day ON scenes(solar_day);

CREATE TABLE IF NOT EXISTS chunks (
    solar_day  TEXT NOT NULL,
    cy         INTEGER NOT NULL,
    cx         INTEGER NOT NULL,
    written_at TEXT NOT NULL,
    PRIMARY KEY (solar_day, cy, cx)
) WITHOUT ROWID;

CREATE TABLE IF NOT EXISTS coverage (
    solar_day  TEXT NOT NULL,
    row0       INTEGER NOT NULL,
    row1       INTEGER NOT NULL,
    col0       INTEGER NOT NULL,
    col1       INTEGER NOT NULL,
    written_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS coverage_by_day ON coverage(solar_day);

CREATE TABLE IF NOT EXISTS searches (
    x0 REAL NOT NULL, y0 REAL NOT NULL,
    x1 REAL NOT NULL, y1 REAL NOT NULL,
    start TEXT NOT NULL, end TEXT NOT NULL,
    searched_at TEXT NOT NULL
);
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


class Index:
    """Thin, transaction-safe wrapper around the cube's SQLite database."""

    def __init__(self, path: str):
        self.path = path
        self.db = sqlite3.connect(path)
        self.db.execute('PRAGMA journal_mode=WAL')
        self.db.executescript(_SCHEMA)
        self.db.commit()
        self._migrate_chunks_to_coverage()

    def _migrate_chunks_to_coverage(self) -> None:
        """One-time, idempotent: convert legacy 256 px chunk marks to rects.

        Each ``(day, cy, cx)`` row is the rect ``(cy*256, cy*256+256,
        cx*256, cx*256+256)``; runs of chunks adjacent along x coalesce
        into one rect per (day, cy) run — exactness matters, minimality
        does not. Runs only when the legacy table has rows and coverage
        is empty, all inside one transaction; legacy rows are kept (the
        table is dead weight, dropped in a later release).
        """
        has_legacy = self.db.execute('SELECT 1 FROM chunks LIMIT 1').fetchone()
        has_coverage = self.db.execute('SELECT 1 FROM coverage LIMIT 1').fetchone()
        if not has_legacy or has_coverage:
            return
        CHUNK = 256
        now = _now()
        rows = self.db.execute(
            'SELECT solar_day, cy, cx FROM chunks ORDER BY solar_day, cy, cx'
        ).fetchall()
        rects = []
        run = None   # (day, cy, cx_start, cx_end_exclusive)
        for day, cy, cx in rows:
            if run is not None and run[0] == day and run[1] == cy and run[3] == cx:
                run = (day, cy, run[2], cx + 1)
                continue
            if run is not None:
                d, y, x0, x1 = run
                rects.append((d, y * CHUNK, (y + 1) * CHUNK, x0 * CHUNK, x1 * CHUNK, now))
            run = (day, cy, cx, cx + 1)
        if run is not None:
            d, y, x0, x1 = run
            rects.append((d, y * CHUNK, (y + 1) * CHUNK, x0 * CHUNK, x1 * CHUNK, now))
        with self.db:
            self.db.executemany(
                'INSERT INTO coverage (solar_day, row0, row1, col0, col1, written_at) '
                'VALUES (?, ?, ?, ?, ?, ?)', rects)
        print(f'pysentinel2: migrated {len(rows)} legacy chunk marks '
              f'to {len(rects)} coverage rects')

    # -- searches ---------------------------------------------------------

    def search_covered(self, bbox6933: tuple, start: date, end: date) -> bool:
        """True iff a previous STAC search fully contains this bbox + range."""
        x0, y0, x1, y1 = bbox6933
        row = self.db.execute(
            'SELECT 1 FROM searches WHERE x0<=? AND y0<=? AND x1>=? AND y1>=? '
            'AND start<=? AND end>=? LIMIT 1',
            (x0, y0, x1, y1, str(start), str(end)),
        ).fetchone()
        return row is not None

    def record_search(self, bbox6933: tuple, start: date, end: date) -> None:
        x0, y0, x1, y1 = bbox6933
        with self.db:
            self.db.execute(
                'INSERT INTO searches (x0, y0, x1, y1, start, end, searched_at) '
                'VALUES (?, ?, ?, ?, ?, ?, ?)',
                (x0, y0, x1, y1, str(start), str(end), _now()),
            )

    # -- scenes -----------------------------------------------------------

    def upsert_scenes(self, scenes: list[tuple[str, str, float, dict]]) -> None:
        """Insert/refresh ``(item_id, solar_day, cloud_cover, item_dict)`` rows."""
        with self.db:
            self.db.executemany(
                'INSERT OR REPLACE INTO scenes (item_id, solar_day, cloud_cover, item_json) '
                'VALUES (?, ?, ?, ?)',
                [(i, d, cc, json.dumps(j)) for i, d, cc, j in scenes],
            )

    def scenes_for_range(self, start: date, end: date, max_cloud_cover: float) -> dict[str, list[dict]]:
        """Item dicts per solar day within ``[start, end]``, cloud-filtered."""
        rows = self.db.execute(
            'SELECT solar_day, item_json FROM scenes '
            'WHERE solar_day >= ? AND solar_day <= ? AND (cloud_cover IS NULL OR cloud_cover <= ?) '
            'ORDER BY solar_day',
            (str(start), str(end), max_cloud_cover),
        ).fetchall()
        by_day: dict[str, list[dict]] = {}
        for day, item_json in rows:
            by_day.setdefault(day, []).append(json.loads(item_json))
        return by_day

    # -- coverage ---------------------------------------------------------

    def covered_rects(self, solar_day: str) -> list[tuple[int, int, int, int]]:
        """Pixel rects ``(row0, row1, col0, col1)`` populated for the day."""
        return self.db.execute(
            'SELECT row0, row1, col0, col1 FROM coverage WHERE solar_day = ?',
            (solar_day,),
        ).fetchall()

    def mark_rect(self, solar_day: str, rect: tuple[int, int, int, int]) -> None:
        """Record a written pixel window. Call strictly after the write."""
        row0, row1, col0, col1 = rect
        with self.db:
            self.db.execute(
                'INSERT INTO coverage (solar_day, row0, row1, col0, col1, written_at) '
                'VALUES (?, ?, ?, ?, ?, ?)',
                (solar_day, row0, row1, col0, col1, _now()),
            )

    def unmark_day(self, solar_day: str) -> None:
        """Drop a day's coverage so the next fill re-downloads it.

        Used by :meth:`pysentinel2.cube.Cube.repair` when stored pixels
        turn out to be download failures rather than real data gaps.
        """
        with self.db:
            self.db.execute('DELETE FROM coverage WHERE solar_day = ?', (solar_day,))

    def close(self) -> None:
        self.db.close()


def _tmp_index():
    import tempfile
    return Index(f'{tempfile.mkdtemp(prefix="pysentinel2_index_test_")}/index.db')


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
    return (same_day == [(100, 200, 50, 300), (200, 250, 50, 300)]
            and other_day == [] and ix.covered_rects('2024-01-03') == [])


def test_chunk_migration():
    """Legacy chunk rows must convert to rects covering identical pixels."""
    import numpy as np
    ix = _tmp_index()
    legacy = [('2024-01-03', 1, 1), ('2024-01-03', 1, 2), ('2024-01-03', 3, 7),
              ('2024-01-08', 2, 2)]
    with ix.db:
        ix.db.executemany(
            "INSERT INTO chunks (solar_day, cy, cx, written_at) VALUES (?, ?, ?, 't')",
            legacy)
    ix.close()
    ix = Index(ix.path)   # reopen triggers migration
    for day in ('2024-01-03', '2024-01-08'):
        want = np.zeros((1200, 2400), dtype=bool)
        for d, cy, cx in legacy:
            if d == day:
                want[cy * 256:(cy + 1) * 256, cx * 256:(cx + 1) * 256] = True
        got = np.zeros_like(want)
        for r0, r1, c0, c1 in ix.covered_rects(day):
            got[r0:r1, c0:c1] = True
        if not (got == want).all():
            return False
    # idempotent: reopening again must not duplicate rects
    n = len(ix.covered_rects('2024-01-03'))
    ix.close()
    ix = Index(ix.path)
    return len(ix.covered_rects('2024-01-03')) == n


def test_scenes_cloud_filter():
    ix = _tmp_index()
    ix.upsert_scenes([
        ('item_a', '2024-01-03', 5.0, {'id': 'item_a'}),
        ('item_b', '2024-01-03', 80.0, {'id': 'item_b'}),
        ('item_c', '2024-01-08', 10.0, {'id': 'item_c'}),
    ])
    by_day = ix.scenes_for_range(date(2024, 1, 1), date(2024, 1, 31), max_cloud_cover=30.0)
    return (
        [i['id'] for i in by_day.get('2024-01-03', [])] == ['item_a']
        and [i['id'] for i in by_day.get('2024-01-08', [])] == ['item_c']
    )


def test():
    return all([
        test_search_coverage(),
        test_coverage_ledger(),
        test_chunk_migration(),
        test_scenes_cloud_filter(),
    ])


if __name__ == '__main__':
    print(test())
