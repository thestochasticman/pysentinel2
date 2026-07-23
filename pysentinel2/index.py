"""SQLite index of what the cube holds and what it has already looked for.

Three tables, no pixels:

- ``scenes`` — every STAC item ever seen (solar day, cloud cover, full item
  JSON), so day selection and re-fills work offline without re-hitting STAC.
- ``chunks`` — which ``(solar_day, cy, cx)`` cells of the grid are populated.
  A cell present here is on disk and complete; absence means never fetched.
  This is the dedup ledger: nothing in it is ever downloaded again.
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

    # -- chunks -----------------------------------------------------------

    def chunks_done(self, solar_day: str, chunks: list[tuple[int, int]]) -> set[tuple[int, int]]:
        """Subset of ``chunks`` already populated for ``solar_day``."""
        rows = self.db.execute(
            'SELECT cy, cx FROM chunks WHERE solar_day = ?', (solar_day,)
        ).fetchall()
        return set(rows) & set(chunks)

    def mark_chunks(self, solar_day: str, chunks: list[tuple[int, int]]) -> None:
        now = _now()
        with self.db:
            self.db.executemany(
                'INSERT OR REPLACE INTO chunks (solar_day, cy, cx, written_at) '
                'VALUES (?, ?, ?, ?)',
                [(solar_day, cy, cx, now) for cy, cx in chunks],
            )

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


def test_chunks_ledger():
    ix = _tmp_index()
    ix.mark_chunks('2024-01-03', [(1, 1), (1, 2)])
    done = ix.chunks_done('2024-01-03', [(1, 1), (1, 2), (1, 3)])
    other_day = ix.chunks_done('2024-01-08', [(1, 1)])
    return done == {(1, 1), (1, 2)} and other_day == set()


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
        test_chunks_ledger(),
        test_scenes_cloud_filter(),
    ])


if __name__ == '__main__':
    print(test())
