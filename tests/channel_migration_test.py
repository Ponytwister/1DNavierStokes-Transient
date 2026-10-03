import importlib.util
from contextlib import closing
from pathlib import Path
import sqlite3
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("migration", ROOT / "tools/migrate_channel_dimensions.py")
migration = importlib.util.module_from_spec(spec)
spec.loader.exec_module(migration)

class MigrationTest(unittest.TestCase):
    def test_defaults_preservation_backup_and_repeat(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "experiments.db"
            backup = Path(folder) / "backup.db"
            with closing(sqlite3.connect(str(path))) as db, db:
                db.executescript((ROOT / "tests/fixtures/workflow.sql").read_text())
                db.execute("INSERT INTO experiments(NAME) VALUES('second')")
            with closing(sqlite3.connect(str(path))) as db, db:
                tables = [r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")]
                before = {t: db.execute(f'SELECT * FROM "{t}"').fetchall() for t in tables}
            self.assertEqual(migration.migrate(path, backup), backup.resolve())
            with closing(sqlite3.connect(str(path))) as db, closing(sqlite3.connect(str(backup))) as old, db, old:
                for table in tables:
                    self.assertEqual(old.execute(f'SELECT * FROM "{table}"').fetchall(), before[table])
                    rows = db.execute(f'SELECT * FROM "{table}"').fetchall()
                    if table == "experiments":
                        self.assertEqual([r[:-3] for r in rows], before[table])
                        self.assertTrue(all(r[-3:] == (5e-4, 4e-5, .025) for r in rows))
                    else:
                        self.assertEqual(rows, before[table])
                for value in (0, -1, None, "bad", float("inf")):
                    with self.assertRaises(sqlite3.IntegrityError):
                        db.execute("UPDATE experiments SET CHANNEL_WIDTH=?", (value,))
                db.execute("UPDATE experiments SET CHANNEL_WIDTH=.001 WHERE NAME='uniform'")
            self.assertIsNone(migration.migrate(path, backup))
            with closing(sqlite3.connect(str(path))) as db, db:
                self.assertEqual(db.execute("SELECT CHANNEL_WIDTH FROM experiments WHERE NAME='uniform'").fetchone()[0], .001)

    def test_partial_schema_rejected_without_writes(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "partial.db"
            with closing(sqlite3.connect(str(path))) as db, db:
                db.execute("CREATE TABLE experiments(NAME TEXT, CHANNEL_WIDTH REAL)")
            with self.assertRaises(ValueError):
                migration.migrate(path)
            with closing(sqlite3.connect(str(path))) as db, db:
                self.assertEqual(len(db.execute("PRAGMA table_info(experiments)").fetchall()), 2)

if __name__ == "__main__":
    unittest.main()
