"""Apply migration 001 to an existing database, with a consistent SQLite backup."""
import argparse
from contextlib import closing
from datetime import datetime, timezone
from pathlib import Path
import sqlite3

COLUMNS = {"CHANNEL_WIDTH", "CHANNEL_HEIGHT", "CHANNEL_LENGTH"}

def migrate(database, backup=None):
    database = Path(database).resolve(strict=True)
    with closing(sqlite3.connect(database.as_uri() + "?mode=rw", uri=True)) as db, db:
        columns = {row[1].upper() for row in db.execute("PRAGMA table_info(experiments)")}
        if not columns:
            raise ValueError("The database has no experiments table")
        present = columns & COLUMNS
        if present == COLUMNS:
            return None  # Never reset previously edited dimensions.
        if present:
            raise ValueError("Partial channel schema; restore the backup or resolve it before migrating")
        backup = Path(backup) if backup else database.with_name(
            database.name + ".before-channel-dimensions-" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ") + ".bak")
        backup = backup.resolve()
        # Exclusive creation prevents overwriting an existing backup or database.
        with backup.open("xb"):
            pass
        with closing(sqlite3.connect(str(backup))) as target, target:
            db.backup(target)
        migration = Path(__file__).resolve().parents[1] / "migrations/001_channel_dimensions.sql"
        try:
            db.executescript(migration.read_text(encoding="utf-8"))
        except Exception:
            db.rollback()
            raise
        if db.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
            raise RuntimeError("Integrity check failed; retain and restore the backup")
        return backup

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path)
    parser.add_argument("--backup", type=Path)
    args = parser.parse_args()
    result = migrate(args.database, args.backup)
    print(f"Migration applied. Backup: {result}" if result else "Already migrated; existing values preserved.")
