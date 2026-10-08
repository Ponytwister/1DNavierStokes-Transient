"""Create a new, empty Navier experiment database with the application schema."""
import argparse
from contextlib import closing
from pathlib import Path
import sqlite3


def create_database(destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    # Exclusive creation ensures an existing experiment database is never replaced.
    with destination.open("xb"):
        pass
    try:
        with closing(sqlite3.connect(str(destination))) as database:
            with database:
                schema = Path(__file__).resolve().parents[1] / "schema/blank_database.sql"
                database.executescript(schema.read_text(encoding="utf-8"))
                database.execute("PRAGMA user_version = 1")
                if database.execute("PRAGMA integrity_check").fetchone()[0] != "ok":
                    raise RuntimeError("The created database failed SQLite integrity_check")
    except Exception:
        if destination.exists():
            destination.unlink()
        raise
    return destination


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("database", type=Path, help="new database path (must not already exist)")
    args = parser.parse_args()
    print(f"Created empty database: {create_database(args.database)}")
