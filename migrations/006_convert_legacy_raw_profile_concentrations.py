#!/usr/bin/env python3
"""Move legacy inlet matrices into raw_profile_entrance_concentrations.

Each line represents an inlet. The first three values map to SPECIES order;
an optional fourth value is discarded. Units come from the experiment's
SPECIE_INLET_CONC_UNITS list. Existing normalized rows are never overwritten.
"""
import math
import re
import sqlite3
import sys
from pathlib import Path
from typing import Tuple


def convert(database: Path) -> Tuple[int, int]:
    db = sqlite3.connect(str(database))
    try:
        db.execute("BEGIN IMMEDIATE")
        profiles = db.execute(
            "SELECT NAME, WT_PERCENT, ENTRANCE_CONC FROM raw_profile "
            "WHERE typeof(ENTRANCE_CONC)='text'"
        ).fetchall()
        insertions = []
        updates = []
        for name, weight, raw in profiles:
            if "#" in raw:
                continue
            lines = [line.strip() for line in re.split(r"\r\n|\r|\n", raw.strip()) if line.strip()]
            tokens = [re.split(r"\s+", line) for line in lines]
            if len(lines) <= 1 and (not tokens or len(tokens[0]) <= 1):
                continue
            experiment = db.execute(
                "SELECT SPECIES, SPECIE_INLET_CONC_UNITS, ENTRANCE_FLOWRATE "
                "FROM experiments WHERE NAME=?", (name,)
            ).fetchone()
            if not experiment:
                raise ValueError(f"No experiment row for raw profile {name!r}.")
            species = (experiment[0] or "").split()
            units = (experiment[1] or "").split()
            entrances = (experiment[2] or "").split()
            if len(species) != 3 or len(units) != len(species):
                raise ValueError(f"{name!r} needs three species and matching SPECIE_INLET_CONC_UNITS.")
            if len(lines) != len(entrances):
                raise ValueError(f"{name!r}/{weight}: {len(lines)} concentration lines for {len(entrances)} inlets.")
            if any(len(row) not in (3, 4) for row in tokens):
                raise ValueError(f"{name!r}/{weight}: each inlet line must have three or four values.")
            existing = db.execute(
                "SELECT count(*) FROM raw_profile_entrance_concentrations "
                "WHERE NAME=? AND WT_PERCENT=?", (name, str(weight))
            ).fetchone()[0]
            if existing:
                raise ValueError(f"{name!r}/{weight}: normalized concentration rows already exist.")
            normalized_lines = []
            for entrance_number, row in enumerate(tokens, start=1):
                normalized_lines.append(" ".join(row[:3]))
                for specie, unit, token in zip(species, units, row[:3]):
                    value = float(token)
                    if not math.isfinite(value) or value < 0:
                        raise ValueError(f"{name!r}/{weight}: concentrations must be finite and non-negative.")
                    insertions.append((name, str(weight), entrance_number, specie, value, unit))
            updates.append(("\n".join(normalized_lines), " ".join(units), name, weight))

        db.executemany(
            "INSERT INTO raw_profile_entrance_concentrations "
            "(NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS) "
            "VALUES (?, ?, ?, ?, ?, ?)", insertions
        )
        db.executemany(
            "UPDATE raw_profile SET ENTRANCE_CONC=?, ENTRANCE_CONC_UNITS=? "
            "WHERE NAME=? AND WT_PERCENT=?", updates
        )
        db.commit()
        return len(updates), len(insertions)
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {Path(sys.argv[0]).name} DATABASE")
    row_count, value_count = convert(Path(sys.argv[1]))
    print(f"Converted {row_count} raw profiles into {value_count} species/inlet values.")
