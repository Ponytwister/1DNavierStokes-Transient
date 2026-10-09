#!/usr/bin/env python3
"""Move remaining raw_profile concentration values to per-inlet rows.

Existing normalized values win on key conflicts. After a successful conversion,
the legacy concentration columns are cleared and PRAGMA user_version is set to 4.
"""
import math
import re
import sqlite3
import sys
from pathlib import Path
from typing import Any, List, Optional, Tuple


def parse_value(token: str, profile: str) -> Tuple[Any, bool]:
    token = token.strip()
    if not token:
        raise ValueError(f"{profile}: empty concentration value.")
    if "#" in token:
        number, alias = token.split("#", 1)
        if not alias or re.search(r"\s|#", alias):
            raise ValueError(f"{profile}: invalid concentration alias {token!r}.")
        value: Optional[float] = None
        if number:
            try:
                value = float(number)
            except ValueError as error:
                raise ValueError(f"{profile}: invalid concentration value {token!r}.") from error
            if not math.isfinite(value) or value < 0:
                raise ValueError(f"{profile}: concentration values must be finite and non-negative.")
        return token, value is not None
    try:
        value = float(token)
    except ValueError as error:
        raise ValueError(f"{profile}: invalid concentration value {token!r}.") from error
    if not math.isfinite(value) or value < 0:
        raise ValueError(f"{profile}: concentration values must be finite and non-negative.")
    return value, True


def as_tokens(value: Any) -> List[str]:
    return re.split(r"\s+", str(value).strip()) if str(value).strip() else []


def convert(database: Path) -> Tuple[int, int]:
    db = sqlite3.connect(str(database))
    inserted = 0
    cleared = 0
    try:
        db.execute("BEGIN IMMEDIATE")
        rows = db.execute(
            "SELECT NAME, CAST(WT_PERCENT AS TEXT), ENTRANCE_CONC "
            "FROM raw_profile WHERE ENTRANCE_CONC IS NOT NULL"
        ).fetchall()
        for name, weight, raw in rows:
            profile = f"{name!r}/{weight}"
            experiment = db.execute(
                "SELECT SPECIES, SPECIE_INLET_CONC_UNITS, ENTRANCE_FLOWRATE "
                "FROM experiments WHERE NAME=?", (name,)
            ).fetchone()
            if not experiment:
                raise ValueError(f"{profile}: no matching experiment row.")
            species = as_tokens(experiment[0] or "")
            inlet_units = as_tokens(experiment[1] or "")
            entrances = as_tokens(experiment[2] or "")
            if not species or len(inlet_units) != len(species):
                raise ValueError(f"{profile}: species and SPECIE_INLET_CONC_UNITS must have matching entries.")

            if isinstance(raw, (int, float)):
                scalar = True
                matrix = []
                scalar_token = str(raw)
            else:
                text = str(raw).strip()
                if not text:
                    scalar = False
                    matrix = []
                    scalar_token = ""
                else:
                    lines = [line.split() for line in re.split(r"\r\n|\r|\n", text) if line.strip()]
                    scalar = len(lines) == 1 and len(lines[0]) == 1
                    scalar_token = lines[0][0] if scalar else ""
                    matrix = [] if scalar else lines

            if scalar:
                if "FITC" not in species:
                    raise ValueError(f"{profile}: scalar ENTRANCE_CONC requires FITC in SPECIES.")
                fitc_index = species.index("FITC")
                token_value, _ = parse_value(scalar_token, profile)
                values = [(1, "FITC", token_value, inlet_units[fitc_index])]
            elif matrix:
                if len(matrix) != len(entrances):
                    raise ValueError(f"{profile}: {len(matrix)} concentration rows for {len(entrances)} entrances.")
                if any(len(line) not in (len(species), len(species) + 1) for line in matrix):
                    raise ValueError(f"{profile}: each concentration row needs {len(species)} values; an optional fourth is ignored.")
                values = []
                for entrance_number, line in enumerate(matrix, start=1):
                    for index, species_name in enumerate(species):
                        stored_value, _ = parse_value(line[index], profile)
                        values.append((entrance_number, species_name, stored_value, inlet_units[index]))
            else:
                values = []

            for entrance_number, species_name, concentration, unit in values:
                cursor = db.execute(
                    "INSERT OR IGNORE INTO raw_profile_entrance_concentrations "
                    "(NAME, WT_PERCENT, ENTRANCE_NUMBER, SPECIES_NAME, CONCENTRATION, UNITS) "
                    "VALUES (?, ?, ?, ?, ?, ?)",
                    (name, weight, entrance_number, species_name, concentration, unit),
                )
                inserted += cursor.rowcount

        cursor = db.execute(
            "UPDATE raw_profile SET ENTRANCE_CONC=NULL, ENTRANCE_CONC_UNITS=NULL "
            "WHERE ENTRANCE_CONC IS NOT NULL OR ENTRANCE_CONC_UNITS IS NOT NULL"
        )
        cleared = cursor.rowcount
        db.execute("PRAGMA user_version=4")
        db.commit()
        return inserted, cleared
    except Exception:
        db.rollback()
        raise
    finally:
        db.close()


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit(f"Usage: {Path(sys.argv[0]).name} DATABASE")
    inserted, cleared = convert(Path(sys.argv[1]))
    print(f"Moved {inserted} concentrations; cleared legacy values on {cleared} raw profiles.")
