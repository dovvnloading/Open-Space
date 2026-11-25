"""After a successful full build, mutate a copy of each evidence table.
Check that Lean rejects the corresponding generated certificate, not merely
that the Python adapter rejects malformed data. Never modify deposited files.
"""
import csv
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import generate as g


def reject_mutation(work, name, attribute, selector, column, replacement, builder):
    path = getattr(g, attribute)
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        fields = reader.fieldnames
        rows = list(reader)
    row = next(row for row in rows if selector(row))
    old = row[column]
    row[column] = replacement(old)
    if row[column] == old:
        raise AssertionError(f"{name}: mutation changed nothing")
    mutated = work / f"{name}.csv"
    with mutated.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with patch.object(g, attribute, mutated):
        lines = builder()
    lean_file = work / f"{name}.lean"
    lean_file.write_text("\n".join([
        "import OpenSpace.Foundations", "set_option maxRecDepth 4096",
        "set_option maxHeartbeats 2000000", "namespace OpenSpace",
        *lines, "end OpenSpace", "",
    ]), encoding="utf-8")
    result = subprocess.run(
        ["lean", "-DwarningAsError=true", f"--root={work}", str(lean_file)],
        capture_output=True, text=True, timeout=60,
    )
    output = result.stdout + result.stderr
    if (result.returncode == 0 or "error:" not in output
            or "decide" not in output or "false" not in output):
        raise AssertionError(f"{name}: expected a false-proposition rejection, got:\n{output}")
    print(f"{name}: Lean rejected altered evidence ({column})")


def main(work):
    reject_mutation(
        work, "NegativeNICER", "NICER",
        lambda row: row["narrow_literal"] != "" and row["narrow_plus69"] != "",
        "narrow_literal", lambda x: str(int(x) + 1), g.nicer_lines,
    )
    reject_mutation(
        work, "NegativeOccultation", "CONTACTS",
        lambda row: row["date"] == "2026-12-24" and row["spacecraft"] == "SolO"
        and row["direction"] == "earth" and row["venus_height_km"] == "0"
        and row["source_radius_in_target_radii"] == "0",
        "ingress_utc", lambda x: x.replace("02.324Z", "02.325Z"), g.contact_lines,
    )
    reject_mutation(
        work, "NegativeOrbit", "ORBITS", lambda row: row["planet"] == "Kepler-445 b",
        "oec_reported_a_au", lambda _: "0.001", g.orbit_lines,
    )


if __name__ == "__main__":
    main(Path(sys.argv[1]).resolve())
