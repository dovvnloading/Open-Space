"""Verify the frozen evidence and cross-check the reported populations without network access."""

import csv
import hashlib
import json
import math
from pathlib import Path


def records(name):
    with Path(name).open(newline="") as f:
        return list(csv.DictReader(f))


def main():
    paths = []
    for line in Path("SHA256SUMS").read_text().splitlines():
        digest, name = line.split("  ", 1)
        p = Path(name)
        with p.open("rb") as f:
            assert hashlib.file_digest(f, "sha256").hexdigest() == digest, name
        paths.append(name)
    depth = json.loads(Path("desi_depth.json").read_text())
    assert depth["input_rows"] == 4140838 and depth["specdis_speed_ge550_rows"] == 2209
    assert len(records("desi_screen.csv")) == 1020
    assert len(records("moderate/desi_screen.csv")) == 1187
    assert len({r["source_id"] for r in records("literature/matches.csv")}) == 10
    for name in ["inference.csv", "moderate/inference.csv"]:
        table = records(name)
        assert len(table) == len({(r["source_id"], r["model"]) for r in table})
        for row in table:
            for quantity in (
                "x",
                "y",
                "z",
                "vx",
                "vy",
                "vz",
                "d",
                "v_gal",
                "vesc_mc17",
                "vesc_gala22",
            ):
                values = [float(row[quantity + "_" + q]) for q in ("16", "50", "84")]
                assert all(math.isfinite(v) for v in values) and values == sorted(
                    values
                ), row["source_id"]
            assert all(
                0 <= float(row[k]) <= 1 for k in ("p_mc17", "p_gala22", "p_v1000")
            )
    lead = [
        r
        for r in records("candidate_models.csv")
        if r["source_id"] == "1717410430709663616"
    ]
    assert len(lead) == 9
    assert all(float(r["p_mc17"]) == float(r["p_gala22"]) == 1 for r in lead)
    ref = next(r for r in lead if r["model"] == "bj_photogeo")
    assert 696 < float(ref["v_gal_50"]) < 697
    mass = next(
        r
        for r in records("heavy_halo.csv")
        if r["source_id"] == ref["source_id"] and r["model"] == ref["model"]
    )
    assert float(mass["p_unbound"]) == 0
    assert (
        len(records("candidates.csv")) == 8
        and len(records("desi/los/line_fits.csv")) == 5
    )
    o = records("positive/orbit.csv")[0]
    assert 4.6 < float(o["first_plane_t_myr_50"]) < 4.9 and float(o["p_gc_1kpc"]) > 0.9
    print(
        f"Verified {len(paths)} frozen files, DESI selections, distance-model tables, literature matches and S5-HVS1 recovery."
    )


if __name__ == "__main__":
    main()
