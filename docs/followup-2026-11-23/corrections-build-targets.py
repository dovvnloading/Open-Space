"""Build a complete, unranked NGC 1261b follow-up list from verified atlas bytes.

This is a selection manifest, not a new membership catalog. Gaia IDs stay strings.
Usage: python corrections-build-targets.py --atlas /path/to/ibata.txt
"""

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SOURCE_SHA256 = "76742daaf44d69e1afe71b56f49ec93c88631a6ae10261cc4adcb6f8485644ba"


def build(atlas, output):
    content = atlas.read_bytes()
    if hashlib.sha256(content).hexdigest() != SOURCE_SHA256:
        raise ValueError("Atlas bytes differ from the pinned primary source")
    rows = [line.split() for line in content.decode().splitlines()
            if line[:19].strip().isdigit() and len(line.split()) == 12]
    if len(rows) != 24540:
        raise ValueError("Unexpected atlas population")
    with (ROOT / "results/ngc1261b-reticulum-ii/matches.csv").open() as stream:
        matched = {row["gaia_source_id"] for row in csv.DictReader(stream)}
    coverage = json.loads((ROOT / "results/ngc1261b-reticulum-ii/external-coverage.json").read_text())
    later = {row["GaiaDR3"]: row for row in coverage["awad_b1_exact_matches"]}
    groups = defaultdict(set)
    for row in rows:
        groups[row[0]].add(int(row[-1]))
    branch = [row for row in rows if row[-1] == "14"]
    selected = [row for row in branch if row[0] not in matched]
    if len(branch) != 70 or len(matched) != 12 or len(selected) != 58:
        raise ValueError("Unexpected selection counts")
    if len({row[0] for row in selected}) != 58 or not set(later) <= {row[0] for row in selected}:
        raise ValueError("Duplicate IDs or missing later-spectrum target")
    records = []
    for row in sorted(selected, key=lambda row: float(row[1])):
        id_, ra, dec, plx, pmra, pmdec, g, color, rv, error, source, _ = row
        if int(source) > 0 and float(error) < 300:
            raise ValueError("Unexpected atlas velocity in selected sample")
        followup = later.get(id_, {})
        records.append(dict(
            gaia_source_id=id_, ra_deg=ra, dec_deg=dec,
            coordinate_epoch_as_published="J2016.5", parallax_mas=plx,
            pm_ra_cosdec_mas_per_yr=pmra, pm_dec_mas_per_yr=pmdec,
            gmag_extinction_corrected=g, bp_rp_extinction_corrected=color,
            atlas_stream_ids=";".join(map(str, sorted(groups[id_]))),
            atlas_velocity_available=0,
            known_later_spectrum="Awad2025_Table_B1" if followup else "not_in_checked_Awad_table",
            later_rv_kms=followup.get("RV", ""),
            later_rv_error_kms=followup.get("e_RV", ""),
            later_feh_dex=followup.get("[Fe/H]", ""),
            later_feh_error_dex=followup.get("e_[Fe/H]", ""),
            later_ngc1261_membership_probability=followup.get("Pmemb", ""),
        ))
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
    print(f"Wrote {len(records)} candidates; {len(later)} have the checked later Awad spectrum.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=HERE / "corrections-ngc1261b-targets.csv")
    args = parser.parse_args()
    build(args.atlas, args.output)
