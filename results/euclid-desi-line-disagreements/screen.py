"""Join positions to a frozen Euclid SPE catalogue and report redshift disagreements."""

from argparse import ArgumentParser
import hashlib
import json
from pathlib import Path

from astropy.io import fits
import numpy as np


def fit_record(row):
    fields = {
        "rank": "SPE_RANK", "z": "SPE_Z", "zerr": "SPE_Z_ERR",
        "probability": "SPE_Z_PROB", "reliability": "SPE_Z_REL",
        "reduced_chi2": "SPE_Z_REDUCED_CHI2",
    }
    record = {}
    for label, name in fields.items():
        value = float(row[name])
        record[label] = value if np.isfinite(value) else None
    record["rank"] = int(record["rank"])
    return record


def quality_record(row):
    return {
        "grism": str(row["SPE_GRISM"]).strip(),
        "valid_pixels": int(row["SPE_NPIX"]),
        "median_dithers": int(row["SPE_N_DITH_MED"]),
        "warning_flag": int(row["SPE_WARNING_FLAG"]),
        "error_flag": int(row["SPE_ERROR_FLAG"]),
        "wave_min": float(row["SPE_W_MIN"]),
        "wave_max": float(row["SPE_W_MAX"]),
    }


def screen(associations, spectra):
    fits_by_id, quality_by_id = {}, {}
    targets = {row["euclid_match"] for row in associations["matches"] if row["euclid_match"] is not None}
    candidate_records = spectra["SPE_GALAXY_CANDIDATES"].data
    candidate_records = candidate_records[np.isin(candidate_records["OBJECT_ID"], list(targets))]
    for row in candidate_records:
        fits_by_id.setdefault(int(row["OBJECT_ID"]), []).append(fit_record(row))
    qa_records = spectra["SPE_QUALITY"].data
    qa_records = qa_records[np.isin(qa_records["OBJECT_ID"], list(targets))]
    for row in qa_records:
        quality_by_id.setdefault(int(row["OBJECT_ID"]), []).append(quality_record(row))

    ranked = []
    for row in associations["matches"]:
        ident = row["euclid_match"]
        if ident is None or not row["desi_quality_gate"]:
            continue
        fits_for_object = sorted(fits_by_id.get(ident, []), key=lambda item: item["rank"])
        best = [solution for solution in fits_for_object if solution["rank"] == 0]
        item = {
            "euclid_id": ident, "desi_targetid": row["desi_targetid"],
            "desi_z": row["desi_z"], "separation_arcsec": row["neighbors"][0]["sep_arcsec"],
            "neighbors": row["neighbors"], "solutions": fits_for_object,
            "quality": quality_by_id.get(ident, []),
            "status": "no-unique-rank-zero", "signed_dz": None,
            "selected": False, "discordant": False,
        }
        if len(best) == 1:
            fit = best[0]
            z = fit["z"]
            item["status"] = "rank-zero-pair"
            if z is not None:
                delta = (z - row["desi_z"]) / (1 + row["desi_z"])
                eligible = (fit["probability"] is not None and fit["probability"] >= 0.99
                            and fit["reliability"] is not None and 0 <= fit["reliability"] < 0.2)
                item.update(signed_dz=delta, selected=eligible, discordant=eligible and abs(delta) > 0.01)
        ranked.append(item)
    return {
        "mer_unique": associations["mer_unique"], "desi_rows": associations["desi_rows"],
        "with_neighborhood": sum(bool(row["neighbors"]) for row in associations["matches"]),
        "provisional_associations": len(targets), "desi_quality_pairs": len(ranked),
        "unique_rank_zero": sum(row["status"] == "rank-zero-pair" for row in ranked),
        "selected": sum(row["selected"] for row in ranked),
        "discordant": sum(row["discordant"] for row in ranked),
        "pairs": ranked,
    }


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("associations", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_bytes())
    item = manifest["spectrum"]
    path = args.manifest.parent / item["file"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
        raise ValueError(f"SPE digest mismatch for {path.name}")
    with fits.open(path, memmap=True) as spec:
        result = screen(json.loads(args.associations.read_bytes()), spec)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(result, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
