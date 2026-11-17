"""Reproduce the Euclid Q1 Discovery Engine and AgileLens catalogue screen (stdlib)."""
import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import Counter
from pathlib import Path
from urllib.request import urlopen

RELEASE = "https://zenodo.org/records/15025832/files/"
AGILE_URL = ("https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
             "-source=J/A%2bA/712/A207/tablee1&-out.all&-out.max=unlimited")
HASHES = {
    "q1_discovery_engine_lens_catalog.csv": "ee5e60cd507413eabf3da8ffa37a3c527212da00feb29083cf516a86d6f26877",
    "modeling_lens_mass.csv": "f2a52616a1ac65137b34abe0892d251e18e2c60e9a64c459914102f44ec49db5",
    "modeling_lens_mge.csv": "b85f7c78aebe898a3d7d77d07a72568cd2d596b525883b358b792e4ce9fa2432",
    "modeling_lens_sersic.csv": "c0b05783ecfdcb4068c3487db11c7eda18e070f8011726fc17d50a04aecd5dc2",
    "modeling_mge_magnitude.csv": "764f786e0cc003602258b724fa96c99f0a95dc084381dde36c5e37b2f98f937f",
    "modeling_sersic_magnitude.csv": "7c7506af33f27d65cc618e998cbbb0868c36a106703f303d6b0b0a74240d94f2",
    "agilelens.tsv": "ebca8eae5e70c9442c887fd15cd00e3a068abb6a779f66c1862be56a93c6c95e",
}
BANDS = {"vis": "vis", "y": "nir_y", "j": "nir_j", "h": "nir_h"}


def get_inputs(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for filename, expected in HASHES.items():
        path = directory / filename
        url = AGILE_URL if filename == "agilelens.tsv" else RELEASE + filename + "?download=1"
        if path.exists():
            data = path.read_bytes()
        else:
            with urlopen(url, timeout=120) as stream:
                data = stream.read()
        payload = data
        if filename == "agilelens.tsv":
            # The mirror adds the current query date to its generated comment header.
            payload = b"".join(line for line in data.splitlines(keepends=True)
                               if line.strip() and not line.startswith(b"#"))
        if hashlib.sha256(payload).hexdigest() != expected:
            raise ValueError(f"Input snapshot does not match SHA256 for {filename}")
        if not path.exists():
            path.write_bytes(data)


def model_table(directory, name):
    with (directory / (name + ".csv")).open() as stream:
        rows = list(csv.DictReader(stream))
    by_id = {r["id_str"]: r for r in rows}
    if len(by_id) != len(rows):
        raise ValueError(f"Duplicate identifier in {name}")
    return by_id


def number(row, key):
    return float(row.get(key) or "nan")


def output_csv(directory, filename, rows):
    with (directory / filename).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys(), lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: (f"{v:.9g}" if math.isfinite(v) else "") if isinstance(v, float) else v
                             for k, v in row.items()})


def annotate_ranks(rows, column):
    ordered = sorted((r for r in rows if math.isfinite(r[column])), key=lambda r: -r[column])
    ref = [r[column] for r in ordered]
    med = statistics.median(ref)
    mad = 1.4826 * statistics.median(abs(x - med) for x in ref)
    for row in rows:
        row[column + "_rank"] = ""
        row[column + "_robust_z"] = math.nan
    for rank, row in enumerate(ordered, 1):
        row[column + "_rank"] = rank
        row[column + "_robust_z"] = (row[column] - med) / mad
    return {"usable": len(ref), "median": med, "scaled_mad": mad, "maximum": max(ref),
            "top_id": ordered[0]["id_str"]}


def analyze(directory, out):
    candidates = model_table(directory, "q1_discovery_engine_lens_catalog")
    mass, light, sersic, phot, sphot = [
        model_table(directory, name) for name in
        ["modeling_lens_mass", "modeling_lens_mge", "modeling_lens_sersic",
         "modeling_mge_magnitude", "modeling_sersic_magnitude"]]
    if any(set(table) != set(mass) for table in [light, sersic, phot, sphot]):
        raise ValueError("Model table row IDs differ")
    selected = [r for r in candidates.values() if r["subset"] in {"discovery_engine", "q1_serendiptous"}
                and r["grade"] in {"A", "B"}]
    rows = []
    for c in selected:
        key = c["id_str"]
        row = {k: c[k] for k in ["id_str", "object_id", "subset", "grade", "right_ascension", "declination",
                                 "flux_vis_1fwhm_aper", "notes", "references"]}
        model = mass.get(key, {})
        theta = number(model, "einstein_radius_effective_median_pdf")
        row["model_state"] = ("no_row" if key not in mass else "missing_radius" if not math.isfinite(theta)
                              else "nonpositive_radius" if theta <= 0 else "positive_radius")
        for summary in ["median_pdf", "lower_1_sigma", "upper_1_sigma"]:
            row["theta_eff_arcsec_" + summary] = number(model, "einstein_radius_effective_" + summary)
        offsets = [number(model, f"centre_{i}_median_pdf") - number(light.get(key, {}), f"centre_{i}_median_pdf")
                   for i in range(2)]
        row["mass_light_offset_arcsec"] = math.hypot(*offsets)
        row["sersic_re_arcsec"] = number(sersic.get(key, {}), "effective_radius_median_pdf")
        row["sersic_n"] = number(sersic.get(key, {}), "sersic_index_median_pdf")
        for band, release_band in BANDS.items():
            field = release_band + "_lens_magnitude_ab_"
            row[band + "_mge_ab"] = number(phot.get(key, {}), field + "median_pdf")
            row[band + "_sersic_ab"] = number(sphot.get(key, {}), field + "median_pdf")
            for suffix in ["lower_1_sigma", "upper_1_sigma"]:
                row[f"{band}_mge_ab_{suffix}"] = number(phot.get(key, {}), field + suffix)
            row[band + "_priority"] = row[band + "_mge_ab"] + 5 * math.log10(theta) if theta > 0 else math.nan
        rows.append(row)
    summary = {"release_rows": len(candidates), "model_rows": len(mass), "accepted": len(rows),
               "model_grade_counts": dict(Counter(candidates.get(k, {}).get("grade", "unmatched") for k in mass)),
               "unmatched_model_ids": sorted(set(mass) - set(candidates)),
               "unassigned_A_B_ids": sorted(r["id_str"] for r in candidates.values()
                                            if not r["subset"] and r["grade"] in {"A", "B"}),
               "state_counts": dict(Counter(r["model_state"] for r in rows)), "bands": {}}
    for band in BANDS:
        values = [r[band + "_mge_ab"] for r in rows if math.isfinite(r[band + "_mge_ab"])]
        summary["bands"][band] = {"finite_magnitudes": len(values), "faintest_mge_mag": max(values),
                                  **annotate_ranks(rows, band + "_priority")}
    output_csv(out, "slde_screen.csv", rows)
    by_obj = {r["object_id"]: r for r in rows if r["object_id"]}
    with (directory / "agilelens.tsv").open() as stream:
        lines = [line for line in stream if line.strip() and not line.startswith("#")]
    independent = []
    for a in csv.DictReader([lines[0], *lines[3:]], delimiter="\t"):
        a = {k: v.strip() for k, v in a.items()}
        object_id = a["cutoutName"].split("_")[-1]
        ra, dec = float(a["RAJ2000"]), float(a["DEJ2000"])
        # Spherical distance; published cutout centres can be displaced from the lens by up to 10".
        distances = [(separation(ra, dec, float(r["right_ascension"]), float(r["declination"])), r["id_str"])
                     for r in rows]
        distance, near_id = min(distances)
        independent.append({"name": a["objName"], "object_id": object_id, "id_str": a["cutoutName"],
                            "ra_deg": ra, "dec_deg": dec, "grade_xu": a["Grade"], "new_xu": a["New"],
                            "vis_detection_total_mag": a["VISmag"], "y_2fwhm_mag": a["Ymag"],
                            "j_2fwhm_mag": a["Jmag"], "h_2fwhm_mag": a["Hmag"],
                            "slde_ab_same_object_id": by_obj.get(object_id, {}).get("id_str", ""),
                            "nearest_slde_ab": near_id, "nearest_slde_ab_arcsec": distance})
    if len({r["object_id"] for r in independent}) != len(independent):
        raise ValueError("Duplicate merged object in AgileLens")
    output_csv(out, "agilelens_screen.csv", independent)
    summary["agilelens"] = {"table_rows": len(independent), "marked_new": sum(r["new_xu"] == "Y" for r in independent),
                            "slde_exact_ab": sum(bool(r["slde_ab_same_object_id"]) for r in independent),
                            "slde_ab_within_2arcsec": sum(r["nearest_slde_ab_arcsec"] < 2 for r in independent),
                            "slde_ab_within_10arcsec": sum(r["nearest_slde_ab_arcsec"] < 10 for r in independent)}
    (out / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")


def separation(ra1, dec1, ra2, dec2):
    dra, ddec = map(math.radians, [ra2 - ra1, dec2 - dec1])
    a = math.sin(ddec/2)**2 + math.cos(math.radians(dec1)) * math.cos(math.radians(dec2)) * math.sin(dra/2)**2
    return math.degrees(2 * math.asin(math.sqrt(min(1, a)))) * 3600


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_dir", type=Path)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    get_inputs(args.input_dir)
    analyze(args.input_dir, args.output)
