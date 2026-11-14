"""Replay the Gaia WD+MS companion screen and the CD−27 2812 control measurements."""
import csv
import gzip
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"
FAINT = "2897976062821487360"
BRIGHT = "2897976170198264064"


def read_rows(filename):
    path = DATA / filename
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def angle(first, second):
    dec1, dec2 = (math.radians(float(row["dec"])) for row in (first, second))
    dra = math.radians(float(first["ra"]) - float(second["ra"]))
    hav = (math.sin((dec1 - dec2) / 2) ** 2
           + math.cos(dec1) * math.cos(dec2) * math.sin(dra / 2) ** 2)
    return math.degrees(2 * math.asin(math.sqrt(min(hav, 1)))) * 3600


def pm_difference(first, second):
    return math.hypot(*(float(first[k]) - float(second[k]) for k in ("pmra", "pmdec")))


def matching_pair(host, other):
    parallax = float(host["parallax"])
    separation = angle(host, other)
    if other["source_id"] == host["source_id"] or separation > min(300, max(20, 25 * parallax)):
        return False
    projected_au = separation * 1000 / parallax
    escape_speed = 42.12 * math.sqrt(2 / projected_au)
    relative_pm = pm_difference(host, other)
    sigma_pm = math.sqrt(sum(float(r[k]) ** 2 for r in (host, other)
                             for k in ("pmra_error", "pmdec_error")))
    delta_p = abs(parallax - float(other["parallax"]))
    sigma_p = math.hypot(*(float(r["parallax_error"]) for r in (host, other)))
    return (delta_p <= min(3 * sigma_p + 0.05, max(2 * sigma_p, 0.05))
            and sigma_p < 0.15 * parallax and float(other["ruwe"]) < 2.5
            and relative_pm <= 3 * sigma_pm + escape_speed * parallax / 4.74047
            and relative_pm * 4.74047 / parallax < 1.3 * escape_speed)


def main():
    hosts = sorted(read_rows("hosts.csv.gz"), key=lambda row: row["source_id"])
    cone_groups = defaultdict(list)
    for row in read_rows("cone-returns.csv.gz"):
        cone_groups[int(row["group"])].append(row)
    pairs = [(host, other) for index, host in enumerate(hosts)
             for other in cone_groups[index // 9] if matching_pair(host, other)]
    nss_pairs = [[a["source_id"], b["source_id"]] for a, b in pairs if int(b["non_single_star"])]
    records = {r["source_id"]: r for r in read_rows("gaia-pair.csv")}
    a, b = records[BRIGHT], records[FAINT]
    sigma = math.hypot(float(a["parallax_error"]), float(b["parallax_error"]))
    plx = float(a["parallax"])
    separation = angle(a, b)
    projected = separation * 1000 / plx
    delta_pm = pm_difference(a, b)
    field = read_rows("control-field.csv.gz")
    def control_ids(dp, dmu):
        return sorted(row["source_id"] for row in field
                      if abs(float(row["parallax"]) - plx) < dp and pm_difference(row, a) < dmu)
    comparisons = read_rows("comparison-colors.csv.gz")
    colors = [float(x["bp_rp"]) for x in comparisons]
    color = float(b["bp_rp"])
    return {
        "screen": {"hosts": len(hosts), "cone_rows": sum(map(len, cone_groups.values())),
                   "pairs": len(pairs), "pairs_with_gaia_nss_wide_source": nss_pairs},
        "outer_pair": {"angular_separation_arcsec": separation,
                       "reference_parallax_mas": plx, "projected_separation_au": projected,
                       "parallax_difference_mas": plx-float(b["parallax"]),
                       "parallax_difference_sigma": (plx-float(b["parallax"])) / sigma,
                       "relative_pm_mas_yr": delta_pm, "relative_vtan_km_s": 4.74047 * delta_pm / plx,
                       "escape_at_projected_sep_bright_pair_only_km_s": 42.12 * math.sqrt(1.9221 / projected)},
        "field_control": {"sources": len(field), "nominal_ids": control_ids(0.5, 1.0),
                          "wider_ids": control_ids(1.0, 2.0)},
        "color_control": {"absolute_rp_at_primary_distance_mag": float(b["phot_rp_mean_mag"]) + 5*math.log10(plx) - 10,
                          "target_bp_rp_mag": color, "nearby_sources": len(colors),
                          "median_bp_rp_mag": statistics.median(colors),
                          "count_as_blue_or_bluer": sum(x <= color for x in colors)}
    }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
