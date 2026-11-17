"""Replay fixed selections on published, directly matched Gaia DR3 controls."""
import csv
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from zero_point import zpt

HERE = Path(__file__).resolve().parent
K = 4.74047
zpt.load_tables()


def keyed(filename):
    with (HERE / filename).open() as source:
        rows = list(csv.DictReader(source))
    out = {r["source_id"]: r for r in rows}
    assert len(out) == len(rows)
    return out


def f(row, name):
    return float(row[name]) if row[name] else np.nan


def passes_quality(row):
    return (f(row, "ruwe") < 1.4 and int(row["visibility_periods_used"]) >= 9
            and row["duplicated_source"].lower() == "false"
            and int(row["ipd_frac_multi_peak"]) <= 2
            and f(row, "ipd_gof_harmonic_amplitude") < 0.1)


controls = keyed("controls.csv")
gaia = keyed("gaia.csv")
bj = keyed("bj.csv")
published = keyed("published_targets.csv")
assert controls.keys() == gaia.keys() == bj.keys()

report = []
for ident, control in controls.items():
    row, distances = gaia[ident], bj[ident]
    parallax = f(row, "parallax")
    sigma = f(row, "parallax_error")
    pm = math.hypot(f(row, "pmra"), f(row, "pmdec"))
    ruwe = f(row, "ruwe")
    colour = f(row, "bp_rp")
    gmag = f(row, "phot_g_mean_mag")
    raw_abs = gmag + 5 * math.log10(parallax) - 10 if parallax > 0 else np.nan
    eb26 = (gmag < 20.5 and colour < 0.5 and pm > 50 and ruwe < 1.4
            and (parallax < 0 or (4.74 * pm > 600 * parallax
                 and 4.74 * pm > 400 * (parallax + sigma))))
    assert eb26 == (ident in published), (ident, "literature mismatch")
    stage1 = pm > 35 and colour < 1.2
    stage2 = stage1 and parallax / sigma > 3
    stage3 = stage2 and K * pm > 1000 * parallax
    fast = stage3 and raw_abs > 3.5
    quality = passes_quality(row)
    # Vector inputs are supported by both NumPy 1.x and 2.x. Invalid calibrations
    # return NaN instead of extending the colour relation beyond its domain.
    zp = zpt.get_zpt([gmag], [f(row, "nu_eff_used_in_astrometry")], [f(row, "pseudocolour")],
                     [f(row, "ecl_lat")], [int(row["astrometric_params_solved"])], _warnings=False)[0]
    corrected = parallax - zp
    corrected_v = K * pm / corrected if corrected > 0 else np.nan
    record = dict(**control, pm_masyr=round(pm, 3), raw_parallax_snr=round(parallax/sigma, 3),
                  zp_mas=round(zp, 6) if math.isfinite(zp) else "",
                  corrected_parallax_snr=round(corrected/sigma, 3) if math.isfinite(zp) else "",
                  corrected_vtan_kms=round(corrected_v, 1) if math.isfinite(corrected_v) else "",
                  bj_geo_vtan_med_kms=round(K*pm*f(distances,"r_med_geo")/1000, 1),
                  bj_geo_vtan_lo_kms=round(K*pm*f(distances,"r_lo_geo")/1000, 1),
                  bj_geo_vtan_hi_kms=round(K*pm*f(distances,"r_hi_geo")/1000, 1),
                  pass_eb26=eb26, pass_pm_colour=stage1, pass_snr=stage2,
                  pass_raw_vtan=stage3, pass_fast=fast, pass_quality=quality,
                  pass_fast_quality=fast and quality,
                  pass_corrected_5sigma_v1000=(quality and corrected/sigma >= 5 and corrected_v > 1000),
                  has_xp=row["has_xp_continuous"] == "true")
    report.append(record)

with (HERE / "recovery.csv").open("w", newline="") as out:
    writer = csv.DictWriter(out, fieldnames=report[0].keys(), lineterminator="\n")
    writer.writeheader()
    writer.writerows(report)
counts = {}
for cls in ["D6", "LP40", "He", "All"]:
    group = [r for r in report if cls == "All" or cls == r["class"]]
    counts[cls] = {"controls": len(group)}
    for key in report[0]:
        if key.startswith("pass_") or key == "has_xp":
            counts[cls][key] = sum(bool(r[key]) for r in group)
summary = {"control_recovery": counts, "published_92_target_verdicts":
           dict(sorted(Counter(r["verdict"] for r in published.values()).items()))}
# The label is visual/observational evidence for blending, not a known true parallax.
validation = keyed("validation_gaia.csv")
fidelity = keyed("fidelity.csv")
assert validation.keys() == published.keys()
assert set(fidelity) == set(validation) | set(gaia)
strata = {"published_blends": [], "published_other": [], "positive_controls": []}
for ident in fidelity:
    row = validation.get(ident, gaia.get(ident))
    frow = fidelity[ident]
    # Blank norm_dg means no reported neighbour statistic, which is left unflagged.
    flags = {
        "fails_gaia_quality": not passes_quality(row),
        "fails_fidelity": f(frow, "fidelity_v2") <= 0.5,
        "neighbour_flag": f(frow, "norm_dg") >= -3,
    }
    flags["any_flag"] = any(flags.values())
    if ident in published:
        key = "published_blends" if published[ident]["verdict"] == "blended" else "published_other"
        strata[key].append(flags)
    if ident in controls:
        strata["positive_controls"].append(flags)
summary["advisory_flags_by_verdict"] = {
    key: dict(total=len(flags), **{field: sum(r[field] for r in flags) for field in flags[0]})
    for key, flags in strata.items()
}
print(json.dumps(summary, indent=2))
(HERE / "recovery_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
