"""Compare published orbital clocks with previously fitted native-cadence alias models."""
import argparse
import gzip
import json
import math
from pathlib import Path


# Printed reference values are rounded. Each source evolves on its own ephemeris.
# The candidate key below only selects the intended daily-survey record.
FALANGA = "Falanga et al. (2015)", "https://ui.adsabs.harvard.edu/abs/2015A%26A...577A.130F/abstract"
HER = "Staubert, Klochkov & Wilms (2009)", "https://arxiv.org/abs/0904.2307"
HEMPHILL = "Hemphill et al. (2019)", "https://ui.adsabs.harvard.edu/abs/2019ApJ...873...62H/abstract"
PRIORS = [
    ("H1538-522", "H 1538−522", 2.1572, 3.72831, 2, HEMPHILL),
    ("SMCX-1", "SMC X-1", 2.05723, 3.891923160, 2, FALANGA),
    ("HerX-1", "Her X-1", 2.42823, 1.7001676, 1, HER),
    ("LMCX-4", "LMC X-4", 3.44882, 1.4083790, 1, FALANGA),
]


def read_lines(path):
    if path.exists():
        return [json.loads(s) for s in path.read_text().splitlines() if s.strip()]
    return []


def compare(details, period, prior, harmonic):
    carrier_centre = abs(1 - 1/period)
    target_frequency = 1/period
    grids = details.get("fits", [])
    if details.get("status") != "ok":
        return {"status":details.get("status", "no native fit")}
    def get(centre):
        return next((x for x in grids if abs(x["predicted_frequency"]-centre) < 1e-11), None)
    carrier, target = get(carrier_centre), get(target_frequency)
    orbit = get(carrier_centre/harmonic)
    if not all((carrier, target, orbit)):
        return {"status":"expected daily, carrier or orbital hypothesis absent"}
    span = details["baseline_days"]
    # Basis degeneracy was stored against the strongest candidate alias only.
    # Do not transfer an unrelated alternative's window diagnostics to this one.
    strongest = details["best_one_day_alias"]
    same = math.isclose(strongest["predicted_frequency"], carrier_centre, abs_tol=1e-11, rel_tol=0.)
    is_orbit = math.isclose(strongest["predicted_frequency"], carrier_centre/harmonic, abs_tol=1e-11, rel_tol=0.)
    return {
        "status":"compared", "native_rows":details["n"], "span_days":span,
        "daily_fit_period_days":target["period_days"], "daily_gain":target["delta_chi2"],
        "carrier_fit_period_days":carrier["period_days"], "carrier_gain":carrier["delta_chi2"],
        "carrier_minus_daily":carrier["delta_chi2"] - target["delta_chi2"],
        "orbital_fit_period_days":orbit["period_days"], "orbital_gain":orbit["delta_chi2"],
        "carrier_centre_offset_from_prior_rayleigh":abs(carrier_centre - harmonic/prior)*span,
        "carrier_fit_offset_from_prior_rayleigh":abs(carrier["frequency"] - harmonic/prior)*span,
        "strongest_alias_period_days":strongest["period_days"], "carrier_is_strongest_alias":same,
        "carrier_window":details.get("winner_window") if same else None,
        "orbital_window":details.get("winner_window") if is_orbit else None,
        "runtime_seconds":details.get("runtime_seconds")
    }


def build(evidence):
    orbit_dir = evidence/"review"/"orbit"
    root = json.loads((orbit_dir/"summary.json").read_text()) if (orbit_dir/"summary.json").exists() else {}
    residuals = read_lines(evidence/"review"/"residuals.jsonl")
    outputs = []
    for ident, name, key, orbit, harmonic, reference in PRIORS:
        result = dict(id=ident, source=name, daily_record_selector_days=key,
                      published_orbit_days=orbit, harmonic=harmonic, reference=reference,
                      carrier_prior_days=orbit/harmonic,
                      expected_daily_alias_days=1/abs(1-harmonic/orbit))
        outputs.append(result)
        # Broad tolerance here is only for selecting the requested table row.
        matches = [r for r in root.get("tests", []) if r["source_id"] == ident and abs(r["period_days"]-key)<.002]
        if len(matches) != 1:
            result["status"] = "native follow-up absent or daily record ambiguous"
            continue
        row = matches[0]
        period = row["period_days"]
        result.update(status="native comparisons", daily_days=period)
        path = orbit_dir/row["all_predicted_fits"]
        if not path.exists():
            result["status"] = "native comparison file unavailable"
            continue
        with gzip.open(path, "rt") as handle:
            fits = json.load(handle)
        result["instruments"] = {
            band:{epoch:compare(fits.get(band, {}).get(epoch, {}), period, orbit, harmonic)
                  for epoch in ("full", "early", "late")} for band in ("bat", "maxi")}
        comparisons = [s for epochs in result["instruments"].values() for s in epochs.values() if s["status"]=="compared"]
        result["comparisons"] = len(comparisons)
        result["carrier_preferred_by_fit"] = sum(s["carrier_minus_daily"]>0 for s in comparisons)
        result["centres_within_one_rayleigh_of_prior"] = sum(
            s["carrier_centre_offset_from_prior_rayleigh"]<=1 for s in comparisons)
        result["fitted_carriers_within_one_rayleigh_of_prior"] = sum(
            s["carrier_fit_offset_from_prior_rayleigh"]<=1 for s in comparisons)
        result["daily_residuals"] = {}
        residual = next((r for r in residuals if r["source_id"] == ident and abs(r["period_days"]-period)<1e-10), {})
        for band in ("bat","maxi"):
            v = residual.get("instruments", {}).get(band, {})
            result["daily_residuals"][band] = {"status":v.get("status", "missing"), "flags":v.get("flags", []),
                                               "after_conditioning":v.get("grid_after", {})}
    return outputs


def write(outputs, dest):
    text = ["# Published orbital periods and one-day aliases", "",
            "For each daily-survey period, the proposed carrier is the published orbital frequency or its second harmonic. Once-per-day sampling folds this carrier onto the lower frequency $|1 - f_{\\rm carrier}|$. The orbit-product fits use their native sampling, quoted flux errors and piecewise means. The daily-frequency and carrier gains reduce the same weighted chi-square, with identical local frequency refinement. A positive difference favors the carrier in that interval. These are descriptive fits to the same underlying monitor observations used in selection.", "",
            "Offsets use the epoch's Rayleigh width, 1 divided by its span. The frequency predicted from the printed orbital period can evolve across the baseline. An offset larger than one resolution element is reported explicitly. The output JSON includes conditional gains and sampling correlations when the carrier or its orbital fundamental was the stored winner. Near-unit correlations make the corresponding sinusoidal components difficult to distinguish.", ""]
    for r in outputs:
        name = r["source"]
        ref, url = r["reference"]
        text.extend([f"## {name}", "", f"The published orbit is {r['published_orbit_days']:.8g} d. Harmonic {r['harmonic']} has a period of {r['carrier_prior_days']:.8g} d and predicts a daily reflection at {r['expected_daily_alias_days']:.8g} d. [{ref}]({url})", ""])
        if r["status"] != "native comparisons":
            text.extend([r["status"] + ".", ""])
            continue
        text.extend([f"The measured daily period is {r['daily_days']:.8g} d. The reflected carrier has a larger gain in {r['carrier_preferred_by_fit']} of {r['comparisons']} available comparisons. Its fitted maximum is within one Rayleigh width of the printed orbital prediction in {r['fitted_carriers_within_one_rayleigh_of_prior']} comparisons.", "",
                     "| Monitor | Epoch | Daily gain | Carrier period (d) | Carrier gain | Carrier − daily | Orbital gain | Carrier offset (Rayleigh) |",
                     "|---|---|---:|---:|---:|---:|---:|---:|"])
        for band, epochs in r["instruments"].items():
            for epoch, s in epochs.items():
                if s["status"] != "compared":
                    text.append(f"| {band.upper()} | {epoch} | {s['status']} | — | — | — | — | — |")
                else:
                    text.append(f"| {band.upper()} | {epoch} | {s['daily_gain']:.5g} | {s['carrier_fit_period_days']:.8g} | {s['carrier_gain']:.5g} | {s['carrier_minus_daily']:+.5g} | {s['orbital_gain']:.5g} | {s['carrier_fit_offset_from_prior_rayleigh']:.2f} |")
        text.append("")
        for band, epochs in r["instruments"].items():
            s = epochs["full"]
            if s["status"] != "compared":
                continue
            if s["carrier_window"]:
                w = s["carrier_window"]
                c = w["canonical_correlations"]
                text.append(f"{band.upper()} full-span target/carrier sampling correlations are {c[0]:.3f} and {c[1]:.3f}. Adding the carrier to the daily fit gains {w['alias_gain_given_target']:.5g}, and adding the daily term to the carrier fit gains {w['target_gain_given_alias']:.5g}.")
            elif s["orbital_window"]:
                w = s["orbital_window"]
                c = w["canonical_correlations"]
                text.append(f"{band.upper()} favors the orbital fundamental at {s['orbital_fit_period_days']:.8g} d over its reflected harmonic and the daily period. The daily/fundamental sampling correlations are {c[0]:.3f} and {c[1]:.3f}. Adding the orbital fundamental to the daily fit gains {w['alias_gain_given_target']:.5g}, and adding the daily term to the orbital fit gains {w['target_gain_given_alias']:.5g}. The specific harmonic's conditional gains were not stored.")
            else:
                text.append(f"{band.upper()} has a stronger fitted alias at {s['strongest_alias_period_days']:.8g} d. Conditional and window comparisons for the specific published carrier were not stored.")
            if s["carrier_centre_offset_from_prior_rayleigh"]>1:
                text.append(f"The {band.upper()} predicted carrier centre differs from the orbital prior by {s['carrier_centre_offset_from_prior_rayleigh']:.2f} full-span Rayleigh widths.")
            v = r["daily_residuals"][band]
            text.append(f"The daily-data nuisance-fit status is “{v['status']}”" + (f" ({', '.join(v['flags'])})." if v["flags"] else "."))
            text.append("")
    dest.parent.mkdir(exist_ok=True, parents=True)
    dest.write_text("\n".join(text).rstrip()+"\n")
    dest.with_suffix(".json").write_text(json.dumps(outputs, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("output", type=Path, help="Markdown filename, or folder for prior-aliases.md and JSON")
    args = p.parse_args()
    destination = args.output if args.output.suffix.lower()==".md" else args.output/"prior-aliases.md"
    write(build(args.evidence), destination)
