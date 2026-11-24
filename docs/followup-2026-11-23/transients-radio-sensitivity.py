"""Recalculate UGC 2653's radio contrast under stated calibration scenarios.

This is uncertainty sensitivity of frozen peak measurements, not raw-image
recalibration, a formal detection probability, or a new observing epoch.
"""
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
SOURCE = REPO / "results/neowise-nuclear-flares/data/radio.json"
SPHEREX = REPO / "results/neowise-nuclear-flares/data/spherex-colours.json"


def calculate():
    target = next(row for row in json.loads(SOURCE.read_text()) if row["id"] == "03183820+3736290")
    early = target["vlass"][:2]
    late = target["vlass"][2]
    ivar = [1 / row["rms_mjy_beam"]**2 for row in early]
    weights = [value / sum(ivar) for value in ivar]
    baseline = sum(w * row["peak_mjy_beam"] for w, row in zip(weights, early))
    delta = late["peak_mjy_beam"] - baseline
    variance = late["rms_mjy_beam"]**2 + 1 / sum(ivar)
    # Same statistical weights in each scenario. If fractional gain errors
    # have common pairwise correlation rho, propagate all three covariances.
    independent_term = late["peak_mjy_beam"]**2 + sum(
        w**2 * row["peak_mjy_beam"]**2 for w, row in zip(weights, early))
    scenarios = []
    for fraction, rho in [(0., 0.), (.10, 0.), (.10, .5), (.10, 1.)]:
        error = math.sqrt(variance + fraction**2 * (
            (1 - rho) * independent_term + rho * delta**2))
        scenarios.append({"assumed_fractional_calibration_error": fraction,
                          "assumed_common_correlation": rho,
                          "error_mjy_beam": error, "contrast_over_error": delta / error})
    spherex = json.loads(SPHEREX.read_text())
    s = next(row for row in spherex["targets"] if row["id"] == target["id"])
    return {"base_commit": "972cd5169eea2b83d7407df84eb60c167b3c084e",
            "sources": [{"path": str(path.relative_to(REPO)),
                         "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in [SOURCE, SPHEREX]],
            "observing_dates": [row["date_obs"] for row in target["vlass"]],
            "baseline_mjy_beam": baseline, "statistical_weights": weights,
            "delta_mjy_beam": delta, "scenarios": scenarios,
            "spherex": {"colour_excess_mag": s["excess_mag_allwise"],
                        "control_scatter_mag": spherex["control_scatter_mag"],
                        "excess_over_scatter": s["excess_mag_allwise"] / spherex["control_scatter_mag"],
                        "controls": spherex["controls"],
                        "controls_at_least_as_red": s["controls_at_least_as_red"]},
            "limits": "10% gain errors and their correlations are assumptions, not measured calibration errors. The original measurements select the maximum pixel within 2.5 arcsec. Peak-selection bias, PSF and registration are not calibrated here. SPHEREx excess/scatter is not a calibrated Gaussian false-alarm significance; 65 comparison galaxies cannot measure a nine-sigma tail."}


if __name__ == "__main__":
    result = calculate()
    (HERE / "transients-radio-sensitivity.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
