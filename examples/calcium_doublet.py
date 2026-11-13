"""Reproduce the first three-galaxy investigation using fetched observation manifests."""
import argparse
import json
from pathlib import Path
from statistics import mean

from open_space.dataio import load_observation
from open_space.followup import measure
from open_space.spectra import analyze


def vacuum(air):
    """Invert the SDSS DR17 optical Ciddor formula by fixed-point iteration."""
    wave = air
    for _ in range(4):
        sigma2 = (10000 / wave) ** 2
        wave = air * (1 + .05792105 / (238.0185 - sigma2) + .00167917 / (57.362 - sigma2))
    return wave


K, H = vacuum(3933.6614), vacuum(3968.4673)
H_EPSILON = 3971.1951


def investigate(path):
    manifest, rows = load_observation(path)
    z = manifest["redshift"]
    hits = analyze(rows, z)["features"]
    candidates = [(i, h) for i, h in enumerate(hits, 1) if h["direction"] == "absorption"
                  and abs(h["rest_wavelength"] - K) < 2]
    if len(candidates) != 1:
        raise ValueError(f"Expected one K candidate for {path}, found {len(candidates)}")
    index, hit = candidates[0]
    requests = [{"feature_id": f"F{index:04d}", "check": check} for check in
                ("mask_neighborhood", "continuum_sensitivity")]
    mask, continuum = measure(requests, hits, rows)
    # Force a measurement at the predicted H position, independent of detector ranking.
    h_target = H * (1 + z)
    pixel = min(range(len(rows)), key=lambda i: abs(float(rows[i]["wavelength"]) - h_target))
    fixed = [{"center_pixel": pixel, "wavelength": float(rows[pixel]["wavelength"])}]
    fixed_h = measure([{"feature_id": "F0001", "check": "continuum_sensitivity"}], fixed, rows)[0]
    h_neighbors = [h for h in hits if h["direction"] == "absorption"
                   and H - 2 <= h["rest_wavelength"] <= H_EPSILON + 2]
    return {"observation": manifest["observation_id"], "z": z, "k_pred": K * (1 + z),
            "k_observed": hit["wavelength"], "k_rest": hit["rest_wavelength"],
            "k_delta_rest": hit["rest_wavelength"] - K, "k_id": mask["feature_id"],
            "k_score": hit["score"], "k_contrast": continuum["fits"],
            "k_masked_neighbors": mask["unusable_pixels"], "k_core_or": mask["core_or_mask"],
            "h_pred": h_target, "h_nearest_pixel": fixed[0]["wavelength"],
            "h_fixed_contrast": fixed_h["fits"],
            "h_region_blind_centroids_rest": [h["rest_wavelength"] for h in h_neighbors]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cohort", type=Path)
    args = parser.parse_args()
    records = [investigate(args.cohort / f"sdss-6497-56329-{fiber:04d}.json")
               for fiber in (88, 96, 100)]
    fixed_observer_wave = mean(r["k_observed"] for r in records)
    print(json.dumps({"vacuum_k": K, "vacuum_h": H, "vacuum_h_epsilon": H_EPSILON,
                      "rms_fixed_observer": mean((r["k_observed"] - fixed_observer_wave) ** 2
                                                 for r in records) ** .5,
                      "rms_ca_k": mean((r["k_observed"] - r["k_pred"]) ** 2
                                       for r in records) ** .5,
                      "observations": records}, indent=2, allow_nan=False))


if __name__ == "__main__":
    main()
