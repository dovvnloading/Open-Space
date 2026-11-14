"""Test whether the shape of [O III] 5008 predicts the separate 4960 line in DESI DR1."""
import json
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import ndtr

HERE = Path(__file__).resolve().parent
C = 299792.458
Z = 0.8416150156358213
LINES = {"strong": 5008.239, "weak": 4960.295}


@dataclass
class Window:
    wave: np.ndarray
    flux: np.ndarray
    ivar: np.ndarray
    response: np.ndarray
    rest: float

    @property
    def velocity(self):
        return (self.wave / (self.rest * (1 + Z)) - 1) * C

    def profile(self, velocity, sigma):
        center = self.rest * (1 + Z) * (1 + velocity / C)
        width = center * sigma / C
        dx = np.median(np.diff(self.wave))
        upper = (self.wave + dx / 2 - center) / width
        lower = (self.wave - dx / 2 - center) / width
        density = (ndtr(upper) - ndtr(lower)) / dx
        return self.response @ density


def window(data, aperture, name, halfwidth):
    key = f"{aperture}_{name}_"
    wave = data[key + "wave"]
    v = (wave / (LINES[name] * (1 + Z)) - 1) * C
    valid = (data[key + "mask"] == 0) & (data[key + "ivar"] > 0) & (abs(v) < halfwidth)
    return Window(wave, data[key + "flux"], np.where(valid, data[key + "ivar"], 0),
                  data[key + "R"], LINES[name])


def fit_amplitudes(w, velocity, sigma, broad=None, fraction=None, order=1):
    coordinate = (w.wave - w.rest * (1 + Z)) / 20
    continuum = np.polynomial.polynomial.polyvander(coordinate, order)
    profiles = [w.profile(velocity, sigma)]
    if broad is not None:
        wings = w.profile(velocity, broad)
        if fraction is None:
            profiles.append(wings)
        else:
            profiles = [(1 - fraction) * profiles[0] + fraction * wings]
    design = np.column_stack([continuum, *profiles])
    weights = np.sqrt(w.ivar)
    weighted = weights[:, None] * design
    params = np.linalg.lstsq(weighted, weights * w.flux, rcond=None)[0]
    model = design @ params
    residual = weights * (w.flux - model)
    return params, residual, model, np.linalg.pinv(weighted.T @ weighted)


def fit_shape(w, double, order):
    if double:
        starts = [[0, 12, 60], [0, 20, 100], [0, 8, 45]]
        bounds = ([-80, 3, 30], [80, 50, 200])
    else:
        starts = [[0, 20], [0, 50]]
        bounds = ([-80, 3], [80, 120])
    def residual(shape):
        return fit_amplitudes(w, *shape, order=order)[1]
    fits = [least_squares(residual, start, bounds=bounds) for start in starts]
    best = min(fits, key=lambda fit: fit.cost)
    amplitudes, _, model, _ = fit_amplitudes(w, *best.x, order=order)
    if np.any(amplitudes[order + 1:] <= 0):
        raise ValueError("Offset line-component fit is not physically admissible.")
    return best.x, amplitudes, 2 * best.cost, model


def comparison(data, halfwidth=680, order=1, retain_models=False):
    strong = window(data, "offset", "strong", halfwidth)
    weak = window(data, "offset", "weak", halfwidth)
    single, p1, chi1, m1 = fit_shape(strong, False, order)
    double, p2, chi2, m2 = fit_shape(strong, True, order)
    fraction = p2[-1] / (p2[-1] + p2[-2])
    pw1, r1, mw1, _ = fit_amplitudes(weak, *single, order=order)
    pw2, r2, mw2, _ = fit_amplitudes(weak, *double, fraction=fraction, order=order)
    result = {
        "halfwidth_kms": halfwidth, "continuum_order": order,
        "strong_pixels": int(np.count_nonzero(strong.ivar)), "weak_pixels": int(np.count_nonzero(weak.ivar)),
        "strong_chi2_one": chi1, "strong_chi2_two": chi2,
        "weak_chi2_transfer_one": float(r1 @ r1), "weak_chi2_transfer_two": float(r2 @ r2),
        "velocity_kms": double[0], "core_sigma_kms": double[1], "wing_sigma_kms": double[2],
        "wing_flux_fraction": fraction, "strong_to_weak_total_flux_ratio": (p2[-1] + p2[-2]) / pw2[-1],
        "strong_flux_one": p1[-1], "strong_flux_two": p2[-1] + p2[-2], "weak_flux_two": pw2[-1],
    }
    if retain_models:
        return result, single, double, fraction, [(strong, m1, m2), (weak, mw1, mw2)]
    return result


def sky_checks(data, shape, fraction):
    results = {}
    for name in LINES:
        w = window(data, "offset", name, 680)
        sky_flux = data[f"sky_{name}_flux"]
        sky_ivar = data[f"sky_{name}_ivar"]
        sky_mask = data[f"sky_{name}_mask"]
        scores = []
        for y, iv, mask in zip(sky_flux, sky_ivar, sky_mask):
            good = (mask == 0) & (iv > 0) & (w.ivar > 0)
            near = abs(w.wave - w.rest * (1 + Z)) < 3
            if good.sum() < 0.8 * (w.ivar > 0).sum() or (good & near).sum() < 4:
                continue
            blank = replace(w, flux=y, ivar=np.where(good, iv, 0))
            p, _, _, covariance = fit_amplitudes(blank, *shape)
            scores.append(p[-1] / np.sqrt(covariance[-1, -1]))
        p, _, _, covariance = fit_amplitudes(w, *shape)
        results[name] = {"eligible": len(scores), "median_wing_score": float(np.median(scores)),
                         "highest_positive_wing_score": float(max(scores)),
                         "target_wing_score_fixed_shape": p[-1] / np.sqrt(covariance[-1, -1])}
    return results


def central_comparison(data, shape, fraction):
    results = {}
    for name in LINES:
        w = window(data, "central", name, 680)
        def residual(v, widening):
            sigmas = np.hypot(widening, shape[1:])
            return fit_amplitudes(w, v, *sigmas, fraction=fraction)[1]
        direct = least_squares(lambda v: residual(v[0], 0), [0], bounds=([-100], [100]))
        broadened = least_squares(lambda p: residual(*p), [0, 50], bounds=([-100, 0], [100, 300]))
        err = np.sqrt(np.linalg.pinv(broadened.jac.T @ broadened.jac).diagonal())
        results[name] = {"chi2_scaled_copy": 2 * direct.cost, "chi2_extra_broadening": 2 * broadened.cost,
                         "extra_sigma_kms": broadened.x[1], "formal_sigma_error_kms": err[1]}
    return results


def plot(models):
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 2, figsize=(9, 5.5), height_ratios=[2.5, 1], sharex="col")
    for col, (w, m1, m2) in enumerate(models):
        keep = (abs(w.velocity) < 220) & (w.ivar > 0)
        v, y = w.velocity[keep], w.flux[keep]
        sig = 1 / np.sqrt(w.ivar[keep])
        a, b = axes[:, col]
        a.errorbar(v, y, sig, fmt="o", color="#222222", ms=3, lw=.9, label="DESI pixels")
        for model, color, label in [(m1, "#a84f26", "Single Gaussian"), (m2, "#00787a", "Core + wings")]:
            a.plot(v, model[keep], color=color, lw=1.7, label=label)
            b.plot(v, (y-model[keep])/sig, ".-", color=color, lw=1.1)
        a.set_title(["[O III] 5008 Å · fitted shape","[O III] 4960 Å · transferred shape"][col],loc="left")
        b.axhline(0,color="#444444",lw=.65)
        b.set_xlabel("Velocity from catalog redshift (km/s)")
    axes[0,0].legend(frameon=False,fontsize=8)
    axes[0,0].set_ylabel(r"Flux ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ Å$^{-1}$)")
    axes[1,0].set_ylabel("Residual / error")
    fig.suptitle("DESI 39628534716566336",fontsize=15,x=.10,ha="left")
    fig.tight_layout()
    fig.savefig(HERE / "oxygen-profile.png",dpi=170)
    plt.close(fig)


def main():
    with np.load(HERE / "windows.npz", allow_pickle=False) as data:
        primary, _, double, fraction, models = comparison(data, retain_models=True)
        result = {"primary":primary, "sky":sky_checks(data, double, fraction),
                  "foreground_position":central_comparison(data, double, fraction),
                  "sensitivities":[comparison(data, width, order)
                                   for width in (430,520,680,820) for order in (0,1,2)]}
        plot(models)
    (HERE / "results.json").write_text(json.dumps(result,indent=2) + "\n")


if __name__ == "__main__":
    main()
