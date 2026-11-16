"""Refit the two leading period families with the archived CFHT exposure 2616044p."""
import csv
import sys

import numpy as np
from scipy.optimize import least_squares

import fit


def main():
    with (fit.BASE.parent / "archive" / "measurements.csv").open() as stream:
        frame = next(row for row in csv.DictReader(stream) if row["id"] == "2616044p")
    mjd = float(frame["epoch"])
    ephemeris = (fit.BASE / "cfht_horizons.txt").read_text()
    data = ephemeris.split("$$SOE")[1].split("$$EOE")[0]
    # Horizons is tabulated at exposure start. Convert to the same midpoint
    # and emission-time reference as the published orbital fit.
    epoch, = [row for row in csv.reader(data.splitlines())
              if row and abs(float(row[0]) - 2400000.5 - mjd) < 1e-6]
    distance = float(epoch[5])
    east_axis, north_axis = fit.tangent_axes([float(epoch[3])], [float(epoch[4])])
    hst = fit.rows[fit.reference]
    reference_time = float(hst["mjd_utc"]) - float(hst["distance_AU"]) * 499.004783836 / 86400
    time = mjd + float(frame["exptime"]) / 172800 - distance * 499.004783836 / 86400
    extra = np.array([float(frame["east"]), float(frame["north"])])
    sigma = 0.04

    def predict_extra(parameters):
        return fit.predict(parameters, np.array([time - reference_time]), east_axis,
                           north_axis, np.array([40 / distance]))[0]

    with (fit.BASE / "fit_all.csv").open() as stream:
        solutions = list(csv.DictReader(stream))[:2]
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(["frame", "sigma_arcsec", "nominal_period_days", "nominal_chi2",
                     "nominal_dof", "heldout_error_arcsec", "refit_period_days",
                     "refit_chi2", "refit_dof"])
    for solution in solutions:
        initial = np.array([float(solution[key]) for key in fit.PARAMETERS])
        sky = fit.signs[int(solution["signindex"]), :, None] * fit.sky
        signs = (-1, 1)

        def residual(parameters, sign):
            nominal = ((fit.predict(parameters) - sky) / fit.sigma[:, None]).ravel()
            return np.concatenate((nominal, (predict_extra(parameters) - sign * extra) / sigma))

        # Keep the other nights' published quadrant assignments and allow
        # either component assignment for the added image within each family.
        frequency = initial[0]
        bounds = ([frequency - .00035, .001, -20, .07, -30, -.999, -30],
                  [frequency + .00035, .96, 30, 10, 30, .999, 30])
        trials = [least_squares(residual, initial, args=(sign,), bounds=bounds) for sign in signs]
        if not all(trial.success for trial in trials):
            raise RuntimeError("CFHT sensitivity fit did not converge")
        best = min(trials, key=lambda trial: trial.cost)
        error = min(np.linalg.norm(predict_extra(initial) - sign * extra) for sign in signs)
        dof = sky.size - len(initial)
        writer.writerow([frame["id"], sigma, f"{1 / frequency:.3f}",
                         f"{float(solution['chi2']):.2f}", dof, f"{error:.5f}",
                         f"{1 / best.x[0]:.3f}", f"{2 * best.cost:.2f}", dof + 2])


if __name__ == "__main__":
    main()
