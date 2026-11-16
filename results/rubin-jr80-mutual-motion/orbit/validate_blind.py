"""Compare frozen forecasts or the final refit with the accepted DECam axes."""
import argparse
import csv
from pathlib import Path
import sys

import numpy as np

BASE = Path(__file__).resolve().parent


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--final-fit", action="store_true",
                        help="evaluate the final 103- and 155-day fits on the same images")
    args = parser.parse_args()
    measured = read(BASE.parent / "archive" / "measurements.csv")
    forecasts = []
    if args.final_fit:
        from historical_predictions import SETS, calculate

        solutions = read(BASE / "fit_all.csv")[:2]
        ref = next(row for row in read(BASE / "observations.csv") if row["id"] == "HST")
        for _, ephemeris, exposure in SETS:
            for row in calculate(ephemeris, exposure, solutions, ref):
                forecasts.append(("final_refit", row))
    else:
        for label in ("blind", "deep", "second", "last"):
            for row in read(BASE / f"{label}_predictions.csv"):
                forecasts.append((label, row))
    writer = csv.DictWriter(sys.stdout, fieldnames=["frame","prediction_batch","mjd","period_days",
                            "measured_sep","measured_pa","predicted_sep","predicted_pa",
                            "axial_pa_error_deg","axial_vector_error_arcsec"], lineterminator="\n")
    writer.writeheader()
    for frame in measured:
        if frame["qc"] != "blind axis passes":
            continue
        projected = np.array([float(frame["east"]),float(frame["north"])])
        for name, model in forecasts:
            if abs(float(model["mjd"])-float(frame["epoch"])) > .00001:
                continue
            vector = np.array([float(model["east"]),float(model["north"])])
            residual = min(np.linalg.norm(projected-vector), np.linalg.norm(projected+vector))
            angle = abs((float(model["pa_mod180"])-float(frame["pa"])+90)%180-90)
            writer.writerow(dict(frame=frame["id"], prediction_batch=name, mjd=frame["epoch"],
                                 period_days=model["P"], measured_sep=frame["sep"], measured_pa=frame["pa"],
                                 predicted_sep=model["separation"], predicted_pa=model["pa_mod180"],
                                 axial_pa_error_deg=round(angle,5), axial_vector_error_arcsec=round(residual,6)))
