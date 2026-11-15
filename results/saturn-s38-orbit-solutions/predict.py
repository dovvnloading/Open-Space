"""Generate illustrative T14 sky positions for custom Julian dates in the frozen force grid."""

import argparse
import json

import numpy as np

from acquire import horizons, jd_list
from reproduce import DAY, EPOCH, GM_JUPITER, GM_SUN, LENGTH, ROOT, ShortArc, table, trajectory


def predict(epochs):
    grid = table("prediction_saturn")[:, 0]
    if min(epochs) < grid.min() + 1 or max(epochs) > grid.max() - 1:
        raise ValueError("Prediction epochs must be inside the frozen 2023–2026 force grid")
    tdb = np.array(epochs) + 69.182 / DAY
    _, observer = horizons(dict(COMMAND="'0'", CENTER="'T14'", EPHEM_TYPE="VECTORS",
                                REF_PLANE="FRAME", REF_SYSTEM="ICRF", VEC_TABLE="2",
                                CSV_FORMAT="YES", OUT_UNITS="KM-S", TLIST=jd_list(tdb)))
    model = ShortArc(table("astrometry")[table("astrometry")[:, 0] == 38])
    model.epochs = tdb - EPOCH
    model.observer = -np.array([[float(x) for x in row[2:5]] for row in observer])
    model.saturn = trajectory("prediction_saturn")
    model.perturbers = [(trajectory("prediction_sun"), GM_SUN), (trajectory("prediction_jupiter"), GM_JUPITER)]
    fitted = json.loads((ROOT / "results.json").read_text())
    states = [row["icrf_saturn_system_state_km_kms"] for row in fitted["solutions"]]
    states += [row["moderate_icrf_saturn_system_state_km_kms"] for row in fitted["night_exclusions"]]
    predictions = []
    for state in states:
        scaled = np.array(state)
        scaled[:3] /= LENGTH
        predictions.append(model.sky(scaled, rtol=2e-10))
    for i, epoch in enumerate(epochs):
        yield dict(jd_utc=epoch, moderate_ra_dec_deg=predictions[0][i].tolist(),
                   extreme_ra_dec_deg=predictions[1][i].tolist(),
                   moderate_night_deletion_ra_dec_deg=[sample[i].tolist() for sample in predictions[2:]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("jd_utc", nargs="+", type=float)
    print(json.dumps(list(predict(sorted(set(parser.parse_args().jd_utc)))), indent=2))
