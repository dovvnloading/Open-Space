"""Leave one Mt. Lemmon station out and extrapolate across the inter-station gap."""
import gzip
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def main():
    archive = HERE.parent / "evidence/neocp_history_initial.json.gz"
    with gzip.open(archive, "rt") as f:
        case = next(r for r in json.load(f)["results"] if r["id"] == "CEXJ3W2")
    lines = case["observations"]
    minutes = []
    xy = []
    for line in lines:
        minutes.append((float(line[23:32]) - 29.4) * 1440)
        hour, minute, sec = map(float, line[32:44].split())
        deg, amin, asec = map(float, line[45:56].split())
        ra = 15 * (hour + minute/60 + sec/3600)
        dec = (1 if line[44] == "+" else -1) * (deg + amin/60 + asec/3600)
        xy.append([(ra - 29.4)*3600*np.cos(np.deg2rad(3.64)), (dec-3.64)*3600])
    minutes, xy = np.array(minutes), np.array(xy)
    stations = np.array([line[77:80] for line in lines])
    retained = np.array([line[64] != "x" for line in lines])
    first_i52 = np.flatnonzero(stations == "I52")[0]
    last_g96 = np.flatnonzero(stations == "G96")[-1]
    print("Inter-station gap (minutes)", minutes[first_i52] - minutes[last_g96])
    for label, sample in [("all", np.ones(len(lines), dtype=bool)), ("retained", retained)]:
        for degree in [1, 2]:
            for station, target in [("G96", first_i52), ("I52", last_g96)]:
                mask = (stations == station) & sample
                coeff = np.polynomial.polynomial.polyfit(minutes[mask], xy[mask], degree)
                predicted = np.polynomial.polynomial.polyval(minutes[target], coeff)
                miss = np.linalg.norm(predicted - xy[target])
                print(label, degree, station, "held-out miss (arcsec)", float(miss))


if __name__ == "__main__":
    main()
