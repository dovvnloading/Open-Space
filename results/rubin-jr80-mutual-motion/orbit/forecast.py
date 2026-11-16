"""Compare the two best orbital families on observable nights around 2027 opposition."""
import csv
from datetime import date, timedelta
from pathlib import Path

import numpy as np
from fit import BASE, PARAMETERS, predict, rows, tangent_axes


def compute(ephemeris):
    """Horizons observer columns: JD UTC, flags, ICRF coordinates, az/el, range, elongation."""
    section = ephemeris.split("$$SOE")[1].split("$$EOE")[0]
    fields = [r for r in csv.reader(section.splitlines()) if len(r) >= 10]
    data = np.asarray([[float(r[k]) for k in (0, 3, 4, 6, 7, 9)] for r in fields])
    jd, ra, dec, elevation, distance, elongation = data.T
    east, north = tangent_axes(ra, dec)
    ref = next(r for r in rows if r["id"] == "HST")
    t0 = float(ref["mjd_utc"])-float(ref["distance_AU"])*499.004783836/86400
    time = jd-2400000.5 - distance*499.004783836/86400 - t0
    solutions = list(csv.DictReader((BASE/"fit_all.csv").open()))[:2]
    predictions = [predict([float(s[k]) for k in PARAMETERS], time, east, north, 40/distance)
                   for s in solutions]
    difference = np.minimum(np.linalg.norm(predictions[0]-predictions[1], axis=1),
                            np.linalg.norm(predictions[0]+predictions[1], axis=1))
    visible = np.where((elevation > 60) & (elongation > 160))[0]
    best = visible[np.argmax(difference[visible])]
    second = [i for i in visible if abs(jd[i]-jd[best]) >= 10]
    choices = sorted([best, second[np.argmax(difference[second])]])
    forecasts = []
    for i in choices:
        utc = date(1858, 11, 17) + timedelta(days=int(jd[i]-2400000.5))
        entry = dict(utc_0500=str(utc), ra_deg=ra[i], dec_deg=dec[i], elevation_deg=elevation[i],
                     range_AU=distance[i], sun_elongation_deg=elongation[i],
                     axial_model_difference_arcsec=round(difference[i], 5))
        for solution, vectors in zip(solutions, predictions):
            days = round(float(solution["period_days"]))
            e, n = vectors[i]
            entry[f"sep_{days}d_arcsec"] = round(np.hypot(e,n), 5)
            entry[f"pa_{days}d_mod180_deg"] = round(np.degrees(np.arctan2(e,n)) % 180, 4)
        forecasts.append(entry)
    return forecasts


if __name__ == "__main__":
    import sys
    result = compute((BASE/"opposition_horizons.txt").read_text())
    writer = csv.DictWriter(sys.stdout, fieldnames=result[0], lineterminator="\n")
    writer.writeheader()
    writer.writerows(result)
