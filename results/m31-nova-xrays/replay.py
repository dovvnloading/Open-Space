"""Regenerate the positional, catalog-history and fixed-aperture tests offline."""
import csv
import gzip
import math
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import binom, chi2

HERE = Path(__file__).resolve().parent
POINTING = (10.684708, 41.26875)
FIELDS = ("w", "b", "u", "s", "m", "h")


def read(name):
    path = HERE / name
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle))


def sky(ra, dec):
    ra, dec = map(math.radians, (float(ra), float(dec)))
    return [math.cos(dec) * math.cos(ra), math.cos(dec) * math.sin(ra), math.sin(dec)]


def separation(chord):
    return 3600 * math.degrees(2 * math.asin(chord / 2))


def screen():
    novae = read("nova_positions.csv")
    epochs = read("epoch_detections.csv")
    names = {r["name"] for r in epochs if r["obsid"] == "23871"} & {
        r["name"] for r in epochs if r["obsid"] == "26163"}
    unique = {r["name"]: r for r in epochs if r["name"] in names}
    sources = list(unique.values())
    catalog = cKDTree([sky(s["ra"], s["dec"]) for s in sources])
    print(f"Optical eruptions {len(novae)}, X-ray sources in both visits {len(sources)}")
    history = read("source_history.csv.gz")
    for nova in novae:
        vec = sky(nova["ra_deg"], nova["dec_deg"])
        for j in catalog.query_ball_point(vec, math.radians(2.5 / 3600)):
            src = sources[j]
            dist = separation(np.linalg.norm(np.array(vec) - sky(src["ra"], src["dec"])))
            early = [r for r in history if r["name"] == src["name"]
                     and float(r["gti_mjd_obs"]) < float(nova["jd"]) - 2400000.5]
            vals = [max((float(r["flux_significance_" + b]) for b in FIELDS
                         if r["flux_significance_" + b]), default=0.) for r in early]
            print(f'{nova["name"]} {src["name"]} {dist:.3f}" prior visits={len(vals)} '
                  f'prior with significance >=3: {sum(v >= 3 for v in vals)}')
    # Translate the optical list as a unit in the tangent plane near the M31 nucleus.
    def project(ra, dec):
        return [(float(ra) - POINTING[0]) * 3600 * math.cos(math.radians(POINTING[1])),
                (float(dec) - POINTING[1]) * 3600]
    tree = cKDTree([project(s["ra"], s["dec"]) for s in sources])
    xy = np.array([project(r["ra_deg"], r["dec_deg"]) for r in novae])
    for threshold in (0.5, 1, 2.5):
        actual = int(np.count_nonzero(tree.query(xy)[0] <= threshold))
        controls = []
        for radius in (5, 10, 20, 40):
            for angle in np.radians(np.arange(0, 360, 10)):
                offset = radius * np.array([np.cos(angle), np.sin(angle)])
                controls.append(int(np.count_nonzero(tree.query(xy + offset)[0] <= threshold)))
        print(f'Radius {threshold:g}" nearest-source matches={actual}, shifted catalogs '
              f'at least as high={sum(n >= actual for n in controls)}/{len(controls)}, '
              f'shifted range={min(controls)}..{max(controls)}')


def counts():
    photons = read("photons.csv.gz")
    observations = read("observations.csv")
    targets = read("targets.csv")
    print("\nFixed optical-position counts, source aperture 2 arcsec, annulus 4-8 arcsec")
    for target in targets:
        name = target["name"]
        print(name)
        seasons = {}
        for obs in observations:
            oid = obs["obsid"]
            year = obs["date_obs"][:4]
            detector = "HRC" if "/hrcf" in obs["url"] else "ACIS"
            radial = [math.hypot(float(p["x_offset_arcsec"]), float(p["y_offset_arcsec"]))
                      for p in photons if p["target"] == name and p["obsid"] == oid]
            on = sum(r <= 2 for r in radial)
            off = sum(4 < r <= 8 for r in radial)
            exposure = float(obs["exposure_s"])
            values = seasons.setdefault((year, detector), np.zeros(3))
            values += [on, off, exposure]
            if year == "2021":
                epoch = datetime.fromisoformat(obs["date_obs"]) - datetime(1858, 11, 17)
                delay = epoch.total_seconds() / 86400 + 2400000.5 - float(target["jd"])
                print(f"  {oid} delay={delay:.2f} d on={on} off={off} net={on-off/12:.2f} "
                      f"local_binomial_tail={binom.sf(on-1,on+off,1/13):.3g}")
        for (year, detector), (on, off, exposure) in sorted(seasons.items()):
            net = on - off / 12
            print(f"  {year} {detector} {exposure/1000:.3f}ks on={on:.0f} off={off:.0f} "
                  f"net={net:.2f} rate={1000*net/exposure:.3f} counts/ks")


        # One-sided Garwood limits for four independent Poisson means. The union
        # bound gives at least 95% joint coverage for the resulting ratio bound.
        early_on, early_off, early_time = seasons["2021", "HRC"]
        late_on, late_off, late_time = seasons["2023", "HRC"] + seasons["2024", "HRC"]
        tail = 0.05 / 4
        def lower(n):
            return 0.5 * chi2.ppf(tail, 2*n) if n else 0.0
        def upper(n):
            return 0.5 * chi2.ppf(1-tail, 2*(n+1))
        first = (lower(early_on) - upper(early_off)/12) / early_time
        last = (upper(late_on) - lower(late_off)/12) / late_time
        print(f"  95% lower bound on decline in net HRC aperture rate: {first/last:.2f}x")


if __name__ == "__main__":
    screen()
    counts()
