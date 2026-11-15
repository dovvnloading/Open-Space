"""Reproduce the HD 73712 orbit consistency checks with Python 3 (standard library)."""
import csv
import hashlib
import io
import sys
import itertools
import math
from pathlib import Path
from statistics import correlation, median

ROOT = Path(__file__).resolve().parent
INNER_DAYS = 48.717
PARALLAX = 5.3522
# P, epoch of periastron, eccentricity, relative semimajor axis in mas,
# argument of periastron of relative orbit (degrees), inclination (degrees).
ORBITS = {
    "2017": (35.45, 1987.86, .169, 140.2, 226.3, 90.4),
    "2024": (20.96, 2011.48, .980, 76.2, 149.2, 89.0),
    "2026": (20.83, 2011.34, .930, 70.7, 162.2, 89.5),
}


def eccentric_anomaly(mean, eccentricity):
    mean = (mean + math.pi) % math.tau - math.pi
    ecc = mean + math.copysign(.85 * eccentricity, mean)
    for _ in range(16):
        ecc -= (ecc - eccentricity * math.sin(ecc) - mean) / (1 - eccentricity * math.cos(ecc))
    return ecc


def separation(orbit, year):
    period, epoch, e, a, omega, inc = orbit
    ecc = eccentric_anomaly(math.tau * (year - epoch) / period, e)
    x = math.cos(ecc) - e
    y = math.sqrt(1 - e * e) * math.sin(ecc)
    cs, sn = math.cos(math.radians(omega)), math.sin(math.radians(omega))
    return a * math.hypot(x * cs - y * sn, (x * sn + y * cs) * math.cos(math.radians(inc)))


def hierarchy(period, e, mass_ratio, mutual_i):
    # Kepler's law gives a_out / a_in independently of distance.
    semimajor_ratio = ((period * 365.25 / INNER_DAYS) ** 2 * (1 + mass_ratio)) ** (1 / 3)
    actual = semimajor_ratio * (1 - e)
    limiting = 2.8 * ((1 + mass_ratio) * (1 + e) / math.sqrt(1 - e)) ** .4
    limiting *= 1 - .3 * mutual_i / 180
    return actual, actual / limiting


def max_e(period, mass_ratio, mutual_i):
    low, high = 0., .9999
    for _ in range(50):
        center = (low + high) / 2
        if hierarchy(period, center, mass_ratio, mutual_i)[1] > 1:
            low = center
        else:
            high = center
    return (low + high) / 2


def replay_screen(filename):
    raw = Path(filename).read_bytes()
    expected = "1f08bd509156c36f33f3be1e29243d73e0cadbbe9d23edd935232338f73e6362"
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError("Input differs from the screened ORB6 snapshot")
    rows = list(csv.reader(io.StringIO(raw.decode("latin1")), delimiter="|"))[2:]
    ranks = []
    for r in rows:
        try:
            if r[12] not in ("y", "d") or r[15] not in ("a", "m") or int(r[31]) > 5:
                continue
            years = float(r[11]) / (365.25 if r[12] == "d" else 1)
            if years < .2:
                continue
            e, w, i = float(r[25]), math.radians(float(r[27])), math.radians(float(r[17]))
            a = float(r[14]) * (1000 if r[15] == "a" else 1)
            proxy = .5 * 10 ** (-(float(r[7]) - 5) / 5)
            separation_proxy = a * abs(math.cos(i)) * (1 - e * e) / (1 + e * abs(math.sin(w)))
            ranks.append((separation_proxy / proxy, r[2], r[3]))
        except ValueError:
            continue
    print("Raw orbits, eligible, first ten", len(rows), len(ranks), sorted(ranks)[:10])
    assert len(rows) == 4051 and len(ranks) == 3339


def main():
    with (ROOT / "radial_velocities.csv").open() as f:
        rv = list(csv.DictReader(f))
    components = {c: {float(r["hjd"]): float(r["rv_kms"]) for r in rv if r["component"] == c}
                  for c in ["p", "s", "b"]}
    p, s = components["p"], components["s"]
    dates = sorted(p.keys() & s.keys())
    pearson = correlation([p[t] for t in dates], [s[t] for t in dates])
    print("RV measurements and simultaneous p/s pairs", len(rv), len(dates))
    print("Sharp-line Pearson correlation", pearson)
    print("Absolute splitting range km/s", min(abs(p[t]-s[t]) for t in dates), max(abs(p[t]-s[t]) for t in dates))
    residuals = []
    for t in dates:
        e = .11
        ecc = eccentric_anomaly(math.tau * (t - 2448302.6) / INNER_DAYS, e)
        true = 2 * math.atan2(math.sqrt(1 + e) * math.sin(ecc / 2),
                              math.sqrt(1 - e) * math.cos(ecc / 2))
        w = math.radians(258)
        model = (31.6 + 36.7) * (math.cos(true + w) + e * math.cos(w))
        residuals.append((t, p[t] - s[t] - model))
    print("Median absolute published SB2 residual, max outlier", median(abs(v) for t, v in residuals),
          max(residuals, key=lambda tv: abs(tv[1])))
    with (ROOT / "later_astrometry.csv").open() as f:
        astrometry = list(csv.DictReader(f))
    for name, orbit in ORBITS.items():
        period, _, e, a, _, _ = orbit
        print("Orbit", name, "periastron au", a * (1 - e) / PARALLAX)
        for r in astrometry:
            yr = float(r["epoch"]); measured = float(r["sep_mas"])
            pred = separation(orbit, yr)
            print(yr, "measured, predicted, residual mas", measured, round(pred, 3), round(measured - pred, 3))
        qs = [.3, 2.]
        tests = [hierarchy(period, e, q, i) for q, i in itertools.product(qs, [0, 180])]
        for j, label in enumerate(["periastron / inner semimajor axis", "periastron / stability limit"]):
            print(label, min(v[j] for v in tests), max(v[j] for v in tests))
        limits = [max_e(period, q, i) for q, i in itertools.product(qs, [0, 180])]
        print("Threshold eccentricities", min(limits), max(limits))
    assert len(rv) == 61 and len(dates) == 16 and pearson < -.94
    assert 60 < separation(ORBITS["2017"], 2024.1538) < 62
    for name in ["2024", "2026"]:
        per, _, e, *_ = ORBITS[name]
        assert hierarchy(per, e, .3, 180)[1] < .5
    print("Verified")


if __name__ == "__main__":
    main()
    if len(sys.argv) > 1:
        replay_screen(sys.argv[1])
