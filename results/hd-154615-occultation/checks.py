"""Summarize input-level alternatives and conditional timing uncertainty."""
import csv
import gzip
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize_scalar

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def timing_intervals():
    timings = read(ROOT / "timings.csv")
    rng = np.random.default_rng(73)
    for kind in ("primary", "secondary"):
        cycles, draws = [], []
        for sector in (12, 39, 66, 93):
            group = [r for r in timings if r["flux_column"] == "detrended"
                     and r["eclipse"] == kind and int(r["sector"]) == sector]
            oc = [float(r["oc_minutes"]) for r in group]
            cycles.append(np.mean([float(r["cycle"]) for r in group]))
            draws.append(np.median(rng.choice(oc, (20000, len(oc))), axis=1))
        for index, sector in ((1, 39), (2, 66)):
            weight = (cycles[index]-cycles[0]) / (cycles[-1]-cycles[0])
            delay = draws[index] - ((1-weight)*draws[0] + weight*draws[-1])
            print(kind, sector, "conditional 95% interval, min", np.percentile(delay, [2.5, 97.5]))


def light_curves():
    for num, intervals in [(32772, [(8.88, 9.24), (10.67, 10.88)]),
                           (1309, [(12.70, 16.06), (14.38, 15.34)])]:
        lines = (INPUT / f"occ_{num}.txt").read_text().splitlines()
        duration = float(lines[0].split(": ")[-2])
        flux = np.array(lines[-1].split(":")[1:], dtype=float)
        t = np.linspace(0, duration, len(flux))
        masks = [(t > lo + 1e-7) & (t < hi - 1e-7) for lo, hi in intervals]
        outside = ~(masks[0] | masks[1])
        base = np.median(flux[outside])
        print(num, "baseline and two event levels (instrumental units)",
              base, np.median(flux[masks[0] & ~masks[1]]), np.median(flux[masks[1]]))


def neighbors():
    for star in read(INPUT / "gaia_sources.csv"):
        path = INPUT / f"neighbors_{star['source_id']}.csv"
        if not path.exists():
            continue
        ra, dec = float(star["ra"]), float(star["dec"])
        measures = []
        for other in read(path):
            if star["source_id"] == other["source_id"] or not other["phot_g_mean_mag"]:
                continue
            east = (float(other["ra"]) - ra) * np.cos(np.deg2rad(dec))
            sep = 3600 * np.hypot(east, float(other["dec"]) - dec)
            delta_mag = float(other["phot_g_mean_mag"]) - float(star["phot_g_mean_mag"])
            measures.append((sep, delta_mag))
        print(star["source_id"], "closest neighbor (arcsec, delta G)", min(measures),
              "brightest neighbor (arcsec, delta G)", min(measures, key=lambda r: r[1]))


def independent_template():
    lightcurves = []
    for sector in (12, 39, 66, 93):
        file = INPUT / f"qlp_5965875226751271040_{sector}.csv.gz"
        with gzip.open(file, "rt") as stream:
            data = np.genfromtxt(stream, names=True, delimiter=",")
        good = np.isfinite(data["detrended"]) & (data["quality"] == 0)
        flux = data["detrended"][good]
        lightcurves.append((data["btjd"][good], flux / np.median(flux)))
    for name, phase in (("primary", 0), ("secondary", 0.4975)):
        def fold(t):
            return (((t-1644.95)/1.851076 - phase + 0.5) % 1 - 0.5) * 1.851076

        template_t, template_y = lightcurves[0]
        grid = np.linspace(-0.16, 0.16, 33)
        values = [np.median(template_y[abs(fold(template_t)-point) < 0.013]) for point in grid]
        template = CubicSpline(grid, gaussian_filter1d(values, 0.7))
        shifts, cycles = [], []
        for t, whole_flux in lightcurves:
            x = fold(t)
            window = abs(x) < 0.12
            x, flux = x[window], whole_flux[window]
            cycles.append(np.mean((t[window]-1644.95)/1.851076))

            def loss(shift):
                basis = np.column_stack([np.ones_like(x), x, template(x-shift)-1])
                residual = flux - basis @ np.linalg.lstsq(basis, flux, rcond=None)[0]
                return np.mean(np.sqrt(1 + (residual/0.004)**2))

            shifts.append(1440*minimize_scalar(loss, method="bounded", bounds=(-0.025, 0.025)).x)
        endpoint = np.interp(cycles, [cycles[0], cycles[-1]], [shifts[0], shifts[-1]])
        print(name, "template delay from endpoint clock (min)", np.array(shifts)-endpoint)


if __name__ == "__main__":
    timing_intervals()
    light_curves()
    neighbors()
    independent_template()
