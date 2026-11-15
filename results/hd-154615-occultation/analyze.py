"""Recompute the occultation/EB cross-match and time individual HD 154615 eclipses."""
import csv
import gzip
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
HD = "5965875226751271040"
EPOCH = 1644.95  # BTJD = BJD(TDB) - 2457000
PERIOD = 1.851076  # plotting reference; curvature is independent of a linear ephemeris
PER_DAY = 1440


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def csv_out(path, data):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(data[0]))
        writer.writeheader()
        writer.writerows(data)


def screen():
    events = rows(INPUT / "occultation_screen.csv")
    eb = {r["source_id"]: r for r in rows(INPUT / "gaia_eclipsing.csv")}
    stars = {r["source_id"]: r for r in rows(INPUT / "gaia_sources.csv")}
    result = []
    for event in events:
        source = event["GAIA_ID"]
        if source not in eb:
            continue
        parallax = float(stars[source]["parallax"])
        period = 1 / float(eb[source]["frequency"])
        for number in range(1, 5):
            separation = event[f"DOUB_SEP_SOLN_{number}"]
            if not separation:
                continue
            projected = float(separation) / parallax
            mass_floor = (projected / 2) ** 3 / (period / 365.25) ** 2
            result.append(dict(source=source, event=event["SEQ_NUM"], period_days=period,
                               solution=number, separation_mas=separation,
                               position_angle_deg=event[f"DOUB_PA_SOLN_{number}"],
                               projected_au=projected, two_body_floor_solar=mass_floor))
    assert len(events) == 136 and len({r["GAIA_ID"] for r in events}) == 135
    assert {r["source"] for r in result} == set(eb) and len(eb) == 3
    return result


def fit_eclipses(sector, data, column):
    t, flux = data["btjd"], data[column]
    good = np.isfinite(t + flux) & (data["quality"] == 0)
    t, flux = t[good], flux[good] / np.median(flux[good])
    output = []
    for name, phase in [("primary", 0), ("secondary", 0.4975)]:
        for cycle in range(int((min(t) - EPOCH) / PERIOD) - 2,
                           int((max(t) - EPOCH) / PERIOD) + 3):
            prediction = EPOCH + (cycle + phase) * PERIOD
            window = abs(t - prediction) < 0.18
            x, y = t[window] - prediction, flux[window]
            if len(y) < 12 or min(x) > -0.082 or max(x) < 0.082 or min(abs(x)) > 0.015:
                continue

            def residual(params):
                level, slope, depth, shift, width = params
                return level + slope*x - depth*np.exp(-0.5*((x-shift)/width)**2) - y

            fit = least_squares(residual, [1, 0, 0.04, 0, 0.041],
                                bounds=([0.75, -0.6, 0, -0.06, 0.012],
                                        [1.25, 0.6, 0.15, 0.06, 0.085]),
                                loss="soft_l1", f_scale=0.003)
            depth, shift = fit.x[2:4]
            if depth >= 0.01 and abs(shift) <= 0.04:
                output.append(dict(sector=sector, flux_column=column, eclipse=name, cycle=cycle,
                                   center_btjd=prediction+shift, oc_minutes=shift*PER_DAY,
                                   fractional_depth=depth, cadences=len(y)))
    return output


def summaries(timings):
    out = []
    for column in ("detrended", "small", "large", "sap"):
        for name in ("primary", "secondary"):
            epochs = []
            for sector in (12, 39, 66, 93):
                samples = [r for r in timings if (r["flux_column"], r["eclipse"], r["sector"])
                           == (column, name, sector)]
                mid = np.mean([r["cycle"] for r in samples])
                oc = np.median([r["oc_minutes"] for r in samples])
                epochs.append((mid, oc))
            cycles, oc = np.array(epochs).T
            original_clock = np.polyfit(cycles[[0, 3]], oc[[0, 3]], 1)
            excess = oc - np.polyval(original_clock, cycles)
            for sector, (cycle, m), delay in zip((12, 39, 66, 93), epochs, excess):
                out.append(dict(flux_column=column, eclipse=name, sector=sector,
                                mean_cycle=cycle, median_oc_minutes=m,
                                delay_from_endpoint_clock_min=delay))
    return out


def main():
    matched = screen()
    timings = []
    for path in sorted(INPUT.glob(f"qlp_{HD}_*.csv.gz")):
        sector = int(path.name.split("_")[-1].split(".")[0])
        with gzip.open(path, "rt") as stream:
            data = np.genfromtxt(stream, delimiter=",", names=True)
        for column in ("detrended", "small", "large", "sap"):
            timings.extend(fit_eclipses(sector, data, column))
    summary = summaries(timings)
    csv_out(ROOT / "separations.csv", matched)
    csv_out(ROOT / "timings.csv", timings)
    csv_out(ROOT / "timing_summary.csv", summary)
    for row in summary:
        if row["sector"] in (39, 66):
            print(row["flux_column"], row["eclipse"], row["sector"],
                  f'{row["delay_from_endpoint_clock_min"]:.2f} min')


if __name__ == "__main__":
    main()
