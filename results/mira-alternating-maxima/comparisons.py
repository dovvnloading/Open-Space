"""Replay six-decade comparison-star checks and V Del's field-matched ASAS observations."""

from collections import Counter, defaultdict
from statistics import mean, median
from replay import HERE, load_rows


def observer_checks(star, observations, maxima):
    times = {int(r["cycle"]): int(r["jd"]) for r in maxima if r["star"] == star}
    ref = {"V Del": ("MOW", "CMG"), "S Vir": ("MOW", "SHS"), "S Cyg": ("GMQA", "DRCA")}[
        star
    ]
    for band in ("Vis.", "V"):
        for width in (15, 25, 40):
            for validated in (False, True):
                samples = defaultdict(lambda: defaultdict(list))
                for obs in observations:
                    if (
                        obs["star"] != star
                        or obs["band"] != band
                        or obs["upper_limit"] != "0"
                    ):
                        continue
                    cycle = int(obs["cycle"])
                    if abs(float(obs["jd"]) - times[cycle]) > width:
                        continue
                    if obs["validation"] == "T" or (
                        validated and obs["validation"] != "V"
                    ):
                        continue
                    ob = obs["observer"]
                    for key in ("all", ob if ob in ref else "others"):
                        samples[key][cycle].append(float(obs["magnitude"]))
                for key in ("all", *ref, "others"):
                    peaks = {
                        e: median(v) for e, v in samples[key].items() if len(v) >= 2
                    }
                    contrasts = [
                        round(peaks[e + 1] - peaks[e], 5)
                        for e in peaks
                        if not e % 2 and e + 1 in peaks
                    ]
                    if contrasts:
                        print(
                            star,
                            band,
                            width,
                            "flag V" if validated else "exclude T",
                            key,
                            len(contrasts),
                            sum(v > 0 for v in contrasts),
                            sum(v == 0 for v in contrasts),
                            round(mean(contrasts), 3),
                        )


def phase_checks(star, observations, maxima):
    times = {int(row["cycle"]): int(row["jd"]) for row in maxima if row["star"] == star}
    points = defaultdict(list)
    for row in observations:
        if (
            row["star"] != star
            or row["band"] != "Vis."
            or row["validation"] != "V"
            or row["upper_limit"] != "0"
        ):
            continue
        cycle = int(row["cycle"])
        offset = float(row["jd"]) - times[cycle]
        if abs(offset) <= 25:
            points[cycle].append((row["observer"], offset, float(row["magnitude"])))
    for tolerance in (1, 2, 3, 5):
        contrasts = []
        for cycle in sorted(points):
            if cycle % 2 or cycle + 1 not in points:
                continue
            by_observer = defaultdict(list)
            for observer, phase, magnitude in points[cycle]:
                for observer2, phase2, magnitude2 in points[cycle + 1]:
                    if observer == observer2 and abs(phase2 - phase) <= tolerance:
                        by_observer[observer].append(magnitude2 - magnitude)
            if by_observer:
                contrasts.append(
                    round(median(median(v) for v in by_observer.values()), 6)
                )
        if contrasts:
            print(
                "matched observer/phase",
                star,
                tolerance,
                "pairs/positive/tied",
                len(contrasts),
                sum(v > 0 for v in contrasts),
                sum(v == 0 for v in contrasts),
                "mean",
                round(mean(contrasts), 3),
            )


def asas(maxima):
    series = []
    for row in (HERE / "v-del-asas3.txt").read_text().splitlines():
        if row.startswith("#dataset"):
            field = row.split(";")[-1].strip()
        elif row.strip() and not row.startswith("#"):
            cols = row.split()
            series.append(
                (float(cols[0]) + 2450000, float(cols[1]), cols[-2], cols[-1], field)
            )
    print("ASAS grades", dict(Counter(row[2] for row in series)))
    peaks = {int(r["cycle"]): int(r["jd"]) for r in maxima if r["star"] == "V Del"}
    field = "1 F2048+24_507"
    working = [r for r in series if r[2] in ("A", "B") and r[4] == field and r[1] < 29]
    for e in (23, 24, 25, 26, 27):
        near = sorted(
            (r for r in working if abs(r[0] - peaks[e]) <= 25),
            key=lambda r: abs(r[0] - peaks[e]),
        )
        print(
            "ASAS closest A/B point within 25 days in shared field",
            e,
            len(near),
            near[:1],
        )
    matches = [
        (abs((a[0] - peaks[26]) - (b[0] - peaks[27])), a, b)
        for a in working
        for b in working
        if abs(a[0] - peaks[26]) <= 30 and abs(b[0] - peaks[27]) <= 30
    ]
    print("ASAS nearest relative-phase match within 30 days, E26/E27", min(matches))


if __name__ == "__main__":
    maxima = load_rows("comparison-maxima.csv")
    obs = load_rows("comparison-near-maximum.csv")
    for star in ("V Del", "S Vir", "S Cyg"):
        observer_checks(star, obs, maxima)
        phase_checks(star, obs, maxima)
    asas(maxima)
