"""Recompute S Vir's first and later parity and inspect the missing 2001 catalogue peak."""

from collections import defaultdict
from statistics import mean, median
from replay import load_rows

maxima = {
    int(r["cycle"]): r
    for r in load_rows("comparison-maxima.csv")
    if r["star"] == "S Vir"
}
older = [r for r in load_rows("comparison-near-maximum.csv") if r["star"] == "S Vir"]
later = load_rows("s-vir-late-maximum.csv")
gap = load_rows("s-vir-gap.csv")


def summarize(points, center, band="Vis.", validated=False):
    by = defaultdict(list)
    for row in points:
        if row["band"] != band or row["upper_limit"] != "0" or row["validation"] == "T":
            continue
        if validated and row["validation"] != "V":
            continue
        cycle = int(row["cycle"])
        if abs(float(row["jd"]) - float(maxima[cycle]["jd"])) <= center:
            by[cycle].append(float(row["magnitude"]))
    return {e: median(vals) for e, vals in by.items() if len(vals) >= 2}


for source, points in [("1976–2000", older), ("2002–2018", later)]:
    for width in (15, 25, 40):
        for band in ("Vis.", "V", "TG"):
            for validated in (False, True):
                peaks = summarize(points, width, band, validated)
                contrasts = [
                    round(peaks[e] - peaks[e + 1], 5)
                    for e in sorted(peaks)
                    if e % 2 == 0 and e + 1 in peaks
                ]
                if contrasts:
                    print(
                        source,
                        width,
                        band,
                        "flag V" if validated else "exclude T",
                        "pairs/even fainter/tied",
                        len(contrasts),
                        sum(v > 0 for v in contrasts),
                        sum(v == 0 for v in contrasts),
                        "even minus odd",
                        round(mean(contrasts), 3),
                    )
midpoint = (int(maxima[18]["jd"]) + int(maxima[20]["jd"])) / 2
bins = defaultdict(list)
for row in gap:
    if row["band"] == "Vis." and row["validation"] == "V" and row["upper_limit"] == "0":
        day = float(row["jd"]) - midpoint
        bins[int(day // 20) * 20].append(float(row["magnitude"]))
print(
    "E18→20 span and midpoint", int(maxima[20]["jd"]) - int(maxima[18]["jd"]), midpoint
)
for lower, values in sorted(bins.items()):
    print("gap bin [day,day+20), N/median", lower, len(values), median(values))

# Repeat the same-observer and relative-phase test on the later era.
from comparisons import phase_checks

phase_checks("S Vir", later, load_rows("comparison-maxima.csv"))
