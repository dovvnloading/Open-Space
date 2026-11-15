"""Check the phase of HD 154615's eclipses in pre-TESS ASAS and Gaia photometry."""
import csv
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
INPUT = ROOT / "inputs"
PERIOD = 1.851076
T0 = 1644.95


def phase_shift(t, magnitude, fields, label):
    phase = (t - T0) / PERIOD
    levels = [fields == field for field in np.unique(fields)]

    def profile(shift_days):
        def eclipse(center, width):
            delta = ((phase - shift_days / PERIOD - center + 0.5) % 1) - 0.5
            return np.exp(-0.5*(delta/width)**2)
        template = eclipse(0, 0.023) + 0.55*eclipse(0.4975, 0.022)
        design = np.column_stack([*levels, template])
        coefficients = np.linalg.lstsq(design, magnitude, rcond=None)[0]
        residual = design @ coefficients-magnitude
        return np.sum(residual**2), coefficients[-1]

    # Search the complete orbital cycle. A positive amplitude means a fading.
    offsets = np.linspace(-PERIOD/2, PERIOD/2, 5289)
    results = np.array([profile(offset) for offset in offsets])
    scores = np.where(results[:, 1] > 0, results[:, 0], np.inf)
    best = np.argmin(scores)
    residual_variance = scores[best] / (len(t)-len(levels)-2)
    within = offsets[(scores-scores[best]) <= 3.84*residual_variance]
    return dict(series=label, count=len(t), median_time=np.median(t), earliest=t.min(), latest=t.max(),
                shift_min=offsets[best]*1440, lower_min=within.min()*1440, upper_min=within.max()*1440,
                confidence_segments=1+int(np.sum(np.diff(within)>1.5*(offsets[1]-offsets[0]))),
                depth_mag=results[best, 1])


def main():
    rows, field, expected_count = [], None, 0
    for line in (INPUT / "asas_170827-4207.5.txt").read_text().splitlines():
        if line.startswith("#ndata="):
            expected_count += int(line.split("=")[1])
        elif line.startswith("#dataset="):
            field = line.split(";")[1].split()[1]
        elif line.strip() and not line.startswith("#"):
            parts = line.split()
            rows.append((float(parts[0])-7000, parts[-2], field, [float(v) for v in parts[1:6]]))
    assert expected_count == len(rows) == 628
    results = []
    for ix, magcol in enumerate((3, 0, 1, 2, 4)):
        sel = [(t, mags[ix], field) for t, grade, field, mags in rows
               if grade in ("A", "B") and 9 < mags[ix] < 10]
        t, magnitude, field = zip(*sel)
        t, magnitude, field = np.array(t), np.array(magnitude), np.array(field)
        label = f"ASAS_MAG_{magcol}"
        results.append(phase_shift(t, magnitude, field, label))
        if magcol == 3:
            for half, mask in (("early", t < np.median(t)), ("late", t >= np.median(t))):
                results.append(phase_shift(t[mask], magnitude[mask], field[mask], label+"_"+half))
    with (ROOT.parent / "inputs" / "gaia_epoch_5965875226751271040.csv").open() as file:
        observations = [r for r in csv.DictReader(file) if r["variability_flag_g_reject"] == "False"
                        and r["rejected_by_photometry"] == "False"]
    t = np.array([float(r["g_transit_time"])-1802.5 for r in observations])
    magnitude = np.array([float(r["g_transit_mag"]) for r in observations])
    results.append(phase_shift(t, magnitude, np.zeros_like(t), "Gaia_G"))
    with (ROOT / "early_phase.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, list(results[0]))
        writer.writeheader()
        writer.writerows(results)
    for result in results:
        print(result["series"], result["count"], round(result["shift_min"], 2),
              round(result["lower_min"], 2), round(result["upper_min"], 2))


if __name__ == "__main__":
    main()
