"""Compare the UX Ari flare with the independent sector-44 catalog entry."""
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

from analyze import bytes_at, table_layout

# Xu et al. (2025), J/A+A/699/A322/tablea1. The three published columns
# are BeginTime, EndTime (BJD_TDB-2457000), and bolometric energy (erg).
CAL_START, CAL_END, PUBLISHED_ENERGY = 2501.147, 2501.183, 3.53e35
ORBIT_DAYS = 6.4372703  # Dümmler & Aarum 2001, Table 5

def compare(measurements):
    result = measurements["UX Ari"]
    data = bytes_at(result["tess_url"])
    h, offset, dtype = table_layout(data)
    d = np.frombuffer(data, dtype=dtype, count=int(h["NAXIS2"]), offset=offset)
    d = d[(d["QUALITY"] == 0) & np.isfinite(d["PDCSAP_FLUX"])]
    t = d["TIME"]
    inside = (t >= CAL_START) & (t <= CAL_END)
    edges = np.r_[CAL_START, t[inside], CAL_END]
    durations = []
    for key in ["SAP_FLUX", "PDCSAP_FLUX"]:
        flux = d[key].astype(float)
        for inner, outer in [(.5, 1.5), (1, 3), (1.5, 4), (2, 5)]:
            left = (t > CAL_START-outer/24) & (t < CAL_START-inner/24)
            right = (t > CAL_END+inner/24) & (t < CAL_END+outer/24)
            for degree in [1, 2]:
                model = np.polyfit(t[left|right]-2501.154, flux[left|right], degree)
                frac = flux/np.polyval(model, t-2501.154) - 1
                integral = np.trapezoid(np.interp(edges, t, frac), edges) * 86400
                durations.append(integral)
                if (key, inner, degree) == ("PDCSAP_FLUX", .5, 1):
                    adopted = integral
    ed = result["optical"]["equivalent_duration_s"]
    low, high = result["optical"]["equivalent_duration_sensitivity_s"]
    ref = datetime.fromisoformat(result["alert_utc"]).replace(tzinfo=timezone.utc).timestamp()/86400 + 2440587.5-2457000
    utc = t - d["TIMECORR"] - 69.184/86400
    phase_control = {}
    for cycle in [-2,-1,0,1,2]:
        hours = (utc - ref - cycle*ORBIT_DAYS) * 24
        bg = (abs(hours) >= 6) & (abs(hours) <= 12)
        sample = (hours >= -5) & (hours <= 3)
        if bg.sum()<100 or sample.sum()<80:
            continue
        flux = d["PDCSAP_FLUX"].astype(float)
        coef = np.polyfit(hours[bg], flux[bg], 1)
        phase_control[cycle] = float(np.max((flux/np.polyval(coef,hours) - 1)[sample]))
    return dict(reference_bjd_tdb_minus_2457000=[CAL_START,CAL_END], reference_published_erg=PUBLISHED_ENERGY,
                reference_ED_s=adopted, reference_ED_sensitivity_s=[min(durations),max(durations)],
                ED_ratio=ed/adopted, conditional_bolometric_erg=PUBLISHED_ENERGY*ed/adopted,
                conditional_baseline_bracket_erg=[PUBLISHED_ENERGY*low/max(durations),PUBLISHED_ENERGY*high/min(durations)],
                orbital_period_d=ORBIT_DAYS, phase_matched_maxima=phase_control)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    out = parser.parse_args().out
    measurements = json.loads((out/"measurements.json").read_text())
    (out/"uxari_comparison.json").write_text(json.dumps(compare(measurements), indent=2)+"\n")
