"""Historical stage amplification independent of modern SWOT data or rating curves."""
import csv
import json
from pathlib import Path
import numpy as np
from scipy.stats import theilslopes

ROOT = Path(__file__).resolve().parent
with (ROOT/'boye_kinsuka_1954.csv').open() as f:
    rows = [r for r in csv.DictReader(f) if r['boye_lower_stage_m']]


def summarize(rows):
    x = np.array([float(r['kinsuka_upper_stage_m']) for r in rows])
    y = np.array([float(r['boye_lower_stage_m']) for r in rows])
    # The datums of the two gauges cancel in changes and slopes.
    slope = theilslopes(y,x).slope
    return dict(n=len(rows),start=rows[0]['date'],end=rows[-1]['date'],upper_change_m=round(x[-1]-x[0],3),
                lower_change_m=round(y[-1]-y[0],3),fall_change_m=round(x[-1]-y[-1]-x[0]+y[0],3),
                lower_stage_response_per_m_upper=round(float(slope),3))


peak = max(range(len(rows)),key=lambda i: float(rows[i]['boye_lower_stage_m']))
result = dict(all_dates=summarize(rows),rising=summarize(rows[:peak+1]),recession=summarize(rows[peak:]))
(ROOT/'historical_summary.json').write_text(json.dumps(result,indent=2)+'\n')
