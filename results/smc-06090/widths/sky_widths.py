"""Compute eclipse widths from the same TESS photons using sky/median detrending."""
import json
from pathlib import Path
from measure_widths import GROUPS, measurement
import numpy as np

HERE=Path(__file__).resolve().parent
with np.load(HERE/'sky_folds.npz') as d:
    result={}
    for group,sectors in GROUPS.items():
        result[group]=measurement([d[f's{s}_sums'] for s in sectors],
                                  [d[f's{s}_counts'] for s in sectors]).tolist()
(HERE/'sky_widths.json').write_text(json.dumps(result,indent=2)+'\n')
print(result)
