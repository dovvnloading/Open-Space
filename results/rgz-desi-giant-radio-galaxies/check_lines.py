"""Measure forced [O II] line-window contrasts in the archived DESI coadds."""
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

root = Path(__file__).resolve().parent
for target in json.loads((root/'measurements.json').read_text()):
    tag = target['rgz'].split('_')[1][:7]
    obs = pd.read_csv(root/'inputs'/(tag+'_line.csv')).query('ivar > 0 and mask == 0')
    offset = np.abs(obs.wavelength-3727.42*(1+target['z']))
    core = obs[offset<4.5]
    background = obs[offset.between(10,28)]
    difference = np.average(core.flux, weights=core.ivar)-np.average(background.flux, weights=background.ivar)
    error = math.sqrt(1/core.ivar.sum()+1/background.ivar.sum())
    print(tag, 'forced contrast / diagonal-noise standard error', round(difference/error,2))
