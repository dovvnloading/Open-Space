"""Report the raw aperture, background and centroid control for TOI-5697."""
import sys
from pathlib import Path

import numpy as np

with np.load(Path(sys.argv[1]) / '416739490-44.npz') as product:
    delta = abs(product['TIME'] - 2514.5623)
    valid = product['QUALITY'] == 0
    event = (delta < .06) & valid
    shoulders = (delta > .23) & (delta < .4) & valid
    print('field,event_median,shoulder_median,difference')
    for name in ('SAP_FLUX', 'PDCSAP_FLUX', 'SAP_BKG', 'MOM_CENTR1'):
        core = np.nanmedian(product[name][event])
        outside = np.nanmedian(product[name][shoulders])
        print(f'{name},{core:.8f},{outside:.8f},{core - outside:.8f}')
