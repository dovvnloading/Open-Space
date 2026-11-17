"""Summarize predeclared spectral windows in microJy, grouping temporally close frames."""
import csv
import sys
from pathlib import Path

import numpy as np

BANDS = [(3.20, 3.45), (4.05, 4.19), (4.20, 4.35), (4.40, 4.62), (4.62, 4.82)]


def weighted(f, sigma):
    weights = sigma**-2
    return np.average(f, weights=weights), np.sum(weights)**-0.5


def main(inputs, output):
    rows = []
    for fn in inputs:
        with Path(fn).open() as f:
            rows.extend(csv.DictReader(f))
    result = []
    for name in sorted({r['target'] for r in rows}):
        accepted = sorted((r for r in rows if r['target']==name and r['status']=='fit'),
                          key=lambda r:float(r['mjd']))
        if not accepted:
            continue
        mjd = np.array([float(r['mjd']) for r in accepted])
        visits = np.cumsum(np.r_[0, np.diff(mjd)>20])
        wave = np.array([float(r['wavelength_um']) for r in accepted])
        flux = np.array([float(r['flux_ujy']) for r in accepted])
        sig = np.array([float(r['sigma_ujy']) for r in accepted])
        for a,b in BANDS:
            use = (a<=wave)&(wave<b)
            if not np.any(use):
                continue
            mean, formal = weighted(flux[use],sig[use])
            vi = np.unique(visits[use]); leave=[]
            for v in vi:
                cut = use & (visits!=v)
                if cut.any():
                    leave.append(weighted(flux[cut],sig[cut])[0])
            jack = np.sqrt((len(leave)-1)*np.var(leave)) if len(leave)>1 else None
            result.append(dict(target=name,lo_um=a,hi_um=b,n_frame=int(use.sum()),
                        n_visits=len(vi),flux_ujy=round(mean,3),formal_ujy=round(formal,3),
                        visit_jackknife_ujy=round(jack,3) if jack is not None else ''))
    if result:
        with Path(output).open('w') as out:
            w=csv.DictWriter(out,fieldnames=list(result[0]));w.writeheader();w.writerows(result)


if __name__=='__main__':
    main(sys.argv[1:-1],sys.argv[-1])
