"""Construct four independent TESS separation measurements and block-bootstrap samples."""
from pathlib import Path
import json
import numpy as np
from tess_template import GROUPS, fit
from orbit import EP, P

HERE = Path(__file__).resolve().parent
DATA = HERE / 'inputs'

def measure(draws=600):
    rng = np.random.default_rng(20261114)
    meta = json.loads((DATA/'tess_inputs.json').read_text())
    alternative = json.loads((DATA/'alternate_reduction.json').read_text())
    results = {}
    with np.load(DATA/'orbital_folds.npz') as data:
        for group, sectors in GROUPS.items():
            flux, count, times, blocks = [], [], [], []
            for sector in sectors:
                n = data[f's{sector}_counts']
                cadence = meta[str(sector)]['cadence_days']
                keep = ((n.sum(axis=1)*cadence > .65*P) &
                        ((n>0).sum(axis=1) >= min(55,int(.65*P/cadence))))
                flux.append(data[f's{sector}_sums'][keep,:,0]); count.append(n[keep])
                cycles = data[f's{sector}_cycles'][keep]
                times.extend((EP + (cycles+.3)*P).tolist())
                # Do not bridge gaps in orbital-cycle number within a TESS sector.
                runs = np.split(np.arange(len(cycles)),np.flatnonzero(np.diff(cycles)!=1)+1)
                blocks.append([run[j:j+3] for run in runs for j in range(0,len(run),3)])
            point = fit(flux,count)[2]
            samples = []
            for _ in range(draws):
                ix = [np.concatenate([b[i] for i in rng.integers(0,len(b),len(b))]) for b in blocks]
                samples.append(float(fit([a[i] for a,i in zip(flux,ix)],
                                         [n[i] for n,i in zip(count,ix)])[2]))
            results[group] = dict(sectors=sectors, phase=point, wide_phase=fit(flux,count,True)[2],
                                  alternate=alternative[group]['centers_and_spacing'][2],
                                  orbits=len(times), hjd_by_orbit=times, bootstrap=samples,
                                  sigma=float(np.std(samples,ddof=1)))
            print(group, len(times), round(point,5), round(results[group]['sigma'],5),flush=True)
    (HERE/'tess_measurements.json').write_text(json.dumps(results,indent=2)+'\n')
if __name__=='__main__':
    measure()
