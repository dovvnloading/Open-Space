"""Compare adjacent covered even/odd dimmings at the reported epoch."""
import csv
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parent
with (root/'events.csv').open() as handle:
    events = list(csv.DictReader(handle))
rng = np.random.default_rng(744)
with (root/'parity.csv').open('w',newline='') as handle:
    out = csv.writer(handle,lineterminator='\n')
    out.writerow(['mpc','complete_even_odd_pairs','even_deeper','median_even_minus_odd_mag',
                  'bootstrap_low','bootstrap_high'])
    for ident in dict.fromkeys(r['mpc'] for r in events):
        depth = {int(r['cycle']):float(r['delta_mag']) for r in events if r['mpc']==ident}
        pairs = np.array([d-depth[k+1] for k,d in depth.items() if k%2==0 and k+1 in depth])
        bootstrap = np.median(rng.choice(pairs,(20000,len(pairs)),replace=True),axis=1)
        out.writerow([ident,len(pairs),sum(pairs>0),np.median(pairs),*np.quantile(bootstrap,[.025,.975])])
