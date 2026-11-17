"""Run a hot-star template difference test on two disjoint coadd pairs from the giant SB1."""
import argparse
import json
from pathlib import Path
import numpy as np
from hot_spectra import HERE, MODELS, PAIRS, ROOT, SID
from inject import OFFSETS, run

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',required=True,type=Path)
    a=p.parse_args()
    results=[]
    for pair in PAIRS:
        for temp in ('hot10000','hot12000'):
            r=run(a.cache,pair,np.unique(np.r_[OFFSETS,np.arange(-150,151,3)]),source=SID,
                  catalog=ROOT/'extension/lamost/exposures.csv',gamma=-51.526,mass_ratio=1.25,
                  model_files={'wave':MODELS['wave'],temp:MODELS[temp]},rotations=[15.,150.],
                  splits=[(0,0)],fractions=[0,.05,.10,.15,.20,.30],write_scores=False)
            r['template']=temp
            results.append(r)
    (HERE/'hot_ccf.json').write_text(json.dumps(results,indent=2)+'\n')
