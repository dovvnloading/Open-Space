"""Use the Swift UV data to score timing curves trained without these photons.

Evaluate the centers at each exposure. Refit brightness scales, yearly levels and
optionally constant UV phase offsets. The optical fit parameters remain frozen.
"""
from pathlib import Path
import csv
import json
import sys
import numpy as np
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parent.parent
sys.path.insert(0,str(ROOT/'joint'))
from joint import geometry,initial
from orbit import EP,P,TC,Y
fits=json.loads((ROOT/'joint/fits.json').read_text())['full']
table=list(csv.DictReader((ROOT/'external/swift-photometry.csv').open()))
templates=json.loads((ROOT/'external/inputs/ogle-model.json').read_text())['params']
def col(k):
    return np.array([float(row.get(k) or np.nan) for row in table])
def main(rad):
    t=col('hjd_mid');ye=np.array([int(r['start_utc'][:4]) for r in table])
    f,e,g,z=(col(f'{obj}_{rad}_{key}') for obj,key in [('target','rate'),('target','err'),('ref1','rate'),('ref1','err')])
    ok=np.isfinite(f+g+e+z+t)&(e>0)&(z>0)&(g>10)&(g<70)&(col(f'target_{rad}_exposure')>20)&(col(f'ref1_{rad}_exposure')>20)
    f,e,g,z,t,ye=[a[ok] for a in (f,e,g,z,t,ye)]
    flux=f/g;err=np.hypot(e/g,f*z/g**2)
    norm=np.median(flux);flux/=norm;err=np.hypot(err/norm,.06)
    year,indx=np.unique(ye,return_inverse=True);annual=np.eye(len(year))[indx]
    phase=(t-EP)/P%1;years=(t-TC)/Y;out={}
    for model in ('kepler','single_trough','cubic','constant'):
        p=np.array(fits[model]['parameters'])
        shift,_=geometry(model,p[2:2+len(initial(model)[0])],years)
        centers=p[0]+p[1]*years/10+shift+np.array([[0],[.5]])
        def evaluate(offsets):
            delta=(phase-centers-np.array(offsets)[:,None]+.5)%1-.5
            G=np.exp(-.5*(delta/np.array(templates[4:6])[:,None])**2)
            design=np.column_stack((annual,-G.T));weighted=design/err[:,None]
            coeff=np.linalg.lstsq(weighted,flux/err,rcond=None)[0]
            return (flux-design@coeff)/err
        base=evaluate([0,0])
        free=least_squares(evaluate,[0,0],bounds=(-.045,.045))
        out[model]=dict(n=len(t),fixed_centers_chi2=float(base@base),free_offsets_chi2=float(free.fun@free.fun),free_uv_center_offset_min=(free.x*P*1440).tolist())
    return out
if __name__=='__main__':
    r={s:main(s) for s in ['2.7','3.5','4.5']}
    (ROOT/'context/uv_prediction.json').write_text(json.dumps(r,indent=2)+'\n')
    print(r['3.5'])
