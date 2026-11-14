"""Profile the center using only Swift's nine snapshots on 2012 September 25."""
import csv,json
from scipy.optimize import minimize_scalar
from common import D,np,PERIOD,PARAMS
rows=[r for r in csv.DictReader((D/'swift-photometry.csv').open()) if r['start_utc'].startswith('2012-09-25')]
def analyze(rad,width,floor):
    def read(s):return np.array([float(r[s]) for r in rows])
    x,e,r,re=(read(s) for s in [f'target_{rad}_rate',f'target_{rad}_err',f'ref1_{rad}_rate',f'ref1_{rad}_err'])
    y=x/r;y0=np.median(y);s=np.hypot(np.hypot(e/r,x*re/r**2),floor*y0)
    ph=read('phase');measured_exp=read(f'target_{rad}_exposure')/86400/PERIOD
    # Model the finite exposure, fit baseline and fractional deficit as linear nuisance parameters.
    offsets=np.linspace(-.5,.5,41)
    def solve(c):
        template=np.exp(-.5*((ph[:,None]-c+measured_exp[:,None]*offsets)/width)**2).mean(axis=1)
        X=np.c_[np.ones(len(y)),-template]
        a=np.linalg.lstsq(X/s[:,None],y/s,rcond=None)[0]
        val=np.sum(((y-X@a)/s)**2)
        if a[1]<0 or a[1]>a[0]:return np.inf,a
        return val,a
    grid=np.arange(.46,.55001,.00005);chi=np.array([solve(c)[0] for c in grid]);m=chi.argmin()
    ok=grid[chi<=chi[m]+3.84146];v,p=solve(grid[m])
    return dict(center=grid[m],conditional_95=[ok[0],ok[-1]],fractional_depth=p[1]/p[0],minimum_chi2=v)
out={f'{rad}_{w}_{fl}':analyze(rad,w,fl) for rad in [2.7,3.5,4.5] for w in [.018,PARAMS[5],.035] for fl in [.0,.03,.06]}
nom=out[f'3.5_{PARAMS[5]}_0.03'];d1=np.array([float(r['phase'])-float(r['delta_secondary_cycles']) for r in rows]);model=float(np.median(d1))
out['nominal']=nom;out['OGLE_secondary_phase']=model;out['primary_not_in_core']=True;out['n_snapshots']=len(rows)
(D/'swift-2012.json').write_text(json.dumps(out,indent=2));print(nom,model)
for w in [.018,PARAMS[5],.035]:print(w,out[f'3.5_{w}_0.03'])
