"""Empirical Gaussian feature through shifted fiducial responses; no molecular identity.
Input is output of swan_response.py. Model outputs remain aperture MJy/sr sums.
"""
import argparse,json
from pathlib import Path
from swan_fit_sample import residual_dof, verify_fit_sample
import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import differential_evolution
ap=argparse.ArgumentParser();ap.add_argument('--input',required=True);ap.add_argument('--output',required=True);ap.add_argument('--calibration-root',required=True);args=ap.parse_args()
rows=json.load(open(args.input));cal=Path(args.calibration_root);epochs=['2025-12-22','2025-12-25','2026-01-02'];rows=[r for r in rows if r['start'][:10] in epochs and r['map_bc_um']>4.1 and r['diagnostic_outlier_inclusive'].get('2.5')]
verify_fit_sample(rows)
wave=np.linspace(4.05,4.48,4301);resp=[]
for r in rows:
 ch=r['channel'];p=cal/f'MiddlePixel/SubChannel17/Band5/FiducialSpec_Band5_SubCh{ch:02d}_17_Cent.csv';s=np.genfromtxt(p,delimiter=',',names=True);w=s['Wavelength_um']+r['shift_um'];R=s['NormTransmission'];resp.append(np.interp(wave,w,R,left=0,right=0)/trapezoid(R/w,w))
resp=np.array(resp);bc=np.array([r['map_bc_um'] for r in rows]);group=np.array([epochs.index(r['start'][:10]) for r in rows]);train=group<2;test=group==2

def modelshape(mu,sigma):
 g=np.exp(-.5*((wave-mu)/sigma)**2)/(sigma*np.sqrt(2*np.pi));return trapezoid(resp*g[None,:]*wave[None,:],wave,axis=1)
def run(aperture='2.5',floor=0.,linear=False):
 y=np.array([r['diagnostic_outlier_inclusive'][aperture]['sum_mjysr'] for r in rows]);e=np.array([r['diagnostic_outlier_inclusive'][aperture]['error_mjysr'] for r in rows]);e=np.sqrt(e**2+(floor*y)**2)
 def regress(t,mask,feature=True):
  X=np.column_stack([np.ones(mask.sum())]+([(bc[mask]-4.26)] if linear else [])+([t[mask]] if feature else [])); beta=np.linalg.lstsq(X/e[mask,None],y[mask]/e[mask],rcond=None)[0];pred=X@beta
  return dict(coefficients=beta.tolist(),chi2=float(np.sum(((y[mask]-pred)/e[mask])**2)),prediction=pred.tolist(),residual=(y[mask]-pred).tolist(),sigma_residual=((y[mask]-pred)/e[mask]).tolist(),n=int(mask.sum()),rank=int(np.linalg.matrix_rank(X)))
 def objective(p):
  t=modelshape(*p);fits=[regress(t,group==i) for i in [0,1]]
  # Require nonnegative feature amplitudes but permit signed local background.
  if any(f['coefficients'][-1]<0 for f in fits):return 1e12
  return sum(f['chi2'] for f in fits)
 opt=differential_evolution(objective,[(4.23,4.29),(.0005,.04)],seed=17,tol=1e-8,popsize=12,maxiter=100,polish=True)
 t=modelshape(*opt.x);tr=[regress(t,group==i) for i in [0,1]];te=regress(t,test);null=[regress(t,group==i,False) for i in range(3)]
 return dict(aperture=aperture,fractional_error_floor=floor,continuum='linear' if linear else 'constant',center_um=float(opt.x[0]),sigma_um=float(opt.x[1]),train_chi2=float(opt.fun),train_dof=residual_dof(tr, shared_parameters=2),train=tr,heldout=te,heldout_dof=residual_dof([te]),continuum_only=null,template=t.tolist(),observed=y.tolist(),errors=e.tolist())
results=[run(),run(floor=.1),run(aperture='3.5'),run(linear=True)]
base=results[0];best=base['train_chi2'];grid=[]
# Feasible shape range for a broad delta chi2 <=6 region; NOT confidence interval (zero df).
y=np.array(base['observed']);e=np.array(base['errors'])
for mu in np.linspace(4.23,4.29,61):
 for sig in np.linspace(.0005,.04,41):
  t=modelshape(mu,sig);chi=0.;valid=True
  for i in [0,1]:
   m=group==i;X=np.column_stack([np.ones(m.sum()),t[m]]);b=np.linalg.lstsq(X/e[m,None],y[m]/e[m],rcond=None)[0];chi+=np.sum(((y[m]-X@b)/e[m])**2);valid=valid and b[1]>=0
  if valid and chi<=best+6:grid.append([float(mu),float(sig),float(chi)])
report=dict(warning='Empirical response-shape test only; no molecular identification. Training has zero degrees of freedom with constant baselines; linear baseline underidentified. Heldout permits amplitude/baseline refit, not wholly parameter-free prediction. Noise floors are sensitivity assumptions, not calibrated errors.',observations=[dict(observationid=r['observationid'],start=r['start'],bc=r['map_bc_um']) for r in rows],fits=results,training_near_fit_grid=grid)
Path(args.output).write_text(json.dumps(report,indent=2));print(args.output)
for r in results:print(r['aperture'],r['continuum'],r['fractional_error_floor'],'center/sigma',r['center_um'],r['sigma_um'],'trainchi',r['train_chi2'],'heldoutchi',r['heldout']['chi2'],'nullheldoutchi',r['continuum_only'][2]['chi2'])
