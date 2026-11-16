"""Generic documented PSG CO2 shape sensitivity; not a SWAN production-rate retrieval."""
import argparse,json,hashlib
from pathlib import Path
from swan_fit_sample import verify_fit_sample
import numpy as np
from scipy.integrate import trapezoid
from scipy.optimize import lsq_linear
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);p.add_argument('--calibration-root',required=True);p.add_argument('--model-root',required=True);a=p.parse_args()
epochs=['2025-12-22','2025-12-25','2026-01-02']; rows=json.load(open(a.input)); rows=[r for r in rows if r['start'][:10] in epochs and 4.1<r['map_bc_um']<4.45 and r['diagnostic_outlier_inclusive'].get('2.5')]
verify_fit_sample(rows)
group=np.array([epochs.index(r['start'][:10]) for r in rows]);bc=np.array([r['map_bc_um'] for r in rows]);cal=Path(a.calibration_root);mr=Path(a.model_root);out=[]; provenance={}
for t in [20,60,100,200]:
 f=mr/f'co2_shape_{t}K_response.txt'
 if not f.exists() or '# Molecules considered: CO2' not in f.read_text():continue
 s=np.loadtxt(f);wave=s[:,0];mol=s[:,-1];mol=mol/trapezoid(mol,wave);shape=[]
 provenance[str(t)]={'spectrum_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'config_sha256':hashlib.sha256((mr/f'co2_shape_{t}K.cfg').read_bytes()).hexdigest()}
 for r in rows:
  ch=r['channel'];z=np.genfromtxt(cal/f'MiddlePixel/SubChannel17/Band5/FiducialSpec_Band5_SubCh{ch:02d}_17_Cent.csv',delimiter=',',names=True);w=z['Wavelength_um']+r['shift_um'];R=z['NormTransmission'];shape.append(trapezoid(np.interp(wave,w,R,left=0,right=0)*mol*wave,wave)/trapezoid(R/w,w))
 shape=np.array(shape)
 for ap in ['2.5','3.5']:
  for floor in [0.,.1]:
   for linear in [False,True]:
    y=np.array([r['diagnostic_outlier_inclusive'][ap]['sum_mjysr'] for r in rows]);err=np.array([r['diagnostic_outlier_inclusive'][ap]['error_mjysr'] for r in rows]);err=np.hypot(err,floor*y);fits=[];null=[]
    for i in range(3):
     m=group==i
     for feature,save in [(True,fits),(False,null)]:
      X=np.column_stack([np.ones(m.sum())]+([bc[m]-4.26] if linear else [])+([shape[m]] if feature else []));low=np.full(X.shape[1],-np.inf)
      if feature:low[-1]=0
      beta=lsq_linear(X/err[m,None],y[m]/err[m],bounds=(low,np.full(X.shape[1],np.inf))).x;pred=X@beta
      save.append(dict(epoch=epochs[i],coefficients=beta.tolist(),chi2=float(np.sum(((y[m]-pred)/err[m])**2)),dof=int(m.sum()-X.shape[1]),prediction=pred.tolist(),observed=y[m].tolist(),sigma_residual=((y[m]-pred)/err[m]).tolist()))
    out.append(dict(temperature_K=t,aperture=ap,fractional_floor=floor,continuum='linear' if linear else 'constant',train_chi2=sum(f['chi2'] for f in fits[:2]),fits=fits,continuum_only=null,template=shape.tolist()))
report=dict(warning='Generic pure CO2 GSFC fluorescence spectrum, arbitrary geometry and normalization; only shape tested. No SWAN temperature or production-rate retrieval. Photon weighting assumes ideal flat Fnu calibration and needs independent radiometric verification. Training selects among four temperatures; heldout refits nonnegative amplitude and signed continuum. Ten percent floor is assumed sensitivity only.',api='https://psg.gsfc.nasa.gov/api.php',official_template='https://psg.gsfc.nasa.gov/data/templates/Comet_67P_from_Rosetta.txt',input_sha256=hashlib.sha256(Path(a.input).read_bytes()).hexdigest(),models=provenance,observations=[dict(id=r['observationid'],start=r['start'],bc=r['map_bc_um']) for r in rows],fits=out)
Path(a.output).write_text(json.dumps(report,indent=2))
for r in out:
 if r['aperture']=='2.5':print(r['temperature_K'],r['fractional_floor'],r['continuum'],r['train_chi2'],r['fits'][2]['chi2'])
