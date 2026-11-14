"""Conditional parameter and prediction uncertainty; cadence-preserving rival-model calibration."""
import json
from pathlib import Path
import numpy as np
from joint import JointFit, load, stats, timed_prediction, HOURS, HERE
from orbit import EP, P, Y, TC, SCALE, circ
from concurrent.futures import ProcessPoolExecutor
import argparse

O,T=load()
F=json.loads((HERE/'fits.json').read_text())
TG=list(T)
DATES=[2461587.5,2461953.5,2462683.5]
HELD=O['h']>=2459580.5
OGROUP=O['survey']
def opt_pred(f,p):
    f.evaluate(p)
    return f.optical_pred.copy()
def model_profile(generator):
    f=JointFit(O,T,generator)
    p=np.asarray(F['full'][generator]['parameters'])
    f.evaluate(p)
    return f.optical_pred.copy(),np.array(f.tess_pred)
def gen(seed,generator):
    rng=np.random.default_rng(seed)
    op,pt=model_profile(generator)
    # A thirty-day random sign per survey and time block keeps correlated optical residuals together.
    blocks=np.c_[OGROUP,np.floor((O['h']-O['h'].min())/30).astype(int)]
    _,labels=np.unique(blocks,axis=0,return_inverse=True)
    syn=dict(O)
    syn['mag']=op+(O['mag']-op)*rng.choice([-1.,1.],labels.max()+1)[labels]
    ts={}
    for k,pred in zip(T,pt):
        draw=np.array(T[k]['bootstrap'])
        ts[k]=dict(T[k]);ts[k]['phase']=float(pred+rng.choice(draw)-np.mean(draw))
    return syn,ts
def rough_loss(r):
    return np.sum(2*SCALE**2*(np.sqrt(1+(r/SCALE)**2)-1))
def simulate(args):
    i,generator=args
    syn,ts=gen(357247+i,generator)
    full={};held={}
    for kind in ['kepler','cubic','sinusoid_offset','single_trough']:
        for cut in ['full','before_2022']:
            f=JointFit(syn,ts,kind,None if cut=='full' else ~HELD,None if cut=='full' else TG[:2])
            sol=f.solve(np.array(F[cut][kind]['parameters']))
            if cut=='full':
                full[kind]={'loss':float(2*sol.cost),'period_years':sol.x[2] if kind in ('kepler','sinusoid_offset') else None,
                            'ecc':sol.x[3] if kind=='kepler' else None,
                            'forecast':(HOURS*timed_prediction(kind,sol.x,DATES)).tolist()}
            else:
                scored=JointFit(syn,ts,kind)
                scored.evaluate(sol.x)
                rr=(scored.optical_pred[HELD]-syn['mag'][HELD])/syn['error'][HELD]
                tr=np.array([(scored.tess_pred[j]-ts[k]['phase'])/ts[k]['sigma'] for j,k in enumerate(TG) if j>=2])
                held[kind]=dict(loss=rough_loss(rr)+sum(tr**2),tess_loss=sum(tr**2),
                                pred_h=HOURS*scored.tess_pred[-1],
                                pred_2023_h=HOURS*scored.tess_pred[-2])
    return full,held
def calc_observed():
    r={}
    for k in ['kepler','cubic','sinusoid_offset','single_trough']:
        row=F['before_2022'][k]
        r[k]=row['held_late_ogle_loss']+sum(row['tess'][s]['standardized_residual']**2 for s in TG[2:])
    return r
def main():
    argv=argparse.ArgumentParser();argv.add_argument('--draws',type=int,default=300);argv.add_argument('--workers',type=int,default=2);args=argv.parse_args()
    quant=lambda x: np.quantile(x,[.025,.5,.975]).tolist()
    result={};observed=calc_observed()
    for generator in ['kepler','cubic','sinusoid_offset','single_trough']:
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            draws=list(pool.map(simulate,[(i,generator) for i in range(args.draws)],chunksize=4))
        diffs=np.array([held['cubic']['loss']-held['kepler']['loss'] for full,held in draws])
        free_diffs=np.array([held['sinusoid_offset']['loss']-held['kepler']['loss'] for full,held in draws])
        full_diffs=np.array([full['cubic']['loss']-full['kepler']['loss'] for full,held in draws])
        plateau_diffs=np.array([held['single_trough']['loss']-held['kepler']['loss'] for full,held in draws])
        out=dict(observed_held_losses=observed,
                 held_trough_minus_kepler_quantiles=quant(plateau_diffs),
                 n_as_large_as_observed_held_trough_advantage=int(np.sum(plateau_diffs>=observed['single_trough']-observed['kepler'])),
                 full_trough_minus_kepler_quantiles=quant([f['single_trough']['loss']-f['kepler']['loss'] for f,h in draws]),
                 n_as_large_as_observed_full_trough_advantage=int(sum(f['single_trough']['loss']-f['kepler']['loss']>=F['full']['single_trough']['loss']-F['full']['kepler']['loss'] for f,h in draws)),
                 held_cubic_minus_kepler_quantiles=quant(diffs),
                 n_as_large_as_observed_held_cubic_advantage=int(np.sum(diffs>=observed['cubic']-observed['kepler'])),
                 held_free_minus_kepler_quantiles=quant(free_diffs),
                 n_as_large_as_observed_held_free_advantage=int(np.sum(free_diffs>=observed['sinusoid_offset']-observed['kepler'])),
                 full_cubic_minus_kepler_quantiles=quant(full_diffs),
                 n_as_large_as_observed_full_advantage=int(np.sum(full_diffs>=F['full']['cubic']['loss']-F['full']['kepler']['loss'])),
                 conditional_forecast_quantiles_h={date:quant([f[generator]['forecast'][j] for f,h in draws])
                        for j,date in enumerate(['2027-07-01','2028-07-01','2030-07-01'])},
                 params_period_years=quant([f[generator]['period_years'] for f,h in draws]) if generator in ('kepler','sinusoid_offset') else None,
                 conditional_before_2022_fit_forecast_2025_h=quant([h[generator]['pred_h'] for f,h in draws]),
                 conditional_before_2022_fit_forecast_2023_h=quant([h[generator]['pred_2023_h'] for f,h in draws]),
                 params_ecc=quant([f[generator]['ecc'] for f,h in draws]) if generator=='kepler' else None,
                 draws=[{'full':a,'held':b} for a,b in draws])
        result[generator]=out
        print(generator,args.draws, [out[k] for k in ['params_period_years','params_ecc','held_cubic_minus_kepler_quantiles','n_as_large_as_observed_held_cubic_advantage','n_as_large_as_observed_full_advantage']],flush=True)
        (HERE/'uncertainty.json').write_text(json.dumps(result,indent=1)+'\n')
if __name__=='__main__':main()
