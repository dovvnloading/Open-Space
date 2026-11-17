"""Pre-review prototype for a conservative two-epoch Stokes comparison.

Consumes summary JSON only; this code cannot acquire or open event files.
No holdout is authorized by its existence. Bounds and gate attestations must
be independently justified and frozen before any holdout summary is made.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution


def vector_and_cov(summary, multiplier):
    z=np.zeros(3);c=np.zeros((3,3))
    for d in summary['du_sufficient']:
        a=multiplier*d['scale'];z+=np.array(d['source_sum'])-a*np.array(d['background_sum'])
        c+=np.array(d['source_outer'])+a*a*np.array(d['background_outer'])
    i,q,u=z
    if not np.isfinite(z).all() or i<=0:raise ValueError('non-positive or invalid net I')
    jac=np.array([[-q/i**2,1/i,0],[-u/i**2,0,1/i]])
    return z[1:]/i,jac@c@jac.T


def certified_epoch_bounds(summary, alpha_bounds, recovery_bounds):
    """Conservative component intervals and covariance trace upper bound.

    Numerator and denominator are affine in the shared multiplier. Each ratio
    is monotone while I remains positive; reciprocal recovery is positive.
    The covariance bound uses ||J||_F**2 * trace(C), avoiding optimization.
    """
    a0,a1=alpha_bounds; f0,f1=recovery_bounds
    src=np.zeros(3); bg=np.zeros(3); cs=np.zeros((3,3)); cb=np.zeros((3,3))
    for d in summary['du_sufficient']:
        k=float(d['scale']);src+=np.asarray(d['source_sum']);bg+=k*np.asarray(d['background_sum'])
        cs+=np.asarray(d['source_outer']);cb+=k*k*np.asarray(d['background_outer'])
    if np.linalg.eigvalsh(cs).min() < -1e-8 or np.linalg.eigvalsh(cb).min() < -1e-8:raise ValueError('invalid sufficient covariance')
    zs=[src-a*bg for a in (a0,a1)];imin=min(z[0] for z in zs)
    if not np.isfinite(imin) or imin<=0:raise ValueError('non-positive net I in nuisance range')
    values=np.array([z[1:]/(z[0]*f) for z in zs for f in (f0,f1)])
    maxq=max(abs(z[1]) for z in zs);maxu=max(abs(z[2]) for z in zs)
    trace_bound=(2/imin**2+(maxq**2+maxu**2)/imin**4)*np.trace(cs+a1*a1*cb)/f0**2
    return values.min(axis=0),values.max(axis=0),float(trace_bound)


def evaluate(baseline,holdout,config):
    required=['input_lineage_disjoint','gtis_valid','response_valid','regions_valid','centroid_rule_valid','spectra_within_family','nuisance_bounds_valid','statistical_approximation_valid']
    missing=[k for k in required if config['gates'].get(k) is not True]
    if missing:raise ValueError('unqualified gates: '+', '.join(missing))
    sys=np.array(config['difference_systematic_cov'],float)
    if sys.shape!=(2,2) or not np.allclose(sys,sys.T) or np.linalg.eigvalsh(sys).min() < -1e-12:raise ValueError('invalid systematic covariance')
    amin,amax=config['background_multiplier_bounds']; fmin,fmax=config['rmf_recovery_bounds']
    if not np.isfinite([amin,amax,fmin,fmax,config['threshold_chisq']]).all() or not (0<=amin<=amax and 0<fmin<=fmax and config['threshold_chisq']>0):raise ValueError('invalid nuisance bounds or threshold')
    results={}
    for key in ('nominal','xy_centroid'):
        def objective(x):
            ab,ah,fb,fh=x
            v1,c1=vector_and_cov(baseline[key],ab);v2,c2=vector_and_cov(holdout[key],ah)
            d=v2/fh-v1/fb;cov=c1/fb**2+c2/fh**2+sys
            if np.linalg.eigvalsh(cov).min()<=0:raise ValueError('singular covariance')
            return float(d@np.linalg.solve(cov,d))
        bounds=[(amin,amax),(amin,amax),(fmin,fmax),(fmin,fmax)]
        m=differential_evolution(objective,bounds,seed=438,polish=True,tol=1e-8)
        if not m.success:raise ValueError('nuisance minimization did not converge')
        results[key]={'minimum_chisq':float(m.fun),'minimizer':m.x.tolist()}
    threshold=config['threshold_chisq']
    for key in ('nominal','xy_centroid'):
        lo1,hi1,v1=certified_epoch_bounds(baseline[key],(amin,amax),(fmin,fmax))
        lo2,hi2,v2=certified_epoch_bounds(holdout[key],(amin,amax),(fmin,fmax))
        distance=np.maximum(np.maximum(lo1-hi2,lo2-hi1),0.)
        var_bound=v1+v2+np.trace(sys)
        if not np.isfinite(var_bound) or var_bound<=0:raise ValueError('invalid covariance bound')
        results[key]['certified_lower_chisq']=float(distance@distance/var_bound)
    return {'alternatives':results,'threshold_chisq':threshold,'all_alternatives_exceed':all(x['certified_lower_chisq']>threshold for x in results.values()),'warning':'certified bound is conservative conditional on declared model, gates, and bounds; numerical minima are diagnostic upper bounds only'}

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('baseline');ap.add_argument('holdout');ap.add_argument('config');args=ap.parse_args()
    def read(path):return json.loads(Path(path).read_text(encoding='utf-8'))
    print(json.dumps(evaluate(read(args.baseline),read(args.holdout),read(args.config)),indent=2))
