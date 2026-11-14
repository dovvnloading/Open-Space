"""Check numerical integration, dilution, direction, and the independent OGLE V photometry."""
import argparse
import json
from pathlib import Path
from hashlib import sha256
import numpy as np
from scipy.optimize import least_squares, minimize_scalar, brentq
from geometry import HERE, KEYS, P, EPOCH, TREF, LightCurve, fit, load, scores

VHASH = {
    'ogle3-V.dat':'6972583288cdb0fdc701041066f0fd5ffca00d88aa5fad528c46ed0861cad6fa',
    'ogle4-V.dat':'4df5182d3cbc75922c598b6f59f223fcd95b0dfc0dd9fb77970a41322c8af2eb'
}

def save(name,value):
    (HERE/(name+'.json')).write_text(json.dumps(value,indent=2)+'\n')

def summary(curve,q):
    out=[]
    for year in (2002.,2011.,2018.5,2023.5,2025.6,2027.5,2029.5):
        # Fractional years here use a Julian year since 2000 Jan 1, 12:00.
        t=(2451545+(year-2000)*365.25-TREF)/365.25
        tmp=q.copy()
        tmp[11:13]=0.
        def mag(ph): return curve.predict_at(tmp,np.array([t]),np.array([ph]))[0]
        stats=[]
        for seed in (0.,.5):
            opt=minimize_scalar(lambda v:-mag(v),bounds=(seed-.114,seed+.13),method='bounded')
            center,depth=opt.x,-opt.fun
            level=depth/2
            left=brentq(lambda x:mag(x)-level,center-.19,center)
            right=brentq(lambda x:mag(x)-level,center,center+.19)
            stats.append([center,depth,(right-left)*P*1440])
        out.append(dict(year=year,primary=stats[0],secondary=stats[1],
                        separation_hours=(stats[1][0]-stats[0][0])*P*24))
    return out

def external(curve,q,u,third):
    chunks=[]
    for g,label in enumerate(VHASH):
        p=HERE/'inputs'/label
        a=np.loadtxt(p)
        assert sha256(p.read_bytes()).hexdigest()==VHASH[label]
        if g==0:a[:,0]+=2450000.
        a=a[np.isfinite(a).all(axis=1)&(a[:,2]>.0)&(a[:,2]<=.1)]
        chunks.append(np.c_[a,np.full(len(a),g)])
    h,m,err,grp=np.concatenate(chunks).T
    c=LightCurve((h,m,np.hypot(err,.012),grp.astype(int)),u,third,curve.sense,order=32)
    sel=np.isin(np.arange(len(q)),[8,9,10])
    start=np.array([1.,18.1,18.1])
    def resid(v):
        fresh=q.copy()
        fresh[sel]=v
        return c.residuals(fresh)
    sol=least_squares(resid,start,bounds=([.5,17.5,17.5],[1.7,18.8,18.8]),
                     loss='soft_l1',f_scale=1.8,ftol=1e-8)
    return dict(loss=2*sol.cost,n=len(m),calibrated=['intensity_ratio','base_iii','base_iv'],
                calib=sol.x.tolist())

def calc():
    data=load(HERE/'inputs')
    base=json.loads((HERE/'full.json').read_text())
    out={}
    for sense in (1,-1):
        fits=[]
        for third in (0.,.1,.2,.25,.30,.35,.40):
            c=LightCurve(data,.4,third,sense,order=32)
            init=base[str(sense)]['pars'] if not fits else fits[-1]['pars']
            rs=fit(c,start=init)
            q=np.array(rs['pars'])
            coarse=LightCurve(data,.4,third,sense,order=16)
            high=LightCurve(data,.4,third,sense,order=64)
            rs.update(third=third,curve_stats=summary(c,q),V_conditional=external(c,q,.5,third),
                      max_integration_step_mmag=1000*float(np.max(abs(c.predict(q)-high.predict(q)))),
                      high_order_loss=float(scores(high.residuals(q)).sum()),
                      coarse_order_loss=float(scores(coarse.residuals(q)).sum()))
            fits.append(rs)
            print(sense,third,rs['loss'],rs['V_conditional']['loss'],rs['inclination_deg'],flush=True)
        out[str(sense)]=fits
    save('dilution',out)

def paired_boot(draws):
    # Nonparametric night-cluster resampling within each of the six disjoint time blocks.
    data=load(HERE/'inputs')
    full=json.loads((HERE/'full.json').read_text())
    curves=[LightCurve(data,.4,0.,s) for s in (1,-1)]
    rng=np.random.default_rng(771493)
    h=data[0]
    labels=np.digitize(h,[2453500,2455000,2457000,2458300,2459700])
    nights=np.floor(h-.5).astype(int)
    groups=[[np.where((labels==g)&(nights==n))[0] for n in np.unique(nights[labels==g])] for g in range(6)]
    res=[]
    for k in range(draws):
        chosen=np.concatenate([a[idx] for a in groups for idx in rng.integers(0,len(a),len(a))])
        v=[fit(c,chosen,full[str(c.sense)]['pars']) for c in curves]
        res.append([v[1]['loss']-v[0]['loss'],v[0]['pars'][2],v[0]['pars'][3],v[0]['inclination_deg']])
        if k%20==0:print(k,flush=True)
    save('bootstrap',dict(keys=['loss_reverse_minus_forward','U_years','e','inclination_deg'],draws=res,
                         caveat='Conditional on dilution=0, common u=0.4, same circularly symmetric stars, fixed six-block target selection.', 
                         quantiles=np.percentile(res,[2.5,50,97.5],axis=0).tolist(),
                         fraction_forward_lower=float(np.mean(np.array(res)[:,0]>0))))
if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--bootstrap',type=int,default=0)
    args=parser.parse_args()
    if args.bootstrap:paired_boot(args.bootstrap)
    else:calc()
