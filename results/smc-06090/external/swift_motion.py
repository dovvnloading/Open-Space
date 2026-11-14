"""Compare predictions trained on OGLE and estimate eclipse motion using Swift alone."""
from common import D,np,PARAMS,PERIOD,shape
from scipy.optimize import least_squares
from astropy.table import Table
from pathlib import Path
import csv,json
def fit_centers(phase,ye,y,e,groups,extra_widths=False):
    years,iy=np.unique(ye,return_inverse=True);gr,ig=np.unique(groups,return_inverse=True);n=len(years);ng=len(gr)
    def fn(p):
        c1,c2,d1,d2=p[n:n+4*ng].reshape(4,ng)
        s1,s2=p[-2:] if extra_widths else PARAMS[4:6]
        k1=(phase-c1[ig]+.5)%1-.5;k2=(phase-c2[ig]+.5)%1-.5
        return (y-p[iy]+d1[ig]*np.exp(-.5*(k1/s1)**2)+d2[ig]*np.exp(-.5*(k2/s2)**2))/e
    start=np.r_[np.ones(n),np.zeros(ng),np.full(ng,.52),np.full(2*ng,.34)]
    lo=np.r_[np.full(n,.55),np.full(ng,-.045),np.full(ng,.46),np.zeros(2*ng)]
    hi=np.r_[np.full(n,1.5),np.full(ng,.07),np.full(ng,.59),np.full(2*ng,.85)]
    if extra_widths:start=np.r_[start,PARAMS[4:6]];lo=np.r_[lo,.014,.014];hi=np.r_[hi,.04,.04]
    f=least_squares(fn,start,bounds=(lo,hi),max_nfev=400)
    return f,n,ng
def main():
    data=list(csv.DictReader((D/'swift-photometry.csv').open()))
    def ft(k):return np.array([float(r.get(k) or np.nan) for r in data])
    out={};rng=np.random.default_rng(20261114)
    for rad in [2.7,3.5,4.5]:
        h=ft('hjd_mid');ye=np.array([int(r['start_utc'][:4]) for r in data]);e1=ft(f'target_{rad}_err');e2=ft(f'ref1_{rad}_err');r=ft(f'ref1_{rad}_rate');x=ft(f'target_{rad}_rate')
        ok=np.isfinite(h+x+e1+e2+r)&(ft(f'target_{rad}_exposure')>20)&(ft(f'ref1_{rad}_exposure')>20)&(r>10)&(r<70)&(ye>=2016)&(e1>0)&(e2>0)
        h,ye,x,e1,e2,r=[s[ok] for s in [h,ye,x,e1,e2,r]]
        phase,delta1,delta2=shape(h);year=(h-2456000)/365.25
        scale=np.median((x/r)[(abs(delta1)>.075)&(abs(delta2)>.075)])
        y=x/r/scale;e=np.hypot(np.hypot(e1/r,x*e2/r**2)/scale,.06)
        years,ys=np.unique(ye,return_inverse=True);eye=np.eye(len(years))[ys];g=(ye>=2021).astype(int)
        nights=np.floor(h-.5)
        def errors(f):
            A=np.linalg.pinv(f.jac.T@f.jac);u=np.unique(nights);scores=np.array([f.jac[nights==i].T@f.fun[nights==i] for i in u]);return A@scores.T@scores@A*(len(y)-1)/(len(y)-len(f.x))*len(u)/(len(u)-1)
        result={}
        for name,free in [('fixed_width',False),('free_width',True)]:
            f,n,ng=fit_centers(phase,ye,y,e,g,free);cov=errors(f);a,b,d1,d2=f.x[n:n+4*ng].reshape((4,ng));v=(b-a)*PERIOD*24;contrast=np.zeros(len(f.x));contrast[n:n+4]=[1,-1,-1,1]
            errs=[]
            for j in range(2):
                row=np.zeros(len(f.x));row[n+j]=-1;row[n+2+j]=1;errs.append(float(np.sqrt(row@cov@row)*PERIOD*24))
            result[name]=dict(groups={'early':sorted(set(ye[g==0].tolist())),'late':sorted(set(ye[g==1].tolist()))},counts=[int((g==i).sum()) for i in [0,1]],intervals_h=v.tolist(),se_h=errs,delta_min=float(np.diff(v)[0]*60),delta_se_min=float(np.sqrt(contrast@cov@contrast)*PERIOD*1440),phase_primary=a.tolist(),phase_secondary=b.tolist(),depth_primary=d1.tolist(),depth_secondary=d2.tolist(),widths_cycles=f.x[-2:].tolist() if free else PARAMS[4:6],robust_loss=float((f.fun**2).sum()))
        if rad==3.5:
            f,n,ng=fit_centers(phase,ye,y,e,np.zeros(len(g),int),False)
            result['constant_fit_chi2']=sum(f.fun**2)
            result['constant_gap_h']=float((f.x[n+1]-f.x[n])*PERIOD*24)
            result['constant_comparison']='Two centers and depths optimized on Swift; note free two-era model additionally frees eclipse depths.'
            # Interpolate between a freely located static pair and the OGLE-predicted differential motion.
            # A shared linear shift accommodates error in the catalog orbital period.
            yr=(h-2458500.)/365.25;yr-=yr.mean()
            raw_c1=phase-delta1;raw_c2=phase-delta2
            raw_c1=(raw_c1+.5)%1-.5;raw_c2=(raw_c2+.5)%1-.5
            raw_c2=np.where(raw_c2<0,raw_c2+1,raw_c2)
            change=(raw_c2-raw_c1);change-=change.mean()
            def with_scale(params):
                c1,c2,slope,ampl,d1,d2=params[-6:]
                shift=ampl*change/2
                sep1=(phase-c1-slope*yr+shift+.5)%1-.5
                sep2=(phase-c2-slope*yr-shift+.5)%1-.5
                return (y-(eye@params[:-6]-d1*np.exp(-.5*(sep1/PARAMS[4])**2)-d2*np.exp(-.5*(sep2/PARAMS[5])**2)))/e
            initial=np.r_[np.ones(len(years)),0.,.53,0.,1.,.4,.32]
            lower=np.r_[np.full(len(years),.55),-.1,.45,-.012,-3,0,0]
            upper=np.r_[np.full(len(years),1.5),.1,.62,.012,4,.85,.85]
            moving=least_squares(with_scale,initial,bounds=(lower,upper))
            sm=errors(moving)
            chi=[]
            for fac in [0.,1.]:
                def fixed(q):
                    return with_scale(np.insert(q,len(q)-2,fac))
                keep=np.arange(len(initial))!=len(initial)-3
                fm=least_squares(fixed,moving.x[keep],bounds=(lower[keep],upper[keep]))
                chi.append(sum(fm.fun**2))
            result['motion_amplitude']={'scale':float(moving.x[-3]),'night_cluster_SE':float(np.sqrt(sm[-3,-3])),'chi2_at_scale_0_and_1':chi,'free_chi2':sum(moving.fun**2),
             'interpretation':'Free primary and secondary zero-points and a common linear drift. Scale 0 fixes separation; scale 1 imposes OGLE predicted variation. Both widths set to OGLE flux template.'}
            vals=[]
            for it in range(999):
                # Matched observation-night bootstrap, stratified by year. Each interval contains observations from several years.
                bidx=[]
                for kk in years:
                    ng=np.unique(nights[ye==kk]);draw=rng.choice(ng,len(ng));bidx.extend(np.concatenate([np.flatnonzero((ye==kk)&(nights==d)) for d in draw]).tolist())
                k=np.array(bidx)
                ff,nn,_=fit_centers(phase[k],ye[k],y[k],e[k],g[k])
                a,b=ff.x[nn:nn+4].reshape(2,2);vals.append(np.diff(b-a)[0]*PERIOD*1440)
            result['boot_delta_95_min']=np.percentile(vals,[2.5,50,97.5]).tolist()
        # Per-OGLE-model comparison with common affine/year Swift baselines and the same two free depths.
        phot=json.load(open(D/'inputs/ogle-timing-models.json'));dev={}
        for model,m in phot['models'].items():
            th=m['theta'];U=m.get('best_period_yr')
            if model=='constant':delta=np.full(len(year),th[2])
            elif model in ['linear','quadratic','cubic']:
                dim=['constant','linear','quadratic','cubic'].index(model)+1
                delta=sum(th[i+2]*(year/10)**i for i in range(dim))
            elif model in ['apsidal','free_offset']:
                delta=th[2]*np.cos(year*2*np.pi/U)+th[3]*np.sin(year*2*np.pi/U)
                if model=='free_offset':delta=th[2]+th[3]*np.cos(year*2*np.pi/U)+th[4]*np.sin(year*2*np.pi/U)
            else:continue
            # Orbit model uses epoch drift in cycles per 10 years.
            center=th[0]+th[1]*year/10
            c1=center-delta/2;c2=.5+center+delta/2
            u=(phase-c1+.5)%1-.5;v=(phase-c2+.5)%1-.5
            X=np.column_stack((eye,-np.exp(-.5*(u/PARAMS[4])**2),-np.exp(-.5*(v/PARAMS[5])**2)))
            p=np.linalg.lstsq(X/e[:,None],y/e,rcond=None)[0]
            sq=((y-X@p)/e)**2;dev[model]=dict(chi2=float(sq.sum()),per_year={str(yy):float(sq[ye==yy].sum()) for yy in years})
        result['ogle_forecasts']=dev;out[str(rad)]=result
    (D/'swift-motion.json').write_text(json.dumps(out,indent=2));print(out['3.5'])
if __name__=='__main__':main()
