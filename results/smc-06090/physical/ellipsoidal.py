"""Compare a held-out orbital phase harmonic with the leading equilibrium-tide estimate."""
import json
import numpy as np
from scipy.optimize import least_squares
from geometry import HERE,P,EPOCH,TREF,load,sky

def wrap(ph):
    return (ph+.5)%1-.5
def robust_regression(B,mag,err):
    W=B/err[:,None]
    return least_squares(lambda v:(B@v-mag)/err,np.linalg.lstsq(W,mag/err,rcond=None)[0],
                         jac=lambda v:W.copy(),loss='soft_l1',f_scale=1.8,ftol=1e-7).x

def measure(data,q,phase_exclusion,order=2,template="circular"):
    h,mag,err,grp=data
    t=(h-TREF)/365.25; ph=(h-EPOCH)/P
    d,_=sky(t,ph%1,q,1)
    orbital_phase=ph-q[0]-q[1]*t/10
    mask=(d>q[5]*phase_exclusion)
    # In magnitudes an ellipsoidal tidal signal has a positive conjunction cosine.
    angle=2*np.pi*orbital_phase[mask]
    B=np.stack([grp[mask]==0,grp[mask]==1,np.cos(angle),np.sin(angle),np.cos(2*angle),np.sin(2*angle)],axis=-1)
    if template=='eccentric':
        # Linear, quasi-static response to the changing separation; eccentric orbit in the sky.
        om=q[4]+2*np.pi*t[mask]/q[2]
        mean=angle+np.pi/2-om
        E=mean.copy()
        for _ in range(5): E-=(E-q[3]*np.sin(E)-mean)/(1-q[3]*np.cos(E))
        cosf=(np.cos(E)-q[3])/(1-q[3]*np.cos(E))
        sinf=np.sqrt(1-q[3]**2)*np.sin(E)/(1-q[3]*np.cos(E))
        lon=om+np.arctan2(sinf,cosf)
        B[:,4]=-np.cos(2*lon)/(1-q[3]*np.cos(E))**3
        B[:,5]=-np.sin(2*lon)/(1-q[3]*np.cos(E))**3
    if order>2:
        harmonics=[ff(k*angle) for k in range(3,order+1) for ff in (np.cos,np.sin)]
        B=np.column_stack((B,*harmonics))
    series=mag[mask];weight=err[mask]
    vec=robust_regression(B,series,weight)
    nnight=np.floor(h[mask]-.5)
    blocks=np.digitize(h[mask],[2453500,2455000,2457000,2458300,2459700])
    groups=[[np.where((blocks==g)&(nnight==n))[0] for n in np.unique(nnight[blocks==g])] for g in range(6)]
    rng=np.random.default_rng(560026)
    draws=[]
    for _ in range(600):
        p=np.concatenate([ns[j] for ns in groups for j in rng.integers(len(ns),size=len(ns))])
        draws.append(robust_regression(B[p],series[p],weight[p])[4:6]*1000)
    injected=(robust_regression(B,series+.006*B[:,4],weight)[4]-vec[4])*1000
    era=[1000*robust_regression(B[g,:][:,[gnum,2,3,4,5]],series[g],weight[g])[3]
         for gnum in (0,1) for g in [grp[mask]==gnum]]
    return dict(n=len(angle),injected_mmag=6.,recovered_mmag=injected,
                separate_survey_cos2_mmag=era,order=order,template=template,night_count=len(np.unique(nnight)),margin=phase_exclusion,
                cos2_mmag=1000*vec[4],sin2_mmag=1000*vec[5],
                conditional_ci95=np.percentile(draws,[2.5,97.5],axis=0).T.tolist())
def expectation(q,third,u,tau):
    rs,k,ci,J=q[5:9]
    ra=rs/(1+k); rb=ra*k
    ratio=J*k*k
    # The unconstrained mass ratio q_M>0 minimizing A*q_M+B/q_M is sqrt(B/A).
    mt=np.sqrt(ratio*rb**3/ra**3)
    alpha=.15*(15+u)*(1+tau)/(3-u)
    def amplitude(qm):
        frac=(1-third)*alpha*(1-ci**2)*(ra**3*qm+ratio*rb**3/qm)/(1+ratio)
        return 1000*2.5/np.log(10)*frac
    return dict(u=u,gravity_darkening=tau,third=third,minimizing_mass_ratio=mt,
                minimum_cos2_mmag=amplitude(mt),equal_mass_cos2_mmag=amplitude(1.))
def main():
    data=load(HERE/'inputs'); fits=json.loads((HERE/'dilution.json').read_text())['1']
    out={}
    for row in fits:
        third=row['third'];q=np.asarray(row['pars'])
        if third not in (0.,.3):continue
        out[str(third)]=dict(empirical=[measure(data,q,margin) for margin in (1.,1.18,1.38)],
                nuisance_harmonics=[measure(data,q,1.,order=o,template=tp) for o,tp in [(2,'eccentric'),(3,'circular'),(4,'circular')]],
                theoretical=[expectation(q,third,u,tau) for u in (.2,.4,.6) for tau in (0.,.32,1.)])
    (HERE/'ellipsoidal.json').write_text(json.dumps(out,indent=2)+'\n')
    for k,r in out.items():print(k,[(x['n'],round(x['cos2_mmag'],3),np.round(x['conditional_ci95'][0],2)) for x in r['empirical']],
                            min(x['minimum_cos2_mmag'] for x in r['theoretical']),max(x['minimum_cos2_mmag'] for x in r['theoretical']))
if __name__=='__main__':main()
