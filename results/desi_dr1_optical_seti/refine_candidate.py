"""Independent common-window LSF refit of a v1 trigger, for follow-up only.

This module does not change the locked survey selection. In one fixed data
window it fits every neighbouring subpixel interval with the same baseline.
"""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2


def refit(e, seed_bary_A):
    w=np.asarray(e.wave,float)
    n=len(w)
    m=e.resolution.shape[0]
    h=m//2
    center=int(np.argmin(abs(w-seed_bary_A)))
    radius=h+5
    if center-radius<0 or center+radius>=n:
        return None
    sl=slice(center-radius,center+radius+1)
    y=np.asarray(e.flux[sl],float)-np.asarray(e.template[sl],float)
    weight=np.asarray(e.ivar[sl],float)
    if not (np.all(np.asarray(e.mask)[sl]==0) and np.all(weight>0) and
            np.all(np.isfinite(y)) and np.all(np.isfinite(weight))):
        return None
    ws=w[sl]
    root=np.sqrt(weight)
    baseline=np.column_stack((np.ones(len(ws)),ws-w[center]))*root[:,None]
    q,_=np.linalg.qr(baseline,mode='reduced')
    def projected(p):
        weighted=p*root
        return weighted-q@(q.T@weighted)
    yy=projected(y)
    widths=np.gradient(w)
    def column(j):
        p=np.zeros(len(ws))
        start=j-h-(center-radius)
        c=np.asarray(e.resolution[:,j],float)
        area=c@widths[j-h:j+h+1]
        if not np.isfinite(area) or area<=0:
            raise ValueError('invalid resolution profile')
        p[start:start+m]=c/area
        return projected(p)
    columns={j:column(j) for j in range(center-1,center+3)}
    best=None
    for j in range(center-1,center+2):
        p=columns[j]
        d=columns[j+1]-p
        def rss(t):
            model=p+t*d
            amp=max(0.,yy@model/(model@model))
            return np.sum((yy-amp*model)**2)
        solution=minimize_scalar(rss,method='bounded',bounds=(0.,1.),options={'xatol':1e-12})
        for phase in (0.,solution.x,1.):
            loss=rss(phase)
            if best is None or loss<best[0]:
                best=(loss,j,phase)
    loss,j,phase=best
    p=columns[j]+phase*(columns[j+1]-columns[j])
    sigma=1/np.linalg.norm(p)
    flux=(yy@p)*sigma**2
    return dict(bary_A=float(w[j]+phase*(w[j+1]-w[j])),flux17=float(flux),snr=float(flux/sigma),
                chi2=float(loss),dof=len(ws)-4,shape_p=float(chi2.sf(loss,len(ws)-4)))
