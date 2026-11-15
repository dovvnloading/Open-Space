"""Measure the orbital clock of KIC 8110754 using the resolved ZTF source alone."""
import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from ground_photometry import read_star, wrap

HERE = Path(__file__).resolve().parent
T_REF = 2459250.0  # BJD TDB; near the center of the observations
P_REF = 0.518059
EPOCH_GUESS = 2454833+131.9425
T_REF += wrap((EPOCH_GUESS-T_REF)/P_REF)*P_REF


def fit(label, use, time, band, oid, mag):
    t, band, oid, y = time[use], band[use], oid[use], mag[use]
    bands, b = np.unique(band, return_inverse=True)
    groups, g = np.unique(oid, return_inverse=True)
    # Every field/filter series has its own zero point. Set a robust magnitude
    # noise scale per band using differences between points close in orbital phase.
    scale = np.ones(len(t))
    for j in range(len(bands)):
        subset = np.flatnonzero(b==j)
        phase = (t[subset]-T_REF)/P_REF % 1
        idx = subset[np.argsort(phase)]
        scale[idx] = max(.025, 1.4826*np.median(abs(np.diff(y[idx])))/np.sqrt(2))
    offsets=[np.median(y[g==j]) for j in range(len(groups))]
    base=[0., 0., 0., .033, .038] + offsets
    v0=base + [.55, .23, 0., 0.] * len(bands)
    lower=[-.06,-1.,-.06,.012,.012]+[-1.]*len(groups)+[.1,.025,-.25,-.25]*len(bands)
    upper=[.06,1.,.06,.09,.12]+[1.]*len(groups)+[1.9,.85,.25,.25]*len(bands)
    for j in range(len(groups)):
        y[g==j] -= offsets[j]
        v0[5+j] = 0.
    def residual(v):
        # Center offset in days and period offset in seconds improve conditioning.
        p = P_REF+v[1]/86400
        phi = wrap((t-T_REF-v[0])/p)
        shape = v[5+len(groups):].reshape(-1,4)[b]
        model = (v[5:5+len(groups)][g]
                 + shape[:,0]*np.exp(-.5*(phi/v[3])**2)
                 + shape[:,1]*np.exp(-.5*(wrap(phi-.5-v[2])/v[4])**2)
                 + shape[:,2]*np.cos(2*np.pi*phi) + shape[:,3]*np.cos(4*np.pi*phi))
        return (model-y)/scale
    solution=least_squares(residual,v0,bounds=(lower,upper),loss='soft_l1',f_scale=1.5,
                          x_scale='jac',max_nfev=300)
    # Cluster-robust covariance for soft-L1 estimating equations. Group all
    # measurements on the same Palomar night, including both filters.
    from scipy.optimize._numdiff import approx_derivative
    v=solution.x
    jac=approx_derivative(residual,v,abs_step=np.repeat(1e-5,len(v)))
    r=residual(v); frac=(1+(r/1.5)**2)
    projected=np.linalg.pinv(jac.T@(jac*frac[:,None]**-1.5))
    nights,night=np.unique(np.floor(t),return_inverse=True)
    scores=np.zeros((len(nights),len(v)))
    np.add.at(scores,night,jac*(r/np.sqrt(frac))[:,None])
    cov=projected@(scores.T@scores)@projected
    return dict(selection=label,n=len(t),nights=len(nights),period_days=P_REF+v[1]/86400,
        period_se_seconds=np.sqrt(cov[1,1]),primary_epoch_bjd_tdb=T_REF+v[0],
        primary_epoch_se_minutes=np.sqrt(cov[0,0])*1440,secondary_phase=.5+v[2],
        first_bjd_tdb=t.min(),last_bjd_tdb=t.max(),success=solution.success)


if __name__=='__main__':
    t,band,oid,mag = read_star('source')
    selections={'all':np.ones(len(t),bool),'r':band=='zr','g':band=='zg',
                '2018-2021':t<2459580.5,'2022-2025':t>=2459580.5}
    rows=[fit(k,z,t,band,oid,mag) for k,z in selections.items()]
    with (HERE/'ztf_periods.csv').open('w') as handle:
        w=csv.DictWriter(handle,rows[0].keys(),lineterminator='\n');w.writeheader();w.writerows(rows)
