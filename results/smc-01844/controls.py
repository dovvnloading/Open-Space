"""Refit external-survey phases and challenge the local OGLE fit and its near harmonics."""
import json

import numpy as np
from scipy.optimize import minimize_scalar, minimize
import replay

freq=4.912246653669684


def main():
    _,t,y,e,h=next(replay.series())
    r={}
    for n in [8,11,15]:
        opt=minimize_scalar(lambda f:replay.fit(t,y,e,f,n)[2],bounds=replay.SEARCH_BOUNDS,
                            method='bounded',options={'xatol':1.e-10})
        r[n]={'freq':opt.x,'peak_to_peak_mag':2*np.linalg.norm(replay.fit(t,y,e,opt.x,n)[0][-2:])}
    ph=(t-replay.EPOCH)/replay.BINARY_PERIOD % 1
    distances=np.abs((ph+.25)%.5-.25)
    out=distances>.12
    # Broad intervals around both eclipses are excised; use the same model order throughout.
    ref=replay.fit(t,y,e,freq)[0][-2:];u=ref/np.linalg.norm(ref)
    b,c,_,_=replay.fit(t[out],y[out],e[out],freq,6)
    off={'n':sum(out).item(),'semiamplitude_mag':u@b[-2:],'se':np.sqrt(u@c[-2:,-2:]@u)}
    # Keep the time-series data and finite sampling fixed as alternate daily alias clocks are refit.
    m0=replay.fit(t,y,e,freq,11)[2]
    alt={}
    for name,start in [('sidereal-day minus',freq-1.00273791),('sidereal-day plus',freq+1.00273791),
                        ('solar-day minus',freq-1),('solar-day plus',freq+1),('9th orbital harmonic',9/replay.BINARY_PERIOD)]:
        if name.startswith('9th'):
            # At exactly the harmonic the design already contains the tested sinusoid.
            test=start
        else:
            opt=minimize_scalar(lambda f:replay.fit(t,y,e,f,11)[2],
                    bounds=(start-.00015,start+.00015),method='bounded',
                    options={'xatol':1.e-10})
            test=opt.x
        alt[name]={'freq':test,'robust_loss_increase':replay.fit(t,y,e,test,11)[2]-m0}
    original_period = replay.BINARY_PERIOD
    def loss(delta):
        replay.BINARY_PERIOD=original_period+delta[0]*1.e-6
        return replay.fit(t,y,e,freq+delta[1]*1.e-5,11)[2]
    joint=minimize(
        loss,[0.,0.],method='Nelder-Mead',options={'initial_simplex':np.array([[0.,0.],[1.,0.],[0.,.5]]),'xatol':1.e-4})
    replay.BINARY_PERIOD=original_period
    print(json.dumps(dict(orders=r,outside_eclipse=off,aliases=alt,joint_fit=dict(binary_period_d=original_period+joint.x[0]*1.e-6,
               short_frequency_per_d=freq+joint.x[1]*1.e-5)),indent=2))


if __name__=='__main__':
    main()
