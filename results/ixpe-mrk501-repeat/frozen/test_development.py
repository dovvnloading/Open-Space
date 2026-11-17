"""Independent synthetic source plus background checks; no observation arrays."""
import numpy as np
from development import estimate, CAL, particle_keep
from ixpeobssim.irf.arf import xEffectiveArea
from ixpeobssim.irf.modf import xModulationFactor


def synthetic(seed, expected):
    rng=np.random.default_rng(seed)
    arf=xEffectiveArea(str(CAL/'arf'/'ixpe_d1_obssim20260701_v013.arf'))
    mod=xModulationFactor(str(CAL/'modfact'/'ixpe_d1_obssim20260701_mfact_v013.fits'))
    def events(n, qu):
        e=rng.uniform(2.02,7.98,n);mu=mod(e);area=arf(e)
        phi=rng.uniform(0,2*np.pi,n);result=np.zeros((n,3));remaining=np.ones(n,dtype=bool)
        while remaining.any():
            k=remaining.sum();phi[remaining]=rng.uniform(0,2*np.pi,k)
            a=phi[remaining];m=mu[remaining]
            prob=(1+m*(qu[0]*np.cos(2*a)+qu[1]*np.sin(2*a)))/(1+m*np.hypot(*qu))
            remaining[np.flatnonzero(remaining)[rng.random(k)<prob]]=False
        result[:,0]=1/area;result[:,1]=2*np.cos(2*phi)/(area*mu);result[:,2]=2*np.sin(2*phi)/(area*mu)
        return result
    source=np.concatenate([events(250000,expected),events(20000,(0,0))]);background=events(200000,(0,0))
    z,c,_=estimate(source,background,0.1);measured=z[1:]/z[0]
    print('expected',expected,'measured',measured,'sigma',np.sqrt(np.diag(c)))
    assert np.max(np.abs(measured-expected))<0.05
    assert np.linalg.eigvalsh(c).min()>0

def coverage():
    rng=np.random.default_rng(9127)
    # Unit-area and constant mu, with independent Poisson source/background counts.
    # Validates joint covariance calibration, not spacecraft background behavior.
    covered=0; checks=[]; mu=0.4
    for _ in range(400):
        ns=rng.poisson(5000); nb1=rng.poisson(600); nb2=rng.poisson(6000)
        def draw(n):
            phi=rng.uniform(0,2*np.pi,n)
            return np.column_stack((np.ones(n),2*np.cos(2*phi)/mu,2*np.sin(2*phi)/mu))
        src=np.concatenate((draw(ns),draw(nb1))); bg=draw(nb2)
        z,c,_=estimate(src,bg,0.1); v=z[1:]/z[0]
        chi=float(v@np.linalg.solve(c,v)); checks.append(chi); covered+=chi<5.991464547
    fraction=covered/len(checks)
    print('400 null realizations: 95% joint ellipse coverage',fraction)
    assert 0.90 < fraction < 0.98

if __name__=='__main__':
    flags=np.zeros((3,16),dtype=bool);flags[0,7]=True;flags[1,8]=True
    assert particle_keep({'XPFLGBGD':True},flags).tolist()==[False,True,True]
    try:particle_keep({},flags)
    except ValueError:pass
    else:raise AssertionError('missing particle flag must fail closed')
    synthetic(764,(0,0));synthetic(765,(0.15,0.2598076211));coverage()
