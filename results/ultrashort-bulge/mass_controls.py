"""Event-rate-weighted Galactic posterior for two published measured-radius controls.

This illustrates the mass-distance degeneracy. It assumes that the lens belongs
to the chosen mass range and gives no probability of planetary or unbound status.
"""
import json
from pathlib import Path
import numpy as np
import fit
import galaxy  # density, projected velocities and quantiles from ogle-dark-lenses

HERE=Path(__file__).resolve().parent
EARTHS_PER_SUN=332946.


def posterior(event,ra,dec,hjd,te,ste,theta_uas,stheta_uas,seed=112917,n=500000):
    rng=np.random.default_rng(seed)
    ra,dec=fit.lens.sexagesimal(ra,dec)
    l,b,pa_l,pa_b=galaxy.galactic(ra,dec)
    ds=rng.uniform(5.,12.,n);dl=rng.uniform(.15,ds,n)
    disk,bar=galaxy.density(dl,l,b)
    src_bar=galaxy.density(ds,l,b)[1]
    isdisk=rng.random(n)<disk/(disk+bar)
    sl,sb=galaxy.mean_velocity(ds,l,b,0.)
    vl,vb=galaxy.mean_velocity(dl,l,b,isdisk.astype(float))
    sl+=rng.normal(0.,100.,n);sb+=rng.normal(0.,100.,n)
    vl+=rng.normal(0.,np.where(isdisk,34.,100.));vb+=rng.normal(0.,np.where(isdisk,18.,100.))
    mu_l=(vl/dl-sl/ds)/galaxy.KMS_PER_MASYR_KPC
    mu_b=(vb/dl-sb/ds)/galaxy.KMS_PER_MASYR_KPC
    pi_rel=1/dl-1/ds
    _,_,sn,se=fit.lens.sun_projected(np.array([hjd-2450000]),ra,dec)
    earth_n,earth_e=-365.25*sn[0],-365.25*se[0]
    mu=np.hypot(mu_l-pi_rel*(earth_n*np.cos(pa_l)+earth_e*np.sin(pa_l)),
                mu_b-pi_rel*(earth_n*np.cos(pa_b)+earth_e*np.sin(pa_b)))
    theta_mas=rng.normal(theta_uas,stheta_uas,n)*.001
    mass=theta_mas**2/(galaxy.KAPPA*pi_rel)*EARTHS_PER_SUN
    te_pred=365.25*theta_mas/np.maximum(mu,1e-12)
    # dGamma ∝ n(Dl) Dl² dDl × theta_E mu × dN/dlogM. Drawing theta
    # from its likelihood instead of drawing logM supplies |dlogM/dtheta|=2/theta.
    weight=(ds-.15)*(ds**2)*src_bar*(dl**2)*(disk+bar)*mu*2.
    weight*=np.exp(-.5*((te_pred-te)/ste)**2)
    weight[(theta_mas<=0)|(mass<.03)|(mass>3300)]=0.
    out={'event':event,'published_input':{'te_d':te,'te_sigma_d':ste,'theta_uas':theta_uas,
        'theta_sigma_uas':stheta_uas},'draws':n,'seed':seed,'prior_mass_range_earth':[.03,3300],'priors':[]}
    for slope in (0.,-.96):
        w=weight*mass**slope
        out['priors'].append({'dN_dlogM_slope':slope,'mass_earth_q16_q50_q84':galaxy.quantiles(mass,w).tolist(),
              'lens_kpc_q16_q50_q84':galaxy.quantiles(dl,w).tolist(),
              'effective_draws':round(w.sum()**2/(w@w)),'disc_probability':float(w[isdisk].sum()/w.sum())})
    return out


if __name__=='__main__':
    a=posterior('OGLE-2016-BLG-1928','18:01:31.25','-29:07:46.2',2457557.83,.0288,.0020,.842,.064)
    b=posterior('KMT-2024-BLG-3237','17:50:31.23','-32:15:48.89',2460541.90,.54,.02,6.30,.48)
    (HERE/'controls/conditional-masses.json').write_text(json.dumps([a,b],indent=2)+'\n')
