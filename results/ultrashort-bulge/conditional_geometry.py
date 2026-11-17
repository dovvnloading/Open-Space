"""Conditional Galactic masses for a specified source identity and lens geometry.

Integrate over the two-dimensional geocentric velocity density at a supplied
Einstein time, which remains numerically stable even for slow proper motions.
Rows for unresolved sources are sensitivity calculations, not mass measurements.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path

import numpy as np
from mass_controls import EARTHS_PER_SUN, HERE, fit, galaxy


def calculate(ra, dec, hjd, te, theta_uas, seed=312741, n=300000,
              mass_limits=(.0001,1e6)):
    rng=np.random.default_rng(seed)
    ra,dec=fit.lens.sexagesimal(ra,dec)
    l,b,angle_l,angle_b=galaxy.galactic(ra,dec)
    ds=rng.uniform(5.,12.,n)
    dl=rng.uniform(.15,ds,n)
    disk,bar=galaxy.density(dl,l,b)
    source=galaxy.density(ds,l,b)[1]
    isdisk=rng.random(n)<disk/(disk+bar)
    sl,sb=galaxy.mean_velocity(ds,l,b,0.)
    vl,vb=galaxy.mean_velocity(dl,l,b,isdisk.astype(float))
    conv=galaxy.KMS_PER_MASYR_KPC
    pil=1/dl-1/ds
    _,_,sn,se=fit.lens.sun_projected(np.array([hjd-2450000]),ra,dec)
    def projected(angle):
        return -365.25*(sn[0]*np.cos(angle)+se[0]*np.sin(angle))
    mean_l=(vl/dl-sl/ds)/conv-pil*projected(angle_l)
    mean_b=(vb/dl-sb/ds)/conv-pil*projected(angle_b)
    sig_l=np.hypot(np.where(isdisk,34.,100.)/dl,100./ds)/conv
    sig_b=np.hypot(np.where(isdisk,18.,100.)/dl,100./ds)/conv
    theta_mas=theta_uas*.001
    mu=365.25*theta_mas/te
    angle=rng.uniform(0.,2*np.pi,n)
    power=((mu*np.cos(angle)-mean_l)/sig_l)**2+((mu*np.sin(angle)-mean_b)/sig_b)**2
    velocity_density=np.exp(-.5*power)/(sig_l*sig_b)
    mass=theta_mas**2*EARTHS_PER_SUN/(galaxy.KAPPA*pil)
    # For fixed theta and tE, dGamma/dtheta/dtE has factor 2*mu^3/tE:
    # theta*mu from the rate, 2/theta from dlogM, and mu^2/tE from d^2mu.
    # That factor is common to all distance draws within a scenario.
    w=source*ds**2*(ds-.15)*(disk+bar)*dl**2*velocity_density
    w[(mass<mass_limits[0])|(mass>mass_limits[1])]=0.
    out=[]
    for slope in [0.,-.96]:
        weights=w*mass**slope
        mass_q=galaxy.quantiles(mass,weights)
        dl_q=galaxy.quantiles(dl,weights)
        out.append(dict(dN_dlogM_slope=slope,mu_geo_mas_yr=mu,
                        mass_earth_16=mass_q[0],mass_earth_50=mass_q[1],
                        mass_earth_84=mass_q[2],dl_kpc_16=dl_q[0],
                        dl_kpc_50=dl_q[1],dl_kpc_84=dl_q[2],
                        effective_draws=round(weights.sum()**2/(weights@weights)),
                        mass_prior_min_earth=mass_limits[0],mass_prior_max_earth=mass_limits[1]))
    return out


def source_cases():
    catalog={r['Event']:r for r in csv.DictReader((HERE/'kmt/catalog-2024.tsv').open(),delimiter='\t')}
    compact_fits=json.loads(gzip.decompress((HERE/'tables/peak-comparisons.json.gz').read_bytes()))
    with gzip.open(HERE/'tables/continuous-fits.csv.gz','rt') as f:
        native={r['event']:r for r in csv.DictReader(f) if r['seed_pass']=='capped KMT grid seed'}
    for row in json.loads((HERE/'kmt/extended/angular.json').read_text()):
        n=row['event']
        short=compact_fits[n]['initial']['models']['fs']
        geometries=[('peak-window',short['t0'],short['te_or_rise_d'],short['rho']),
                    ('full-curve',native[n]['fs_t0'],native[n]['fs_te'],native[n]['rho'])]
        for obj in row['source_options']:
            for label,t0,te,rho in geometries:
                yield dict(event=n,assumed_VIRAC2_source=obj['srcid'],offset_arcsec=obj['sep_arcsec'],
                           source_hypothesis='all J/Ks flux; dwarf or subgiant; red-clump extinction',
                           theta_star_uas=obj['theta_star_uas_nominal'],geometry=label,
                           ra=catalog[n]['RA'],dec=catalog[n]['Dec'],t0_jd=t0,te_d=te,rho=rho)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--scenarios',type=Path,help='optional user-supplied CSV using the output columns')
    p.add_argument('--out',type=Path,default=HERE/'tables/conditional-source-masses.csv')
    a=p.parse_args()
    if a.scenarios:
        with a.scenarios.open() as f:cases=list(csv.DictReader(f))
    else:
        cases=list(source_cases())
    rows=[]
    for case in cases:
        theta=float(case['theta_star_uas'])/float(case['rho'])
        rows.extend({**case,'theta_e_uas':theta,**values} for values in
                    calculate(case['ra'],case['dec'],float(case['t0_jd']),float(case['te_d']),theta))
    if rows:
        with a.out.open('w',newline='') as f:
            w=csv.DictWriter(f,list(rows[0]),lineterminator='\n')
            w.writeheader();w.writerows(rows)
