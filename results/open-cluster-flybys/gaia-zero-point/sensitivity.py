"""Replay the star-by-star Gaia parallax and membership sensitivity check."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
# Jacobi radii in pc from Hunt & Reffert 2024, J/A+A/686/A42.
PAIRS = [
    ('ascc114', 'ubc168', 12.32+10.16), ('upk90','upk93',6.49+6.70),
    ('ascc99','upk20',7.23+7.52), ('ngc6633','upk41',11.193856+9.108565)
]

def vectors(ra, dec, parallax, pmra, pmdec, rv):
    a, b = np.deg2rad(ra), np.deg2rad(dec)
    r = np.stack([np.cos(a)*np.cos(b), np.sin(a)*np.cos(b), np.sin(b)], axis=-1)
    east = np.stack([-np.sin(a), np.cos(a), np.zeros_like(a)], axis=-1)
    north = np.stack([-np.cos(a)*np.sin(b),-np.sin(a)*np.sin(b),np.cos(b)],axis=-1)
    return (r * (1000/parallax)[:,None], r * rv[:,None] +
            4.74047*(east*pmra[:,None] + north*pmdec[:,None])/parallax[:,None])

def sample(g, floor, rng, nboot, probability, use_zpt, max_ruwe):
    g = g[g.probs >= probability]
    q = g[g.radial_velocity_error < 12]
    v = q.radial_velocity
    keep = abs(v - v.median()) < max(6, 4.5 * (v-v.median()).abs().median())
    rv = v[keep].to_numpy()
    m = g[(g.parallax_error < .18) & (g.ruwe < max_ruwe) & g.zpt_mas.notna()]
    parallax = m.parallax.to_numpy() - (m.zpt_mas.to_numpy() if use_zpt else -.017)
    ix = rng.integers(len(m),size=(nboot,len(m)))
    iz = rng.integers(len(rv),size=(nboot,len(rv)))
    def med(x):
        return np.median(np.asarray(x)[ix],axis=1)
    return (vectors(med(m.ra),med(m.dec),med(parallax)+rng.normal(0,floor,nboot),
                   med(m.pmra),med(m.pmdec),np.median(rv[iz],axis=1)),len(m),len(rv))

def replay(data, pairs=PAIRS, draws=12000):
    records=[]
    for a,b,radius in pairs:
        for label,p,floor,corrected,ruwe in [
            ('individual',.5,.010,True,1.4),('fixed_offset',.5,.010,False,1.4),
            ('p80',.8,.010,True,1.4),('20uas',.5,.020,True,1.4),
            ('relaxed_ruwe',.5,.010,True,np.inf)]:
            rng=np.random.default_rng(20261114)
            (x,v),na,va=sample(data[data.name==a],floor,rng,draws,p,corrected,ruwe)
            (y,w),nb,vb=sample(data[data.name==b],floor,rng,draws,p,corrected,ruwe)
            r,u=x-y,v-w
            tau=-np.sum(r*u,axis=1)/np.sum(u*u,axis=1)
            dist=np.linalg.norm(r+tau[:,None]*u,axis=1)
            row=dict(cluster_a=a,cluster_b=b,variant=label,ast_a=na,ast_b=nb,rv_a=va,rv_b=vb,
                     radius_pc=radius,fraction_inside=np.mean(dist<radius))
            for name,col in [('t_myr',tau/1.022712165),('d_pc',dist),('speed_kms',np.linalg.norm(u,axis=1))]:
                for suf,val in zip(['p16','p50','p84'],np.percentile(col,[16,50,84])):
                    row[f'{name}_{suf}']=val
            records.append(row)
    return pd.DataFrame(records)

if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument('--members',type=Path,default=ROOT/'members.csv.gz')
    args.add_argument('--output',type=Path,default=ROOT/'sensitivity.csv')
    args=args.parse_args()
    replay(pd.read_csv(args.members)).to_csv(args.output,index=False,float_format='%.6f')
