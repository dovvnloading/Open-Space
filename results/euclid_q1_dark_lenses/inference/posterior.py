"""Conditional Monte Carlo mass-to-light distributions for the two photometric leads and spectroscopic controls.

python posterior.py MODEL_DIRECTORY PHZ_CSV
"""
import ast
import csv
import io
import json
import sys
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid
from infer import HERE, C_KMS, M_SCALE, BANDS, read, write

LEADS = ('102158887_2655306065654633189','102158888_2659753895653871052')


def run(models, phz_file):
    rng = np.random.default_rng(16032025)
    n = 48000
    reference = np.genfromtxt(HERE/'sun_shifted_band_ab.csv',delimiter=',',names=True)
    zaxis = np.linspace(0,8,64001)
    dc = C_KMS/70 * cumulative_trapezoid(1/np.sqrt(.3*(1+zaxis)**3+.7), zaxis, initial=0)
    def comoving(z):
        if np.any(np.asarray(z)>zaxis[-1]) or np.any(np.asarray(z)<0):
            raise ValueError('redshift outside cosmology grid')
        return np.interp(z,zaxis,dc)
    geometry = {r['id_str']:r for r in read(HERE/'redshift_geometry.csv')}
    radii = {r['id_str']:r for r in read(models/'modeling_lens_mass.csv')}
    q1 = {r['object_id']:r for r in read(phz_file)}
    ls = {r['id_str']:r for r in read(HERE/'legacy_matches.csv') if float(r['delta_arcsec'])<1}
    physical_rows=csv.DictReader(io.StringIO(json.loads((HERE/'lead_physical_query.json').read_text())['raw_csv']))
    pp = {r['object_id']:float(r['phz_pp_median_redshift']) for r in physical_rows}
    phot = read(HERE/'mass_light.csv')
    selected = [r for r in phot if r['id_str'] in LEADS or geometry[r['id_str']]['zs']]
    result, ranks = [], []
    for row in selected:
        key, band = row['id_str'],row['band']
        if not row['L_E_Lsun']:
            continue
        g = geometry[key]
        t = radii[key]
        theta0 = float(g['theta_eff_arcsec'])
        norm = rng.standard_normal(n)
        low = theta0 - float(t['einstein_radius_effective_lower_1_sigma'])
        high = float(t['einstein_radius_effective_upper_1_sigma']) - theta0
        theta = theta0 + np.where(norm < 0, low, high)*norm
        # Flat positive intensity prior + Gaussian likelihood for measured aperture flux.
        f_ratio = 1 + rng.normal(size=n)/float(row['SNR_E'])
        z_u = rng.uniform(.025,.975,n)
        if key in LEADS:
            p = q1[key.split('_')[1]]
            lo70, hi70 = ast.literal_eval(p['phz_70_int'])
            lo95, hi95 = ast.literal_eval(p['phz_95_int'])
            legacy = ls[key]
            scenarios = [
                ('Euclid_PHZ', (.025,.15,.5,.85,.975), [lo95,lo70,float(p['phz_median']),hi70,hi95]),
                ('Legacy_DR9', (.025,.16,.5,.84,.975), [float(legacy[k]) for k in (
                    'z_phot_l95','z_phot_l68','z_phot_median','z_phot_u68','z_phot_u95')])]
            if key == LEADS[1]:
                scenarios.append(('PHZ_PP_median_fixed',[.5],[pp[key.split('_')[1]]]))
            source_z = [1.5,2.,3.]
        else:
            scenarios = [('spectroscopic_fixed', [.5],[float(g['zl'])])]
            source_z = [float(g['zs'])]
        for label, levels, knots in scenarios:
            if np.any(np.diff(knots)<0):
                raise ValueError('Nonmonotonic redshift reference for '+key)
            zl = np.interp(z_u,levels,knots)
            cl = comoving(zl)
            dl = cl*(1+zl)*1e6
            absolute = float(row['mag_E'])-5*np.log10(dl/10)+2.5*np.log10(1+zl)
            msol = np.interp(zl,reference['z'],reference[band])
            lum = 10**(.4*(msol-absolute))*f_ratio
            for zs in ([3.,4.,6.] if label=='PHZ_PP_median_fixed' else source_z):
                mass = M_SCALE*theta**2*cl/(1+zl)/(1-cl/comoving(zs))
                ok = (theta>0)&(f_ratio>0)&(zl<zs)
                ml = mass[ok]/lum[ok]
                q = np.percentile(ml,[2.5,16,50,84,97.5])
                result.append(dict(id_str=key,band=band,z_l_model=label,zs=zs,
                                   n=n,n_physical=int(np.sum(ok)),M_E_med_Msun=np.median(mass[ok]),
                                   L_E_med_Lsun=np.median(lum[ok]),ML_p025=q[0],ML_p16=q[1],
                                   ML_p50=q[2],ML_p84=q[3],ML_p975=q[4]))
                if key not in LEADS or band!='h' or zs!=2:
                    continue
                # Other grade A/B systems with robust redshift inputs; use same band and source plane.
                peers=[]
                for other in phot:
                    i=other['id_str']
                    if i in LEADS or other['band']!='h' or not other['ML_E_beta1']:
                        continue
                    p2=geometry[i]
                    z=float(p2['zl'])
                    if p2['grade'] not in ('A','B') or not (0 < z < zs):
                        continue
                    if p2['zl_basis']=='merged_PHZ':
                        rz=q1[i.split('_')[-1].replace('NEG','-')]
                        bounds=ast.literal_eval(rz['phz_95_int'])
                        if rz['phz_flags']!='0' or (bounds[1]-bounds[0])/(1+z)>.3:
                            continue
                    peers.append((abs(z-knots[len(knots)//2]),z,
                                  float(other['ML_E_beta1'])/(1-comoving(z)/comoving(zs))))
                near=sorted(peers)[:25]
                v=np.array([x[2] for x in near])
                med=np.median(v)
                ranks.append(dict(id_str=key,z_l_model=label,zs=zs,n_peers=len(v),
                                  peer_z_min=min(x[1] for x in near),peer_z_max=max(x[1] for x in near),
                                  peer_ML_p50=med,peer_ML_p95=np.percentile(v,95),peer_ML_max=np.max(v),
                                  candidate_ML_p50=q[2],n_peers_at_or_above=int(np.sum(v>=q[2])),
                                  A_H_to_peer_median_mag=2.5*np.log10(q[2]/med)))
    write(HERE/'conditional_ML_posterior.csv',result)
    write(HERE/'redshift_matched_peers.csv',ranks)

if __name__=='__main__':
    run(Path(sys.argv[1]),Path(sys.argv[2]))
