"""Reference source ratios in good fields to flag spatially coherent LoDeSS flux loss.

Usage: python local_controls.py DATA_DIRECTORY COMPARISON_RESULTS_DIRECTORY [--targets row ...]
Ratios compare a straight spectrum through three external integrated fluxes with measured
LoDeSS flux. Intrinsic curvature and survey-selection biases contribute to the ratios,
so they are diagnostics of measurement reliability and cannot justify a flux correction.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.spatial import cKDTree


def unit_xyz(ra,dec):
    a, d = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)])


def prepare(root, matches):
    s = fits.getdata(root/'sources.fits')
    c = pd.read_csv(root/'source_fields.csv')
    c['ra_deg'] = s['RA']
    c['dec_deg'] = s['DEC']
    c['peak_snr'] = s['Peak_flux']/s['Isl_rms']
    c['flux23_jy'] = s['Total_flux']
    c['integrated_to_peak'] = s['Total_flux']/s['Peak_flux']
    sky = unit_xyz(c.ra_deg,c.dec_deg)
    near = cKDTree(sky).query(sky,k=2)[0][:,1]
    c['near_arcsec'] = np.rad2deg(2*np.arcsin(near/2))*3600
    c['simple'] = s['S_Code']=='S'
    for table in ['vlssr','lotss','nvss']:
        x = pd.read_csv(matches/f'match_{table}.csv.gz')
        c = c.merge(x,on='row',how='left',validate='one_to_one')
    obs = c[['vlssr_flux_jy','lotss_flux_jy','nvss_flux_jy']].values
    # Coefficients of an unweighted linear log-log fit of the external measurements.
    X = np.column_stack([np.ones(3), np.log([73.8,144.,1400.])])
    coef = np.array([1, np.log(23.)]) @ np.linalg.pinv(X)
    with np.errstate(invalid='ignore',divide='ignore'):
        c['ref_ratio'] = np.exp(np.log(obs)@coef)/c.flux23_jy
    return c


def controls(c):
    return c[(c.quality=='good') & c.simple & (c.peak_snr>=10)
             & (c.near_arcsec>90) & c.integrated_to_peak.between(.6,2.5)
             & (c.lotss_sep_as<=20) & (c.vlssr_sep_as<=25) & (c.nvss_sep_as<=25)
             & (c.lotss_sum75_jy <= 1.05*c.lotss_flux_jy) & (c.lotss_dc_maj<=12)
             & (c.lotss_flux_jy/c.lotss_err_jy>10) & (c.vlssr_flux_jy/c.vlssr_err_jy>6)
             & (c.nvss_flux_jy/c.nvss_err_jy>8) & np.isfinite(c.ref_ratio) & (c.ref_ratio>0)]


def nearby(ctrl, field, ra, dec, radius_deg=2., excluded_row=None):
    same = ctrl[(ctrl.field==field) & (ctrl.row!=excluded_row)].copy()
    if same.empty:
        return same
    chord = np.linalg.norm(unit_xyz(same.ra_deg,same.dec_deg)-unit_xyz([ra],[dec]),axis=1)
    same['target_dist_deg'] = np.rad2deg(2*np.arcsin(np.clip(chord/2,0.,1.)))
    return same[same.target_dist_deg<=radius_deg].sort_values('target_dist_deg')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path);p.add_argument('matches',type=Path);p.add_argument('--targets',type=int,nargs='*',default=[])
    args=p.parse_args()
    c=prepare(args.root,args.matches)
    ctrl=controls(c)
    cols=['row','field','ra_deg','dec_deg','separation_deg','flux23_jy','ref_ratio']
    ctrl[cols].to_csv(args.root/'reference_sources.csv',index=False,float_format='%.7g')
    out=[]
    for row in args.targets:
        t=c.loc[c.row==row].iloc[0]
        z=nearby(ctrl,t.field,t.ra_deg,t.dec_deg,excluded_row=row)
        r=z.ref_ratio
        out.append(dict(row=row,field=t.field,target_ratio=t.ref_ratio,
                        matched_refs=len(z),median=r.median(),q25=r.quantile(.25),q75=r.quantile(.75),
                        largest=r.max(),neighbour_radius_deg=2.0))
    if out:
        pd.DataFrame(out).to_csv(args.root/'target_controls.csv',index=False)
    print(ctrl.groupby('field').ref_ratio.agg(['size','median']))
    print(pd.DataFrame(out).round(2))
