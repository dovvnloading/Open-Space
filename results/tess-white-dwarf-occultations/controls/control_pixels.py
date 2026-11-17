"""Test the last 120-s sector for each control on independent 200-s FFI extractions."""
import argparse
import hashlib
import json
import warnings
from pathlib import Path
from urllib.parse import urlencode
import numpy as np
from astropy.io import fits
from astropy.time import Time
from ffi_acquire import download, tsv_rows, propagated, vizier, TESSCUT, GAIA_FIELDS
from ffi_photometry import load_pixel_tools, gaia_t, deposit
from pixel_evidence import measure_pixels
from wd_search import search

TARGETS = {
    267574918: ('2146576589564898688', 86, 284.4167, 53.5089),
    219244444: ('4788741548375134336', 97, 65.2735, -48.6558),
    60040774: ('2943496741862496256', 108, 94.98496, -17.9723),
}

def run(tic, meta, prior, src, cache, px):
    gid,sector,ra,dec=meta
    field=src/f'control-gaia-{tic}.tsv'
    if not field.exists():
        receipt=field.with_suffix('.tsv.json')
        fallback=vizier('I/355/gaiadr3',GAIA_FIELDS,**{'-c':f'{ra} {dec}','-c.rs':'130'})
        download(field, json.loads(receipt.read_text())['url'] if receipt.exists() else fallback)
    rows=tsv_rows(field)
    ffp=cache/f'ffictrl-{tic}-{sector}.fits'
    download(ffp,TESSCUT+urlencode(dict(ra=ra,dec=dec,sector=sector,x=9,y=9,units='px')))
    if not ffp.exists():
        raise RuntimeError(f'could not fetch {ffp}')
    header=json.loads((ffp.with_suffix('.fits.json')).read_text())
    if hashlib.sha256(ffp.read_bytes()).hexdigest()!=header['sha256']:
        raise ValueError('cached FFI hash differs from receipt')
    with fits.open(ffp) as h:
        if h[0].header['SECTOR']!=sector:
            raise ValueError('wrong FFI sector')
        t, data, quality=h[1].data['TIME'], h[1].data['FLUX'], h[1].data['QUALITY']
        year=Time(np.nanmedian(t)+2457000,format='jd',scale='tdb').jyear
        target=next(r for r in rows if r['Source']==gid)
        xy=px.target_pixel(h[1].header,*propagated(target,year))
        yy,xx=np.indices(data.shape[1:])
        ap=np.hypot(xx-xy[0],yy-xy[1])<=1.5
        edge=(xx==0)|(yy==0)|(xx==xx.max())|(yy==yy.max())
        bg=np.nanmedian(data[:,edge],axis=1)
        f=np.sum(data[:,ap],axis=1)-ap.sum()*bg
        q=(quality==0)&np.isfinite(t+f)
        cadence=round(h[1].header['TIMEDEL']*86400)
        if cadence!=200:
            raise ValueError('unexpected FFI cadence')
        res=search(t[q],f[q],cadence)
        search_match=None
        # Apply the same one-template accumulated drift and phase tolerances to
        # test a separately measured, near-contemporaneous two-minute ephemeris.
        p0,e0=prior['matching_peak']['period_days'],prior['matching_peak']['epoch_btjd']
        for rank,r in enumerate(res['periodic'],1):
            wid=r['duration_min']/1440.
            e=r['epoch_btjd']+round((np.nanmedian(t[q])-r['epoch_btjd'])/r['period_days'])*r['period_days']
            if abs(r['period_days']-p0)*np.ptp(t[q])/p0<=wid and abs((e-e0+p0/2)%p0-p0/2)<=wid:
                search_match=dict(rank=rank,**r)
                break
        ref=prior['matching_peak']
        pixels=measure_pixels(t,data,quality,xy,ref['epoch_btjd'],ref['duration_min'],ref['period_days'])
        # Finite cone and colour proxy. A wide, non-gaussian PRF can carry flux
        # from a bright source far beyond this simple width range.
        positions=[(r,px.target_pixel(h[1].header,*propagated(r,year))) for r in rows if r['Gmag'] and r['Source']!=gid]
        fractions=[]; primary_t=gaia_t(target)
        for sigma in (.65,.85,1.05):
            thru=deposit(xx[ap],yy[ap],*xy,sigma)
            ratio=sum(10**(-.4*(gaia_t(r)-primary_t))*deposit(xx[ap],yy[ap],*nxy,sigma)/thru for r,nxy in positions)
            fractions.append(float(1/(1+ratio)))
        nearest=[]
        for r,nxy in sorted(positions,key=lambda z: np.linalg.norm(z[1]-xy))[:3]:
            nearest.append(dict(gaia_id=r['Source'],sep_pix=float(np.linalg.norm(nxy-xy)),gmag=float(r['Gmag']),xy=nxy.tolist()))
        return dict(tic=tic,source=gid,sector=sector,cadence_s=cadence,ffi_receipt=header,
                    field_sha256=hashlib.sha256(field.read_bytes()).hexdigest(), quality0_finite=int(q.sum()),
                    median_flux=float(np.nanmedian(f[q])),gaia_radius_arcsec=130,gaia_rows=len(rows),
                    target_fraction_gaussian_065_085_105_pix=fractions,nearest_gaia=nearest,
                    top_blind_ffi=res['periodic'][0] if res['periodic'] else None,matching_blind_ffi=search_match,
                    blind_ffi_peak_count=len(res['periodic']), localization_ephemeris='SPOC 120-s matching blind peak',
                    pixel_diagnostics=pixels)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('control_results',type=Path)
    p.add_argument('fields',type=Path)
    p.add_argument('cache',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--pixel-tools',type=Path,required=True,help='Repository results/tess-asymmetric-dips/pixels.py')
    a=p.parse_args()
    records=json.loads(a.control_results.read_text())
    px=load_pixel_tools(a.pixel_tools)
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore',message='All-NaN slice encountered')
        results=[run(tic,meta,next(r for r in records if r['tic']==tic and r['sector']==meta[1] and r['cadence_s']==120),
                     a.fields,a.cache,px) for tic,meta in TARGETS.items()]
    a.output.write_text(json.dumps(results,indent=1,allow_nan=False)+'\n')
