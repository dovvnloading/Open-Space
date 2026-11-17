"""Diagnostic ePSF wing fits for saturated stellar candidates and one quiet-core control.

Native IMAGE-ZODI, only source-only or clear wing pixels; fit two predeclared
wing annuli separately at the catalog position. No fit to the saturated core.
This control does not deblend or provide production-quality saturated-star fluxes.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates
from audit_area import local_area
from measure import analyse, SOURCE

TARGETS = {
    'nep': ['J180309.59+662811.5','J180348.24+664107.4'],
    'wide': ['J135641.67+301025.8','J135722.25+294153.7','J135745.42+295804.4']}


def wing(hd, ra, dec, inner):
    h=hd['IMAGE'].header
    w=WCS(h)
    x,y=w.all_world2pix([[ra,dec]],0)[0]
    ax=WCS(h,key='A')
    ax.sip=None
    a,b=ax.all_pix2world([[x,y]],0)[0]
    psftab=hd['EPSF'].data
    k=np.argmin((psftab['XCENTER']-a)**2+(psftab['YCENTER']-b)**2)
    e=psftab['EPSF'][k]
    oy,ox=hd['EPSF'].header['OVSMPY'],hd['EPSF'].header['OVSMPX']
    signal=hd['IMAGE'].data-hd['ZODI'].data
    flags=hd['FLAGS'].data
    var=hd['VARIANCE'].data
    dy,dx=np.indices(signal.shape,dtype=float)
    dx-=x;dy-=y
    radius=np.hypot(dx,dy)
    valid=np.isfinite(signal+var) & (var>0) & ((flags & ~SOURCE)==0)
    # The mission SOURCE mask extends throughout the local sky of these bright
    # targets. Use its clean subset if present, otherwise the median of clean
    # or SOURCE-only pixels in the larger annulus, as for aperture photometry.
    sky=valid & (radius>=6.) & (radius<=11.)
    clear_sky=sky & (flags==0)
    if clear_sky.sum()>=20:
        sky=clear_sky
    use=valid & (radius>=inner) & (radius<=3.)
    quadrants=set(zip(dx[use]>0,dy[use]>0))
    if use.sum()<10 or len(quadrants)<3 or sky.sum()<20:
        return {'status':'insufficient_wing'}
    psf=map_coordinates(e,[dy*oy+(e.shape[0]-1)/2,dx*ox+(e.shape[1]-1)/2],order=1,mode='constant')
    if not np.isfinite(psf).all() or psf.sum() <= 0:
        return {'status':'invalid_epsf'}
    psf/=psf.sum()
    p=psf[use]
    yv=signal[use]-np.median(signal[sky])
    weight=1/var[use]
    den=np.sum(weight*p**2)
    amp=np.sum(weight*p*yv)/den
    chi=np.sum(weight*(yv-amp*p)**2)/(len(p)-1)
    sigma=np.sqrt(max(1,chi)/den)
    scale=1.e9*local_area(w,x,y)
    return dict(status='ok',flux=float(amp*scale),sigma=float(sigma*scale),n_wing=int(use.sum()),
                kept_fraction=float(p.sum()),epsf_bin=int(k),reduced_chi2=float(chi))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('visits',type=Path)
    parser.add_argument('cache',type=Path)
    args=parser.parse_args()
    output=[]; summaries=[]
    for key in TARGETS[args.visits.name]:
        summary=json.loads((args.visits/(key+'.json')).read_text())
        for inner in (1.25,1.75):
            rr=[]
            for row in csv.DictReader((args.visits/(key+'.csv')).open()):
                if row['status'] in ('off_bin','edge') or 'read_error' in row['status']:
                    continue
                p=args.cache/(hashlib.sha256(row['url'].encode()).hexdigest()+'.fits')
                with fits.open(p) as hd:
                    test=wing(hd,summary['ra'],summary['dec'],inner)
                r={name:row[name] for name in ('release','bin','mjd','wave_um','date_obs','sha256')}
                r['mjd']=float(r['mjd']);r['wave_um']=float(r['wave_um'])
                r.update(test);rr.append(r)
                output.append(dict(source_id=key,inner=inner,**r))
            groups=analyse(rr)
            for release in ('dr1_all','qr3_all'):
                g=groups[release]
                summaries.append(dict(source_id=key,inner=inner,release=release,blue=g['bands'].get('blue',{}).get('n',0),
                                      core=g['bands'].get('core',{}).get('n',0),red=g['bands'].get('red',{}).get('n',0),
                                      fractional_contrast=g.get('fractional_contrast'),signed_score=g.get('signed_score')))
    for name,rows in [('wing_fits',output),('wing_summary',summaries)]:
        with (args.visits/(name+'.csv')).open('w') as f:
            writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
            writer.writeheader();writer.writerows(rows)


if __name__=='__main__':
    main()
