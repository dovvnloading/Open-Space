"""Build sources.csv and long-form registered DR1 aperture spectra for the fixed field."""
import csv
import io
import json
import math
import sys
import warnings
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from acquire import ROOT, REPO, TILE, healpix_nest
sys.path.insert(0,str(REPO/'results/spherex-waste-heat/spectra'))
from read_cubes import unpack
from photometry import extract as aperture_photometry, dist_arcsec

def read(name):
    with (ROOT/'inputs'/f'{name}.csv').open() as f:
        return list(csv.DictReader(f))

def number(x): return float(x) if x not in (None,'','null') else np.nan

def write(path, rows):
    with (ROOT/path).open('w',newline='') as f:
        out=csv.DictWriter(f,fieldnames=rows[0],lineterminator='\n')
        out.writeheader()
        for r in rows:
            out.writerow({k:('' if not np.isfinite(v) else format(v,'.9g')) if isinstance(v,(float,np.floating)) else v for k,v in r.items()})

def associate():
    members,wise,gaia=read('allwise'),read('neighbors'),read('gaia')
    nr,nd=(np.array([number(x[k]) for x in wise]) for k in ('ra','dec'))
    gr,gd,pr,pd=(np.array([number(x[k]) for x in gaia]) for k in ('ra','dec','pmra','pmdec'))
    # Match at approximate AllWISE epoch, then use Gaia at coadd mid-epoch.
    # AllWISE catalog coordinates are profile-fit stationary solutions.
    d_ra=np.nan_to_num(pr)/np.maximum(np.cos(np.deg2rad(gd)),1e-6)/3600000
    d_dec=np.nan_to_num(pd)/3600000
    gr2010,gd2010=gr-5.5*d_ra,gd-5.5*d_dec
    for t in members:
        ra,dec=number(t['ra']),number(t['dec'])
        dist=dist_arcsec(nr,nd,ra,dec)
        dist[np.array([s['designation']==t['designation'] for s in wise])]=np.inf
        ii=int(np.argmin(dist))
        gdists=dist_arcsec(gr2010,gd2010,ra,dec)
        jj=int(np.argmin(gdists)); match=gdists[jj]<=2 and np.sum(gdists<=2)==1
        gaia_teff=number(gaia[jj]['teff_gspphot']) if match else np.nan
        ar=(gr[jj]+9.9*d_ra[jj])%360 if match else ra
        ad=gd[jj]+9.9*d_dec[jj] if match else dec
        tile=healpix_nest(5,ra,dec)
        yield dict(source_id=t['designation'],ra=ra,dec=dec,extraction_ra=ar,extraction_dec=ad,
                   tile=tile,disposition='off_tile' if tile!=TILE else 'unmeasured',pixel_x=np.nan,pixel_y=np.nan,
                   w1mag=number(t['w1mpro']),w2mag=number(t['w2mpro']),cc_flags=t['cc_flags'],
                   ext_flg=t['ext_flg'],nb=t['nb'],nearest_wise_id=wise[ii]['designation'],
                   wise_sep_arcsec=dist[ii],wise_neighbor_w1=number(wise[ii]['w1mpro']),
                   wise_n19=int(np.sum(dist<19)),wise_n60=int(np.sum(dist<60)),
                   # At the 0.5-degree boundary, the 0.51-degree neighbor cap
                   # completely covers 19 arcsec, but can truncate larger radii.
                   gaia_id=gaia[jj]['source_id'] if match else '',gaia_sep_arcsec=gdists[jj] if match else np.nan,
                   gaia_n19=int(np.sum(gdists<19)),teff=gaia_teff,
                   gaia_g=number(gaia[jj]['phot_g_mean_mag']) if match else np.nan,
                   parallax=number(gaia[jj]['parallax']) if match else np.nan,strict_n=0,permissive_n=0)

def main():
    sources = list(associate())
    path = ROOT/'localcache'/f'Npix{TILE}_cube.fits'
    with fits.open(path,do_not_scale_image_data=True) as hd:
        hdr=hd['IMAGE'].header.copy()
        for k in list(hdr):
            if k in ('NAXIS3','CTYPE3','CNAME3','CUNIT3','CRPIX3','CDELT3','CRVAL3','PC3_3') or k.startswith(('PS3_','PV3_')): del hdr[k]
        hdr['NAXIS']=hdr['WCSAXES']=2
        w=WCS(hdr); w.wcs.bounds_check(False,False)
        print('cube',hd[0].header['BLANKCH'],flush=True)
    spec=[]
    byid={s['source_id']:s for s in sources}
    jobs=[]
    for s in sources:
        if s['tile']!=TILE: continue
        x,y=map(float,w.world_to_pixel_values(s['extraction_ra'],s['extraction_dec']))
        if not(-.5<=x<511.5 and -.5<=y<511.5):
            wrap=np.linalg.solve(w.pixel_scale_matrix,[360.,0.])
            valid=[(x+sign*wrap[0],y+sign*wrap[1]) for sign in (-1,1) if -.5<=x+sign*wrap[0]<511.5 and -.5<=y+sign*wrap[1]<511.5]
            if len(valid)!=1:
                s['disposition']='propagated_outside_tile';continue
            x,y=valid[0]
        s['pixel_x'],s['pixel_y']=x,y
        if min(x,y,511-x,511-y)<16:
            s['disposition']='edge_lt16px';continue
        s['disposition']='extracted'
        jobs.append(dict(label=s['source_id'],ra_2025_9=s['extraction_ra'],dec_2025_9=s['extraction_dec']))
    def record(row,stamp):
        s=byid[row['label']]
        with warnings.catch_warnings():
            warnings.simplefilter('ignore',category=RuntimeWarning)
            f,loose,cov,blank,err,nb=aperture_photometry(stamp)
        s['strict_n']=np.isfinite(f).sum();s['permissive_n']=np.isfinite(loose).sum()
        x,y=np.rint(stamp['xy']-stamp['origin']).astype(int)
        for i in range(len(f)):
            spec.append(dict(source_id=s['source_id'],channel=i+1,wave_um=stamp['WAVELENGTH'][i],
                         wave_lo_um=stamp['WL_MIN'][i],wave_hi_um=stamp['WL_MAX'][i],
                         detector=stamp['DETECTOR'][i],subchannel=stamp['SUBCHAN'][i],
                         flux=f[i],sigma=err[i],nhits=stamp['NHITS'][i,y,x],flagged=int(not np.isfinite(f[i])),
                         center_flags=int(stamp['FLAGS'][i,y,x]),coverage=cov[i],n_blank=nb[i],
                         flux_allow_outlier=loose[i],ra=s['ra'],dec=s['dec'],teff=s['teff']))
        # Keep the full blank covariance locally for independent QA, if needed.
        np.savez_compressed(ROOT/'localcache'/(s['source_id']+'_blanks.npz'),blank=blank)
        if len(spec)%510==0: print('sources',len(spec)//102,flush=True)
    write('sources.csv',sources)
    print('apertures',len(jobs),flush=True)
    meta=unpack(path.read_bytes(),jobs,ROOT/'localcache',TILE,sink=record)
    write('sources.csv',sources);write('spectra.csv',spec)
    (ROOT/'receipts'/'extraction.json').write_text(json.dumps(meta,indent=2)+'\n')
    print('complete',len(jobs),len(spec),flush=True)
if __name__=='__main__': main()
