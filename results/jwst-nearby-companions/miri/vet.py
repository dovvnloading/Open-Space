"""Forced photometry, stationary-source counts and finite atmospheric colour checks."""
import argparse
import csv
import io
import json
import tarfile
from pathlib import Path
from urllib.request import urlopen

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from fetch_inputs import HERE, digest, download
from search import load, peaks, refine, write_csv, STEP, HALF

BANDS=("F560W","F770W","F1500W","F2100W")
REFS={
    "bobcat.tar.gz":"https://zenodo.org/api/records/5063476/files/evolution_and_photometery.tar.gz/content",
    "jwst_miri_apcorr_0014.fits":"https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_miri_apcorr_0014.fits",
}

def models(path):
    with tarfile.open(path) as tar:
        data=tar.extractfile("photometry_tables/flux_table_JWST+0.0").read().decode().splitlines()
    names=data[6].split()[8:]
    assert names[0]=="F070W"
    a=np.loadtxt(io.StringIO("\n".join(data[8:]).replace("*"," ")))
    chosen=a[:,0]<=1200
    cols=[6+names.index(b) for b in BANDS]
    return a[chosen,:2], 10**a[chosen][:,cols]


def aperture(image,x,y,cal):
    core,si,so=(float(cal[k]) for k in ("radius","skyin","skyout"))
    h=int(np.ceil(so))+2;cx,cy=int(np.rint(x)),int(np.rint(y))
    yy,xx=np.mgrid[cy-h:cy+h+1,cx-h:cx+h+1]
    v=image[cy-h:cy+h+1,cx-h:cx+h+1]
    r=np.hypot(xx-x,yy-y)
    frac=np.clip(core+.5-r,0,1)
    sky=v[(r>=si)&(r<so)]
    if not np.isfinite(v[frac>0]).all() or np.isfinite(sky).mean()<.9:return np.nan
    return float(np.sum((v[frac>0]-np.nanmedian(sky))*frac[frac>0])*cal['apcorr'])


def fluxes(row,src,cal,cache):
    with fits.open(cache/row["filename"]) as hd:
        h=hd['SCI'].header
        # Image units become microJy per native mosaic pixel.
        im=hd['SCI'].data.astype(float)*h["PIXAR_SR"]*1e12
        wt=hd['WHT'].data
        im[wt<=0]=np.nan
        w=WCS(h,naxis=2)
        # Correct the inter-filter translation through the high-S/N host.
        px,py=w.world_to_pixel_values(hd[0].header["TARG_RA"],hd[0].header["TARG_DEC"])
        ix,iy=int(np.rint(px)),int(np.rint(py))
        patch=im[iy-8:iy+9,ix-8:ix+9]; y,x=np.unravel_index(np.nanargmax(patch),patch.shape)
        x,y=refine(ndi.gaussian_filter(patch,.8),x,y)+[ix-8,iy-8]
        from search import sky_grid
        grid=sky_grid(*w.pixel_to_world_values(x,y))
        result=[]
        for s in src:
            east,north=float(s['east_arcsec']),float(s['north_arcsec'])
            def at(e,n):
                u,v=w.world_to_pixel_values(*grid.pixel_to_world_values(HALF+e/STEP,HALF+n/STEP))
                return aperture(im,float(u),float(v),cal)
            flux=at(east,north)
            # Nearby apertures include confusion, diffraction structure, and correlated noise.
            controls=[at(east+radius*np.cos(ang),north+radius*np.sin(ang))
                      for radius in (3.5,6.5) for ang in np.arange(24)*2*np.pi/24]
            controls=np.array(controls)
            sig=1.4826*np.nanmedian(abs(controls-np.nanmedian(controls)))
            result.append((flux,sig))
        return result


def main():
    par=argparse.ArgumentParser(description=__doc__)
    par.add_argument('cache',type=Path);par.add_argument('measured',type=Path)
    args=par.parse_args();out=args.measured
    extra=json.loads((HERE/"aux_inputs.json").read_text())["inputs"]
    standard=json.loads((HERE/"inputs.json").read_text())["inputs"]
    for r in extra:download(r,args.cache)
    refs=json.loads((HERE/'calibration_inputs.json').read_text())
    for r in refs:download(r,args.cache)
    ap=fits.getdata(args.cache/'jwst_miri_apcorr_0014.fits',1)
    model_info,model_flux=models(args.cache/'bobcat.tar.gz')
    with (out/'sources.csv').open() as f:sources=list(csv.DictReader(f))
    reports=[];fields=[]
    meta=json.loads((out/'measurement.json').read_text())
    for target,details in meta.items():
        sel=[r for r in sources if r['target']==target]
        images={r['filters']:r for r in standard+extra if r['target_name']==target}
        got={}
        for band in BANDS:
            cal=ap[(ap['filter']==band)&(ap['subarray']=='FULL')&np.isclose(ap['eefraction'],.7)]
            if len(cal)!=1:raise ValueError(band)
            got[band]=fluxes(images[band],sel,cal[0],args.cache)
        for i,s in enumerate(sel):
            flux,err=np.array([got[b][i] for b in BANDS]).T
            # Free nonnegative scale tests shape only. Keep non-detections as signed fluxes.
            # A 10% term allows for PSF/colour and local calibration systematics.
            var=err**2+(.10*abs(flux))**2
            amp=np.maximum(0,np.sum(model_flux*flux/var,axis=1)/np.sum(model_flux**2/var,axis=1))
            loss=np.sum((flux-model_flux*amp[:,None])**2/var,axis=1)
            idx=np.argmin(loss)
            row=dict(target=target,id=s['id'],astrometry=s['status'])
            for j,b in enumerate(BANDS):
                row[b+'_uJy']=round(flux[j],4);row[b+'_local_MAD_uJy']=round(err[j],4)
            row.update(nbands_above3MAD=int(np.sum(flux>3*err)),best_grid_T_K=model_info[idx,0],
                       best_grid_logg=model_info[idx,1],minimum_weighted_residual=round(loss[idx],3))
            reports.append(row)
        early=load(args.cache/details['epochs'][0]['file'])
        late=load(args.cache/details['epochs'][1]['file'])
        shift=np.array(details['registration']['offset_east_north_arcsec'])/STEP
        ra=early['rad']
        mask=early['good']&(ra>=15)&(ra<29)
        y,x=np.where(mask)
        # Intersect with late coverage at predicted stationary positions.
        xi=np.rint(x+shift[0]).astype(int);yi=np.rint(y+shift[1]).astype(int)
        mask[y,x]&=late['good'][yi,xi]
        p,_=peaks(early['score'],mask)
        q,_=peaks(late['score'],late['good'])
        d,_=cKDTree(q).query(p+shift)
        count=int(np.sum(d*STEP<=.20))
        area=mask.sum()*STEP**2
        fields.append(dict(target=target,control_inner_arcsec=15,control_outer_arcsec=29,
                           common_area_arcsec2=area,n_early_peaks=len(p),stationary_matches=count,
                           stationary_per_arcmin2=3600*count/area,
                           expected_stationary_in_search_annulus=count*np.pi*(12**2-4**2)/area))
    write_csv(out/'colors.csv',reports);write_csv(out/'background_counts.csv',fields)

if __name__=="__main__":main()
