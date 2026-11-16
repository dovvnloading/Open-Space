"""Exploratory motion-aware SPHEREx cutout photometry; not discovery inference."""
import argparse, concurrent.futures, hashlib, json, math, re, time, urllib.request, subprocess
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

def parse_most(path):
    lines=Path(path).read_text().splitlines()
    header=next(x for x in lines if x.startswith('|'))
    bounds=[m.start() for m in re.finditer(r'\|',header)]
    names=[header[a+1:b].strip() for a,b in zip(bounds,bounds[1:])]
    return [{n:line[a+1:b].strip() for n,a,b in zip(names,bounds,bounds[1:])} for line in lines if line.strip() and not line.startswith(('\\','|'))]

def aperture(im,var,flags,x,y,r,allow_outlier=False):
    if x-r<0 or y-r<0 or x+r>=im.shape[1] or y+r>=im.shape[0]: return None
    yy,xx=np.indices(im.shape);rr=np.hypot(xx-x,yy-y)
    good=np.isfinite(im)&np.isfinite(var)&(var>0)&((flags & ~(2097152 | (524288 if allow_outlier else 0)))==0)
    ap=rr<=r; sky=(rr>=8)&(rr<=12)&good
    if np.sum(sky)<50 or not np.all(good[ap]) or not np.any(ap):return None
    b=float(np.median(im[sky])); n=int(ap.sum())
    # Median uncertainty approximate normal-efficiency factor; no claim of full calibration.
    sigma_b=1.4826*np.median(np.abs(im[sky]-b))/np.sqrt(sky.sum())*1.2533
    flux=float(np.sum(im[ap]-b)); error=float(np.sqrt(np.sum(var[ap])+n*n*sigma_b*sigma_b))
    return dict(sum_mjysr=flux,error_mjysr=error,snr=flux/error,npix=n,background=b,sky_n=int(sky.sum()))

def process(row,cache):
    url='https://irsa.ipac.caltech.edu/'+row['uri']+'?center='+row['ra_obj']+','+row['dec_obj']+'&size=0.1'
    # Include target center and cutout size: an exposure can contain many targets.
    key=hashlib.sha256(url.encode()).hexdigest(); path=cache/(key+'.fits')
    out=dict(observationid=row['observationid'],detector=row['productid'],start=row['start_date'],ra=float(row['ra_obj']),dec=float(row['dec_obj']),url=url)
    try:
        if not path.exists():
            data=subprocess.check_output(['curl','--fail','--silent','--show-error','--location','--max-time','90','--retry','1',url])
            if not data.startswith(b'SIMPLE'): raise ValueError('not FITS: '+data[:120].decode(errors='replace'))
            path.write_bytes(data)
        out['sha256']=hashlib.sha256(path.read_bytes()).hexdigest()
        with fits.open(path) as hdus:
            h=hdus['IMAGE'].header; im=hdus['IMAGE'].data;var=hdus['VARIANCE'].data;flags=hdus['FLAGS'].data
            x,y=WCS(h,hdus).all_world2pix([[out['ra'],out['dec']]],0)[0]
            # Spectral TAB coordinates operate on undistorted detector pixels.
            hdr=h.copy()
            for k in list(hdr):
                if re.match(r'^(A|B|AP|BP)_(ORDER|\d+_\d+)$',k):del hdr[k]
            wave,width=WCS(hdr,hdus,key='W').wcs_pix2world([[x,y]],0)[0]
            out.update(x=float(x),y=float(y),wavelength_um=float(wave),bandwidth_um=float(width),mjd=float(h['MJD-AVG']),image_quality=h.get('L2DQAFLG'),pipeline=hdus[0].header.get('VERSION'))
            out['apertures']={str(r):aperture(im,var,flags,x,y,r) for r in (2.5,3.5,5.)}
            out['diagnostic_outlier_inclusive']={str(r):aperture(im,var,flags,x,y,r,True) for r in (2.5,3.5,5.)}
            yy,xx=np.indices(im.shape);core=np.hypot(xx-x,yy-y)<=2.5
            out['core_flag_counts']={str(int(v)):int(np.sum(flags[core]==v)) for v in np.unique(flags[core])}
            out['controls']=[aperture(im,var,flags,x+dx,y+dy,2.5) for dx,dy in [(15,0),(-15,0),(0,15),(0,-15)]]
    except Exception as e:out['error']=str(e)
    return out

def main():
    p=argparse.ArgumentParser();p.add_argument('manifest');p.add_argument('--cache',default='/tmp/spherex-cache');p.add_argument('--output',required=True);p.add_argument('--detectors',nargs='+',default=['D4','D5','D6']);p.add_argument('--limit',type=int);a=p.parse_args()
    rows=[r for r in parse_most(a.manifest) if r['productid'] in a.detectors]
    if a.limit is not None: rows=rows[:a.limit]
    cache=Path(a.cache);cache.mkdir(parents=True,exist_ok=True)
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as ex:
        for r in ex.map(lambda r:process(r,cache),rows):
            results.append(r);Path(a.output).write_text(json.dumps(results,indent=2,allow_nan=False)+'\n');print(len(results),r['observationid'],r.get('wavelength_um'),r.get('error'),flush=True)
if __name__=='__main__':main()
