"""Compare fixed-sky and fixed host-relative loci in independent survey images."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen
import warnings
import numpy as np
from scipy import ndimage as ndi
from astropy.io import fits
from astropy.wcs import WCS,FITSFixedWarning
from astropy.time import Time

HERE=Path(__file__).resolve().parent
warnings.filterwarnings("ignore",category=FITSFixedWarning)
SKY=(335.782384,-17.612404)
REL=(3.1945,.2542)
HOST=(335.7791533548229,-17.607314733412778)
PM=(.308712,-.718387)
COS=np.cos(np.deg2rad(HOST[1]))

def digest(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def ensure(row,cache):
    p=cache/row['file']
    if not p.exists():
        p.parent.mkdir(exist_ok=True,parents=True)
        with urlopen(row['url'],timeout=120) as s,p.open('wb') as f:
            while b:=s.read(1<<20):f.write(b)
    if digest(p)!=row['sha256']:raise ValueError("Frozen input changed: "+str(p))
    return p

def aperture(im,x,y,r):
    h=int(np.ceil(3.5*r));cx,cy=int(np.rint(x)),int(np.rint(y))
    if min(cx-h,cy-h)<0 or cx+h>=im.shape[1] or cy+h>=im.shape[0]:return np.nan
    dy,dx=np.mgrid[cy-h:cy+h+1,cx-h:cx+h+1]
    rr=np.hypot(dx-x,dy-y)
    w=np.clip(r+.5-rr,0,1)
    z=im[cy-h:cy+h+1,cx-h:cx+h+1];sky=z[(rr>2.3*r)&(rr<3.4*r)]
    use=(w>0)&np.isfinite(z)
    area=w[w>0].sum()
    if w[use].sum()<.75*area or np.isfinite(sky).mean()<.75:return np.nan
    return np.sum((z[use]-np.nanmedian(sky))*w[use])*area/w[use].sum()

def local_host(im,x,y,pix):
    # Within 1.4 arcsec of the catalog prediction. The 95th-percentile-clipped
    # centroid reduces the leverage of a saturated core in VISTA images.
    h=int(np.ceil(1.4/pix));cx,cy=int(np.rint(x)),int(np.rint(y))
    yy,xx=np.mgrid[cy-h:cy+h+1,cx-h:cx+h+1]
    v=im[cy-h:cy+h+1,cx-h:cx+h+1]; floor=np.nanmedian(v)
    wt=np.clip(v-floor,0,np.nanpercentile(v,95)-floor)
    return np.sum(xx*wt)/np.sum(wt),np.sum(yy*wt)/np.sum(wt)

def test_locus(im,w,pos,star,radius,pix):
    x,y=map(float,pos);hx,hy=star
    z=aperture(im,x,y,radius/pix)
    if not np.isfinite(z):return dict(coverage="outside_or_masked")
    dx,dy=x-hx,y-hy;dist=np.hypot(dx,dy);angles=np.linspace(0,2*np.pi,72,endpoint=False)
    # Same-radius control positions sample the stellar halo including diffraction.
    controls=[aperture(im,hx+dist*np.cos(phi),hy+dist*np.sin(phi),radius/pix)
              for phi in angles if np.hypot(dist*np.cos(phi)-dx,dist*np.sin(phi)-dy)>3*radius/pix]
    v=np.array(controls);v=v[np.isfinite(v)];med=np.median(v);scatter=1.4826*np.median(abs(v-med))
    return dict(coverage='measured',x=x,y=y,host_separation_arcsec=dist*pix,
                aperture_minus_azimuth_median=z-med,azimuth_MAD=scatter,
                score=(z-med)/scatter,n_control=len(v),control_percentile=np.mean(v<z))

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('cache',type=Path)
    parser.add_argument('out',type=Path);a=parser.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    rows=[];cats=[]
    peaks=list(csv.DictReader((HERE/'positions.csv').open()))
    for r in json.loads((HERE/'archive_inputs.json').read_text()):
        path=ensure(r,a.cache)
        if r['kind']=='tractor':
            tab=fits.getdata(path,1)
            for s in peaks:
                if s['field']!=r['field']:continue
                ra,dec=float(s['ra']),float(s['dec'])
                d=3600*np.hypot((tab['ra']-ra)*np.cos(np.deg2rad(dec)),tab['dec']-dec)
                i=np.argmin(d)
                cats.append(dict(field=r['field'],peak=s['peak'],catalog='Legacy Surveys DR10 South',separation_arcsec=d[i],
                    id=f"{tab['brickname'][i]}-{tab['objid'][i]}",model=tab['type'][i]))
            continue
        if r['kind']=='ps1':
            t=list(csv.DictReader(path.open()))
            for s in peaks:
                if s['field']!=r['field']:continue
                ra,dec=float(s['ra']),float(s['dec'])
                d=[3600*np.hypot((float(o['raMean'])-ra)*np.cos(np.deg2rad(dec)),float(o['decMean'])-dec) for o in t]
                i=np.argmin(d)
                cats.append(dict(field=r['field'],peak=s['peak'],catalog='Pan-STARRS DR2 mean',separation_arcsec=d[i],
                    id=t[i]['objID'],model=f"nDetections={t[i]['nDetections']}"))
            continue
        with fits.open(path) as hd:
            h=hd[1];w=WCS(h.header,naxis=2);im=h.data.astype(float)
            hdr=dict(hd[0].header);hdr.update(dict(h.header))
            dt=hdr['DATE-OBS'];yr=Time(dt).jyear-2000.
            expected=(HOST[0]+PM[0]*yr/(3600*COS),HOST[1]+PM[1]*yr/3600)
            pix=np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600
            sx,sy=local_host(im,*w.world_to_pixel_values(*expected),pix)
            hRA,hDE=w.pixel_to_world_values(sx,sy)
            moving=(hRA+REL[0]/(3600*COS),hDE+REL[1]/3600)
            if r['kind']=='HST':
                # Filter standard DQ failures while preserving the stellar center determined above.
                im[hd['DQ'].data!=0]=np.nan
            for tag,coord in [('stationary',SKY),('fixed_host_offset',moving)]:
                found=test_locus(im,w,w.world_to_pixel_values(*coord),(sx,sy),.25 if r['kind']=='HST' else .85,pix)
                rows.append(dict(dataset=r['file'],date=dt,filter=hdr['FILTER'],hypothesis=tag,ra=coord[0],dec=coord[1],
                    units=hdr.get('BUNIT',''),host_x=sx,host_y=sy,**found))
    for filename,records in [('historical_loci.csv',rows),('historical_catalog_matches.csv',cats)]:
        keys=list(dict.fromkeys(k for row in records for k in row))
        with (a.out/filename).open('w',newline='') as f:
            w=csv.DictWriter(f,keys);w.writeheader();w.writerows(records)
if __name__=="__main__":main()
