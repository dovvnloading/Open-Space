"""Audit MIRI detector footprints and measure Neptune's known satellite Triton.

Requires the frozen observation/product CSVs from archive.py. Output visibility
is evaluated on the optical projection of the rectangular detector; data
quality and surface-brightness tests apply to the eight direct Triton images.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import io
from pathlib import Path
import time
from urllib.request import Request, urlopen
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib.path import Path as Polygon
from scipy.ndimage import gaussian_filter, map_coordinates
from geometry import ARCSEC, AU, PLANET, basis, direction, query, tangent_basis

RADII=[50000,70000,100000,150000,250000,400000,600000]


def read_csv(path):
    with path.open() as f:return list(csv.DictReader(f))


def write_csv(path, fields, rows):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fields,lineterminator='\n');w.writeheader();w.writerows(rows)


def source_range(row, cache, offset, size, receipts):
    n=row['productFilename']
    p=cache/f"{n}.{offset}.{size}.bin"
    if not p.exists():
        url=f"https://stpubdata.s3.amazonaws.com/jwst/public/jw{n[2:7]}/{n.split('_')[0]}/{n}"
        end=offset+size-1
        for retry in range(4):
            try:
                req=Request(url+f'?bytes={offset}-{end}',headers={'Range':f'bytes={offset}-{end}'})
                with urlopen(req,timeout=90) as r:
                    if r.status!=206 or not r.headers.get('Content-Range','').startswith(f'bytes {offset}-{end}/'):
                        raise OSError("server refused byte range")
                    data=r.read(size+1)
                if len(data)!=size:raise OSError("short FITS byte range")
                tmp=p.with_suffix('.part');tmp.write_bytes(data);tmp.replace(p);break
            except OSError:
                if retry==3:raise
                time.sleep(2**retry)
    b=p.read_bytes()
    if len(b)!=size:raise OSError("invalid cached byte range")
    receipts[(n,offset,size)]=hashlib.sha256(b).hexdigest()
    return b


def headers(row,cache,receipts):
    f=io.BytesIO(source_range(row,cache,0,65536,receipts))
    primary=fits.Header.fromfile(f);sci=fits.Header.fromfile(f)
    return row,primary,sci,f.tell()


def perimeter(sci,wcs,planet):
    nx,ny=sci['NAXIS1'],sci['NAXIS2']
    t=np.linspace(0,1,129)
    xx=np.r_[t*(nx-1),np.full_like(t,nx-1),t[::-1]*(nx-1),t*0]
    yy=np.r_[t*0,t*(ny-1),np.full_like(t,ny-1),t[::-1]*(ny-1)]
    ra,dec=wcs.all_pix2world(xx,yy,0);v=direction(ra,dec)
    l,e,n=tangent_basis(*planet[:2])
    pts=ARCSEC*np.column_stack((v@e,v@n))/(v@l)[:,None]
    poly=Polygon(pts)
    # Distance from the planet origin to the piecewise-linear perimeter.
    a,b=pts[:-1],pts[1:];u=b-a
    f=np.clip(-np.sum(a*u,axis=1)/np.maximum((u*u).sum(axis=1),1.e-30),0.,1.)
    gap=np.min(np.linalg.norm(a+f[:,None]*u,axis=1))
    return poly, (0. if poly.contains_point((0.,0.)) else gap)


def plane_coverage(planet,mjd,coords,poly):
    u,v,h=basis(planet,mjd+2400000.5)
    l,e,n=tangent_basis(*coords[:2]); dist=coords[4]*AU
    phase=np.arange(720)*2*np.pi/720
    q=np.cos(phase)[:,None]*u+np.sin(phase)[:,None]*v
    result={}
    for radius in RADII:
        xy=ARCSEC*radius*np.column_stack((q@e,q@n))/(dist+radius*(q@l))[:,None]
        result['phase_fraction_'+str(radius)]=round(float(np.mean(poly.contains_points(xy))),5)
    return result


def triton_pixels(item, cache, receipts, xy, planet_xy):
    row,ph,h,off=item
    nx,ny=h['NAXIS1'],h['NAXIS2'];x,y=xy
    r=26
    x0,x1=int(x)-r,int(x)+r+1;y0,y1=int(y)-r,int(y)+r+1
    if min(x0,y0)<0 or x1>nx or y1>ny:raise RuntimeError("Triton stamp truncated by detector boundary")
    def plane(head,origin,dtype):
        b=source_range(row,cache,origin+y0*nx*4,(y1-y0)*nx*4,receipts)
        return np.frombuffer(b,dtype).reshape(y1-y0,nx)[:,x0:x1]
    img=plane(h,off,'>f4')
    off+=((h.data_size+2879)//2880)*2880
    for kind in ('ERR','DQ'):
        head=fits.Header.fromfile(io.BytesIO(source_range(row,cache,off,2880,receipts)))
        assert head['EXTNAME']==kind
        off+=2880
        if kind=='DQ':
            dq=plane(head,off,'>i4').astype('uint32')^np.uint32(head.get('BZERO',0))
        else:
            err=plane(head,off,'>f4')
        off+=((head.data_size+2879)//2880)*2880
    good=(dq==0)&np.isfinite(img)&np.isfinite(err)&(err>0)
    yy,xx=np.indices(img.shape)
    rr=np.hypot(xx+x0-x,yy+y0-y)
    psep=np.hypot(xx+x0-planet_xy[0],yy+y0-planet_xy[1])
    # Mask Neptune's immediate vicinity for the local fit. Subtract a quadratic
    # surface to allow curvature in the planet and instrumental background.
    ring=(rr>=9)&(rr<=20)&(psep>12)&good
    X=(xx+x0-x)/12;Y=(yy+y0-y)/12
    design=np.stack([np.ones_like(X),X,Y,X*X,X*Y,Y*Y],axis=-1)
    sel=ring.copy()
    for _ in range(3):
        coeff=np.linalg.lstsq(design[sel],img[sel],rcond=None)[0]
        res=img-design@coeff
        noise=1.4826*np.nanmedian(abs(res[sel]-np.nanmedian(res[sel])))
        sel=ring&(np.abs(res)<max(3*noise,1.e-6))
    aperture=(rr<=4)&good
    centers=np.arange(-12,13)
    kernel_sigma=2.8 if ph['FILTER']=='F2100W' else 3.3
    norm=gaussian_filter(good.astype(float),kernel_sigma)
    sm=gaussian_filter(np.where(good,res,0),kernel_sigma)/np.maximum(norm,.01)
    cn=(rr<=4)&(norm>.8)
    peak=np.nanargmax(np.where(cn,sm,np.nan));dy,dx=np.unravel_index(peak,sm.shape)
    local=(rr>=9)&(rr<=17)&(norm>.9)&(psep>12)
    med=np.nanmedian(sm[local]);sigma=1.4826*np.nanmedian(abs(sm[local]-med))
    if not np.isfinite(sigma) or sigma<=0:raise RuntimeError("local control scatter undefined")
    flux=np.sum(res[aperture])*1e6*h['PIXAR_SR']*1000
    return dict(file=row['productFilename'],visit=ph['VISIT_ID'],filter=ph['FILTER'],
                mjd_utc=ph['EXPMID'],triton_x=float(x),triton_y=float(y),planet_sep_px=np.hypot(x-planet_xy[0],y-planet_xy[1]),
                good_aperture_fraction=aperture.sum()/np.sum(rr<=4),aperture_mJy=flux,
                peak_dx=float(dx+x0-x),peak_dy=float(dy+y0-y),peak_local_MAD=(sm[dy,dx]-med)/sigma,
                expected_local_MAD=(float(map_coordinates(sm,[[y-y0],[x-x0]],order=1)[0])-med)/sigma)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("manifests",type=Path);p.add_argument("cache",type=Path)
    p.add_argument("horizons",type=Path);p.add_argument("out",type=Path)
    a=p.parse_args();a.cache.mkdir(parents=True,exist_ok=True);a.out.mkdir(parents=True,exist_ok=True)
    obs={r['obsid']:r for r in read_csv(a.manifests/'observations.csv')}
    rows=[r for r in read_csv(a.manifests/'products.csv') if r['productSubGroupDescription']=='CAL'
          and any(obs[i]['instrument_name']=='MIRI/IMAGE' for i in r['parent_obsids'].split('|'))]
    receipts={}
    with ThreadPoolExecutor(12) as pool:
        info=list(pool.map(lambda r:headers(r,a.cache,receipts),rows))
    footprints=[];targets=[]
    for pl in ('Uranus','Neptune'):
        group=[item for item in info if item[0]['planet']==pl]
        ep=query(PLANET[pl],'JWST',np.array([g[1]['EXPMID'] for g in group])+2400000.5,a.horizons)
        for item,state in zip(group,ep):
            row,ph,sci,_=item;w=WCS(sci,naxis=2)
            poly,gap=perimeter(sci,w,state);frac=plane_coverage(pl,ph['EXPMID'],state,poly)
            # Only interpret inverse pixels if the point is inside the detector.
            xx,yy=(w.all_world2pix(*state[:2],0) if gap==0 else ('',''))
            footprints.append(dict(planet=pl,file=row['productFilename'],proposal=ph['PROGRAM'],visit=ph['VISIT_ID'],
                                   target=ph['TARGNAME'],planned_aperture=ph['PPS_APER'],exp_type=ph['EXP_TYPE'],
                                   filter=ph['FILTER'],exp_mid_mjd_utc=ph['EXPMID'],exp_secs=ph['EFFEXPTM'],
                                   planet_on_chip=int(gap==0),planet_x=xx,planet_y=yy,planet_to_chip_arcsec=round(float(gap),3),**frac))
            if ph['PPS_APER']=='MIRIM_ILLUM' and ph['TARGNAME'].upper()=='TRITON' and gap==0:
                targets.append((item,(xx,yy)))
    write_csv(a.out/'footprints.csv',list(footprints[0]),footprints)
    eph=query(801,'JWST',np.array([item[1]['EXPMID'] for item,s in targets])+2400000.5,a.horizons)
    tasks=[]
    for (item,s),e in zip(targets,eph):
        xy=WCS(item[2],naxis=2).all_world2pix(*e[:2],0)
        tasks.append((item,a.cache,receipts,xy,s))
    with ThreadPoolExecutor(6) as pool:controls=list(pool.map(lambda t:triton_pixels(*t),tasks))
    write_csv(a.out/'triton_control.csv',list(controls[0]),controls)
    write_csv(a.out/'byte_ranges.csv',['file','offset','length','sha256'],
              [dict(file=k[0],offset=k[1],length=k[2],sha256=v) for k,v in sorted(receipts.items())])


if __name__=='__main__':main()
