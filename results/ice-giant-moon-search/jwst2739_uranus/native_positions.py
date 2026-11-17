"""Measure native trailed PSF centers at the 2023 known-moon longitude maximum."""
import argparse
import csv
import importlib.util
import json
from pathlib import Path
import sys
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from screen import geometry,search,root
from prepare import inverse_sky

# Import the January fit and rotation without colliding with this directory's longitude driver.
sys.path.insert(0,str(root/'jwst2768'))
from centroids import fit_blob
from perdita import rotate_tracks


def run(path, files, sweep, out, body, fixed_offset=None):
    c=search.load_npz(path);vis=np.unique(c['visit'])
    with sweep.open() as f: table=list(csv.DictReader(f))
    angles=[]
    for name in ['feb','sep']:
        near=[r for r in table if abs(float(r['offset_deg']))<=5 and int(r[name+'_n'])>=2]
        angles.append(float(max(near,key=lambda r:float(r[name+'_snr']))['offset_deg']))
    if fixed_offset is not None:
        angles=[fixed_offset]*len(vis)
    x,y=rotate_tracks(c,np.array(angles),body)
    rows=[];mid=x.shape[-1]//2
    planet=c['center_xy'];scale=float(c['scale_arcsec'])/geometry.ARCSEC
    for i,r in enumerate(json.loads(Path(str(path)+'.inputs.json').read_text())):
        j=list(vis).index(c['visit'][i])
        qx=x[j,i];qy=y[j,i]
        with fits.open(files/r['filename']) as h:
            w=WCS(h['SCI'].header,naxis=2)
            cr,cd=w.all_pix2world(r['center_native_x'],r['center_native_y'],0)
            wr,wd=inverse_sky(cr,cd,-(qx-planet[0])*scale,(qy-planet[1])*scale)
            xx,yy=w.all_world2pix(wr,wd,0)
            res=fit_blob(xx,yy,h['SCI'].data,h['ERR'].data,h['DQ'].data,float(c['psf_sigma_px']))
            if res is None:continue
            dx,dy,snr,err,scatter,n=res
            wr,wd=w.all_pix2world(dx,dy,0)
            sx,sy=geometry.sky_xy(wr,wd,cr,cd,float(c['scale_arcsec']),planet)
            ra,dec=inverse_sky(c['path_ra'][i,mid],c['path_dec'][i,mid],-(sx-planet[0])*scale,(sy-planet[1])*scale)
            p0=c['known_xy'][i,list(c['known_id']).index(body)]
            pos=c['known_xy'][i];d=np.hypot(pos[:,0]-sx,pos[:,1]-sy);d[list(c['known_id']).index(body)]=np.inf
            # Estimate local longitude by the sky position of the nearest rotated track.
            grid=np.arange(angles[j]-4,angles[j]+4,.01)
            ax,ay=rotate_tracks(c,grid,body)
            ll=grid[np.argmin(np.hypot(ax[:,i,mid]-sx,ay[:,i,mid]-sy))]
            rows.append(dict(frame=r['exposure_id'],visit=c['visit'][i],body=body,mjd_utc=c['mjd'][i],
              ra_deg=float(ra),dec_deg=float(dec),sky_x=float(sx),sky_y=float(sy),peak_longitude_deg=float(ll),
              native_x=float(dx),native_y=float(dy),sigma_x_fit_px=float(err[0]),sigma_y_fit_px=float(err[1]),
              snr=float(snr),nominal_separation_px=float(np.linalg.norm([sx-p0[0],sy-p0[1]])),
              nearest_other_id=int(c['known_id'][np.argmin(d)]),nearest_other_separation_px=float(min(d))))
    columns=['frame','visit','body','mjd_utc','ra_deg','dec_deg','sky_x','sky_y','peak_longitude_deg','native_x','native_y','sigma_x_fit_px','sigma_y_fit_px','snr','nominal_separation_px','nearest_other_id','nearest_other_separation_px']
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(rows)
    for v in vis:
        group=[r for r in rows if r['visit']==v]
        print(v,len(group),[round(float(np.median([r[k] for r in group])),3) if group else None for k in ['peak_longitude_deg','snr','nominal_separation_px']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','files','sweep','out']:p.add_argument(k,type=Path)
    p.add_argument('--body',type=int,default=725)
    p.add_argument('--offset',type=float,help='use the same pre-specified longitude offset in every visit')
    a=p.parse_args();run(a.cube,a.files,a.sweep,a.out,a.body,a.offset)
