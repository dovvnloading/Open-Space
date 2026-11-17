"""Independent W1 forced difference apertures for WISE 0855 on alternating scan visits."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import math
import sys
from pathlib import Path

import numpy as np
from astropy.wcs import WCS
from scipy.ndimage import map_coordinates

def main():
    ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);ap.add_argument('out',type=Path);args=ap.parse_args()
    sys.path.insert(0,str(args.root/'control0855'))
    from native import cutout
    frames=list(csv.DictReader((args.root/'control0855/image-manifest.csv').open()))
    track=next(r for r in csv.DictReader((args.root/'fast/tracks.csv').open()) if r['source_name']=='track1342m076_1')
    a,d=float(track['ra2015_deg']),float(track['dec2015_deg']);cs=math.cos(math.radians(d))
    xx,yy=np.meshgrid(np.arange(-110,111,1.375),np.arange(-80,81,1.375))
    sky_a=a+xx/(3600*cs);sky_d=d+yy/3600
    def get(r):
        u=r['url'].replace('-w2-int-1b','-w1-int-1b');hdu=cutout(u,a,d,half_width=57)
        h,im=hdu.header,hdu.data;wcs=WCS(h).celestial
        x,y=wcs.world_to_pixel_values(sky_a,sky_d); val=map_coordinates(im.astype(float),[y,x],order=1,cval=np.nan)
        # Native DN to Vega nMgy per native pixel.
        val=(val-np.nanmedian(val))*10**(.4*(22.5-float(h['MAGZP'])))
        return r,u,str(h['INHDRSHA']),str(h['INRAWSHA']),val
    with ThreadPoolExecutor(5) as ex: got=list(ex.map(get,frames))
    visits=sorted({int(r['visit']) for r in frames})
    stacks={v:np.nanmedian([im for r,u,hs,rs,im in got if int(r['visit'])==v],axis=0) for v in visits}
    rows=[]; obs=[]
    for r,u,hs,rs,im in got:
        obs.append(dict(visit=r['visit'],mjd=r['mjd'],url=u,input_header_sha256=hs,input_strip_sha256=rs))
    for v in visits:
        t=np.mean([float(r['mjd']) for r in frames if int(r['visit'])==v])
        x0=(t-57170)/365.25*float(track['pmra']);y0=(t-57170)/365.25*float(track['pmdec'])
        other=[stacks[k] for k in visits if k!=v]
        diff=stacks[v]-np.nanmedian(other,axis=0)
        def flux(xc,yc):
            rr=np.hypot(xx-xc,yy-yc); ap=rr<8.25;ann=(rr>=13.75)&(rr<22)
            values=diff[ap]; back=diff[ann]
            if np.isfinite(values).mean()<.95 or np.isfinite(back).mean()<.8:
                return np.nan
            return np.nanmean(values-np.nanmedian(back))*len(values)*.25
        f=flux(x0,y0)
        # Non-overlapping comparison apertures sampled on a staggered grid in this same difference image.
        vals=[flux(x,y) for i,y in enumerate(np.arange(-48,49,17.5)) for x in np.arange(-85,86,17.5)+(i%2)*8.75
              if np.hypot(x-x0,y-y0)>24]
        vals=np.asarray(vals);vals=vals[np.isfinite(vals)]
        print(v,'finite diff',np.isfinite(diff).mean(),'ap',f,'blank',len(vals),flush=True)
        if len(vals)<10:
            continue
        med=np.median(vals);sigma=1.4826*np.median(abs(vals-med))
        rows.append(dict(visit=v,mjd=t,w1_ap_flux_nmgy=f-med,blank_mad_nmgy=sigma,blank_p95_nmgy=np.quantile(vals,.95)-med,n_blank=len(vals)))
    for name,data in [('control_pixel_inputs.csv',obs),('control_pixel_epochs.csv',rows)]:
        with (args.out/name).open('w') as fp:
            wr=csv.DictWriter(fp,fieldnames=list(data[0]));wr.writeheader();wr.writerows(data)
if __name__=='__main__':main()
