"""Ray-traced artificial dark deflectors on the same real VIS backgrounds."""
import argparse
import csv
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi
from search import HERE,MARGIN,RADII,decode_header,digest,find_rings,get_range,sigma
from followup import cut,ring_profile


def source_image(b,offset,angle,amplitude):
    # Supersample 4x per 0.2 arcsec bin, convolve with 0.17 arcsec Gaussian VIS PSF.
    # Surface brightness is conserved along the SIS + external shear ray trace.
    t=(np.arange(65*4)+.5)/4-32.5
    yy,xx=np.meshgrid(t,t,indexing='ij')
    r=np.maximum(np.hypot(xx,yy),1e-8)
    g1=.06*np.cos(2*angle);g2=.06*np.sin(2*angle)
    bx=xx-b*xx/r-g1*xx-g2*yy-offset*b*np.cos(angle+.5)
    by=yy-b*yy/r-g2*xx+g1*yy-offset*b*np.sin(angle+.5)
    u=(bx*np.cos(angle)+by*np.sin(angle))/(.15*b)
    v=(-bx*np.sin(angle)+by*np.cos(angle))/(.11*b)
    src=amplitude*np.exp(-.5*(u*u+v*v))
    conv=ndi.gaussian_filter(src,.17/(.05*2.355))
    return conv.reshape(65,4,65,4).mean(axis=(1,3))


def do_strip(path,field):
    q=json.loads(path.read_text())
    if q['eligible_binned_centers']==0:
        return []
    width=19200;start,end=q['first_byte'],q['last_byte']+1
    step=48*width*4
    with ThreadPoolExecutor(4) as pool:
        bins=list(pool.map(lambda a:get_range(field['url'],a,min(end,a+step)-1), range(start,end,step)))
    raw=b''.join(bins)
    assert digest(raw)==q['sha256']
    a=np.frombuffer(raw,dtype='>f4').reshape(q['height'],width)
    out=[]
    # Different positions across neighbouring chunks; eight injections per realization.
    for ix,x0 in enumerate([3600,9600,14400]):
        arr=a[:,x0:x0+1200].astype('f4')
        mask=np.isfinite(arr)&(arr!=0)&(np.abs(arr)<1e5)
        arr[~mask]=0.
        arr=arr.reshape(q['height']//2,2,600,2).mean(axis=(1,3))
        mask=mask.reshape(q['height']//2,2,600,2).all(axis=(1,3))
        if not mask.any():continue
        _,noise=sigma(arr[mask])
        if not np.isfinite(noise) or noise<=0:continue
        natural,n=find_rings(arr,mask)
        if n==0:continue
        for scenario,offset in enumerate([.02,.18,.45]):
            ai=arr.copy();test=[]
            for j in range(8):
                r=RADII[(j+ix+q['row']//2048)%4]
                amp=(3,8)[j//4]
                cx,cy=75+150*(j%4),48+96*(j//4)
                if not mask[cy-32:cy+33,cx-32:cx+33].all():
                    continue
                phi=.73*j+.43*ix+.29*(q['row']//2048)
                add=source_image(r,offset,phi,amp*noise)
                ai[cy-32:cy+33,cx-32:cx+33]+=add
                test.append((r,amp,cx,cy,add.max()/noise,phi))
            if not test:continue
            peaks,_=find_rings(ai,mask)
            for r,amp,x,y,maxsnr,phi in test:
                matches=[v for v in peaks if np.hypot(v[1]-x,v[2]-y)<=2 and v[3]==r]
                pre=any(np.hypot(v[1]-x,v[2]-y)<=2 and v[3]==r for v in natural)
                recovered=any(ring_profile(cut(ai,v[1],v[2]),v[3])['pass_vis'] for v in matches)
                out.append(dict(tile=q['tile'],row=q['row'],x0_native=x0+2*x+.5,y0_native=q['row']+2*y+.5,
                    theta_arcsec=r*.2,source_offset_over_theta=offset,orientation_rad=round(phi,4),
                    intrinsic_peak_sigma=amp,convolved_peak_sigma=round(maxsnr,4),preexisting=int(pre),
                    stage1=int(bool(matches)),recovered=int(recovered)))
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--tile',required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    a.out.mkdir(exist_ok=True,parents=True)
    field=next(f for f in json.loads((HERE/'fields.json').read_text()) if f['tile']==a.tile)
    results=[]
    for f in sorted((HERE/'runs').glob(a.tile+'_*.json')):
        results.extend(do_strip(f,field))
    if results:
        with (a.out/(a.tile+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
    print(a.tile,len(results),flush=True)
if __name__=='__main__':main()
