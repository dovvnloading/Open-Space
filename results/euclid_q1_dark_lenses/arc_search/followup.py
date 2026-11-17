"""Recheck every VIS trigger against the sky on *both* sides of the putative ring."""
import argparse
import csv
import json
import math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np
from scipy import ndimage as ndi
from PIL import Image,ImageDraw
from search import HERE, decode_header, digest, find_rings, get_range, sigma

HALF = 24
ANGLES=np.arange(36)*2*np.pi/36


def ring_profile(arr,radius):
    # arr is 49x49, half native resolution; radius here is in binned pixels.
    yy,xx = np.indices(arr.shape)-HALF
    rr=np.hypot(xx,yy)
    sm=ndi.gaussian_filter(arr,1)
    sky = (rr > 1.65*radius+1) & (rr < HALF-1)
    bg,noise=sigma(sm[sky])
    if not np.isfinite(noise) or noise<=0:
        return None
    def polar(rad):
        return ndi.map_coordinates(sm,[HALF+rad*np.sin(ANGLES),HALF+rad*np.cos(ANGLES)],order=1)
    def avg(rads):
        return np.mean([polar(k*radius) for k in rads],axis=0)
    rim=avg([.90,1.,1.10])
    inner=avg([.4,.55])
    outer=avg([1.5,1.65])
    # Bright structure also entering inside or outside an arc fails this comparison.
    contrast=(rim-np.maximum(np.maximum(inner,outer),bg))/noise
    center=(sm[HALF,HALF]-bg)/noise
    signs=contrast>=2.5
    quadrant=signs.reshape(4,9).sum(axis=1)
    # Demand extended structure in at least three quadrants, and cap the share of
    # flux from one patch. A negative/bright core is outside this test.
    pos=np.maximum(rim-bg,0)
    share=max(np.roll(pos,-j)[:6].sum() for j in range(36))/max(pos.sum(),1e-12)
    keep=(signs.sum()>=19 and (quadrant>=4).sum()>=3 and share<.45 and abs(center)<=3.)
    return dict(pass_vis=int(keep),ring_bins=int(signs.sum()),quadrants=int((quadrant>=4).sum()),
                contrast_med=round(float(np.median(contrast)),3),
                peak_share=round(float(share),4),local_center_z=round(float(center),3))


def render(stamp, path):
    finite=stamp[np.isfinite(stamp)]
    bg,s=sigma(finite)
    gray=(255*np.clip((np.arcsinh((stamp-bg)/(2*max(s,1e-12)))+1)/4,0,1))
    Image.fromarray(gray.astype('u1')[::-1]).resize((147,147)).save(path)


def cut(a,x,y):
    return a[y-HALF:y+HALF+1,x-HALF:x+HALF+1]


def check(path,item,out):
    q=json.loads(path.read_text())
    h=get_range(item['url'],0,5759)
    assert digest(h)==q['header_sha256']
    header,offset=decode_header(h)
    start,end=q['first_byte'],q['last_byte']+1
    step=48*int(header['NAXIS1'])*4
    with ThreadPoolExecutor(4) as ex:
        pieces=list(ex.map(lambda a:get_range(item['url'],a,min(a+step,end)-1),range(start,end,step)))
    raw=b''.join(pieces)
    del pieces
    assert digest(raw)==q['sha256']
    width=int(header['NAXIS1'])
    pixels=np.frombuffer(raw,dtype='>f4').reshape(q['height'],width)
    binned=pixels.reshape(q['height']//2,2,width//2,2).mean(axis=(1,3))
    rows=[]
    for j,t in enumerate(q['triggers']):
        x=round((t['x0']-.5)/2)
        y=round((t['y0']-q['row']-.5)/2)
        stamp=cut(binned,x,y)
        measures=ring_profile(stamp,t['radius_arcsec']/.2)
        name=f"{q['tile']}_{q['row']}_{j}"
        if measures is None:
            raise ValueError(f'noise unavailable for {name}')
        rows.append(dict(id=name,**t,**measures))
        if measures['pass_vis']:
            render(stamp,out/(name+'.png'))
    injections=[]
    for j,(r,amp) in enumerate([(r,s) for r in (4,6,9,12) for s in (2,5)]):
        cx,cy = 75+(j%4)*150,48+(j//4)*96
        # Exactly the same injected systems as search.py. Only assess those tested in stage 1.
        hit=next((entry for entry in q['injections'] if abs(entry['r']-.2*r)<1e-6 and entry['peak_sigma']==amp),None)
        if hit is None:
            continue
        sec=np.array(pixels[:,9600:10800],dtype='f4')
        valid=np.isfinite(sec) & (sec!=0) & (np.abs(sec)<1e5)
        good=valid.reshape(q['height']//2,2,600,2).all(axis=(1,3))
        _,s=sigma(binned[:,4800:5400][good])
        stamp=cut(binned,4800+cx,cy)
        yy,xx=np.indices(stamp.shape)-HALF
        add=amp*s*np.exp(-.5*(np.hypot(xx,yy)-r)**2)
        p=ring_profile(stamp+add,r)
        # Re-run at the measured peak with identical masking and a float32 injected chunk.
        sec[~valid]=0
        ai=sec.reshape(q['height']//2,2,600,2).mean(axis=(1,3))
        y1,x1=np.ogrid[-32:33,-32:33]
        ai[cy-32:cy+33,cx-32:cx+33]+=amp*s*np.exp(-.5*(np.hypot(x1,y1)-r)**2)
        again,_=find_rings(ai,good)
        matches=[peak for peak in again if peak[3]==r and math.hypot(peak[1]-cx,peak[2]-cy)<=2]
        if bool(matches)!=hit['recovered']:
            raise ValueError('injection first-stage replay mismatch')
        recovered=False
        for peak in matches:
            res=ring_profile(cut(ai,peak[1],peak[2]),r)
            recovered=recovered or bool(res and res['pass_vis'])
        injections.append(dict(r=.2*r,amplitude=amp,stage1=hit['recovered'],pre=hit['preexisting'],
                               profile_at_injected_center=bool(p and p['pass_vis']),joint_recovered=recovered))

    return rows,injections


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runs',type=Path,default=HERE/'runs')
    p.add_argument('--out',type=Path,default=HERE/'profiles')
    p.add_argument('--stamps',type=Path,default=HERE/'stamps')
    p.add_argument('--tile')
    args=p.parse_args()
    for d in [args.stamps,args.out]:d.mkdir(exist_ok=True,parents=True)
    fields={v['tile']:v for v in json.loads((HERE/'fields.json').read_text())}
    for f in sorted(args.runs.glob('*.json')):
        tile=f.stem.split('_')[0]
        if args.tile and tile!=args.tile:
            continue
        target=args.out/f.name
        if target.exists():
            continue
        rows,inj=check(f,fields[tile],args.stamps)
        tmp=target.with_suffix('.tmp')
        tmp.write_text(json.dumps({'profiles':rows,'injections':inj},indent=2)+'\n')
        tmp.replace(target)
        print(f.stem,len(rows),sum(x['pass_vis'] for x in rows),flush=True)
if __name__=='__main__':
    main()
