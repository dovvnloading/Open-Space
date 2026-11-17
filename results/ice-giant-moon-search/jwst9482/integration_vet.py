"""Resolve the 9482 CAL trail into 60 independent readout integrations."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import warnings
import numpy as np
from link_visits import geometry,search,trials
import positive_control as pc
from prepare import inverse_sky
from independent_vet import save
ROOT=Path(__file__).resolve().parent


def run(args):
    args.cache.mkdir(parents=True,exist_ok=True);args.out.mkdir(parents=True,exist_ok=True)
    c=search.load_npz(args.cube)
    native=json.loads((ROOT/'inputs.json').read_text())
    roots=[str(s).removesuffix('_cal') for s in c['exposure_id']]
    with ThreadPoolExecutor(6) as pool:meta=list(pool.map(lambda s:pc.metadata(args.cache,s),roots))
    config=json.loads((ROOT/'independent_vet.json').read_text())['F150W2']
    a,phi,epoch=(config[k] for k in ['radius_km','phase_deg','emission_epoch_tdb_jd'])
    with TemporaryDirectory() as tmp:
        d=Path(tmp);el=np.array([[a,0,0,0,0,phi]])
        orbit=trials(d/'t.csv',args.cube,d/'o.npz',el,epoch)
    patches={k:[] for k in ['trial','Perdita','Belinda','S2025U1']}
    output=[]
    for j,(m,record) in enumerate(zip(meta,native,strict=True)):
        times=m['times'];orig=c['path_mjd'][j]
        if record['exposure_id']!=roots[j]+'_cal':raise ValueError('frame ordering')
        paths=dict(trial=np.stack([orbit['x_path'][0,j],orbit['y_path'][0,j]],axis=-1))
        for k,b in [('Perdita',725),('Belinda',714),('S2025U1',75052)]:
            paths[k]=c['known_path'][np.flatnonzero(c['known_id']==b)[0],j]
        # These are detector pixels relative to the measured Uranus limb center.
        ref=m['wcs'].all_pix2world(record['center_native_x'],record['center_native_y'],0)
        tracks={}
        for name,track in paths.items():
            x,y=[np.interp(times,orig,coord) for coord in track.T]
            east=-(x-c['center_xy'][0])*float(c['scale_arcsec'])/geometry.ARCSEC
            north=(y-c['center_xy'][1])*float(c['scale_arcsec'])/geometry.ARCSEC
            tracks[name]=np.array(m['wcs'].all_world2pix(*inverse_sky(*ref,east,north),0,quiet=True)).T
        residuals,x0,y0=pc.wide_images(args.cache,m)
        images=pc.patches_for_positions(residuals,x0,y0,tracks)
        for name,p in images.items():
            patches[name].extend(p)
            for i,(px,py) in enumerate(tracks[name]):
                output.append(dict(body=name,file=m['file'],integration=i+1,mjd_utc=times[i],x=round(px,4),y=round(py,4),
                                  filtered_DN_s=p[i][pc.RADIUS,pc.RADIUS]))
        print(m['file'],len(times),flush=True)
    save(args.out/'integration_positions.csv',output)
    fig_summary=[]
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(4,3,figsize=(8.7,9.3),layout='constrained')
    for j,(name,arr) in enumerate(patches.items()):
        stack=np.array(arr); n=len(stack)//2
        combined=np.nanmedian(stack,axis=0)
        dx,dy=pc.summary(combined)[2:4]
        for k,(label,subset) in enumerate([('First 30',stack[:n]),('Last 30',stack[n:]),('All 60',stack)]):
            im=np.nanmedian(subset,axis=0);score,peak,ox,oy,mad=pc.summary(im)
            yy,xx=np.indices(im.shape)-pc.RADIUS
            bg=np.nanmedian(im[(xx**2+yy**2>=64)&(xx**2+yy**2<=400)])
            z=(im-bg)/mad
            fig_summary.append(dict(body=name,subset=label,central_score=round(float(score),3),peak_3px=round(float(peak),3),
                               peak_dx=int(ox),peak_dy=int(oy),score_at_full_peak=round(float(z[pc.RADIUS+dy,pc.RADIUS+dx]),3)))
            ax=axes[j,k];ax.imshow(z,origin='lower',vmin=-3,vmax=7,cmap='gray',extent=(-24.5,24.5,-24.5,24.5));ax.plot(0,0,'+',color='#f08445',ms=9)
            if j==0:ax.set_title(label+' integrations')
            if k==0:ax.set_ylabel(name+'   Δy (px)')
            if j==3:ax.set_xlabel('Δx (px)')
            ax.set_xlim(-14,14);ax.set_ylim(-14,14)
    save(args.out/'integration_summary.csv',fig_summary)
    fig.savefig(args.out/'integration_stacks.png',dpi=125)
    receipts=[]
    for p in sorted(pc.USED):
        name,start,length,_=p.name.rsplit('.',3)
        receipts.append(dict(file=name,offset=int(start),length=int(length),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    save(args.out/'integration_ranges.csv',receipts)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ['cube','cache','out']:ap.add_argument(k,type=Path)
    run(ap.parse_args())
