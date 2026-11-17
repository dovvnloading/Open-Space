"""Adjudicate annual aliases with the second 1272 visit, per-filter persistence and native morphology."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from matplotlib import pyplot as plt
from coverage import OUT,frame,get_states,pixels
from fetch import products,read_product
from native_search import SIGMA,maps,visit_scores
from depth import get_xy,sample
from miri_audit import write_csv
from vet import photometry


def run(cache,head):
    rr=list(csv.DictReader((OUT/'neptune_all_leads.csv').open()));src=json.loads((OUT/'neptune_all_summary.json').read_text())
    radii=np.array([float(r['a_km']) for r in rr]);ph=np.deg2rad([float(r['phase_deg']) for r in rr]);states=get_states()
    fs=[];native=[]
    for f,p in products():
        if f['planet']!='Neptune' or f['planet_on_chip']!='1':continue
        g=frame(f,p,head,states[f['file']],{});g['p']=p;fs.append(g)
    fs.sort(key=lambda g:g['t'])
    for g in fs:native.append(maps(g,g['p'],cache))
    xy=get_xy(fs,radii,ph,src['epoch_jd_tdb']);z,ok=sample(fs,native,xy)
    # Require a full additional diffraction FWHM of support to test annular glare at the mask boundary.
    eroded=[(sc,ndi.distance_transform_edt(good)>2.355*SIGMA[g['f']['filter']]) for (sc,good),g in zip(native,fs)]
    ze,oke=sample(fs,eroded,xy);es,eligible,_,_,_=visit_scores(ze,oke,fs,'Neptune')
    scored,_,v,vc,vs=visit_scores(z,ok,fs,'Neptune');assert np.allclose(scored,[float(r['rank']) for r in rr])
    # Independent negative evidence is informative only when the model is actually covered.
    obs_visit=np.array([g['f']['visit'] for g in fs]);bands=np.array([g['f']['filter'] for g in fs]);allrows=[]
    for i,r in enumerate(rr):
        q=dict(a_km=radii[i],phase_deg=np.rad2deg(ph[i]),rank=scored[i],PSF_clearance_eligible=int(eligible[i]),PSF_clearance_rank=es[i])
        for j,visit in enumerate(vs):
            q['visit_'+visit+'_covered']=int(vc[i,j]);q['visit_'+visit+'_rank']=v[i,j]
            for b in sorted(set(bands[obs_visit==visit])):
                take=ok[i]&(obs_visit==visit)&(bands==b)
                q[visit+'_'+b+'_good']=sum(take)
                q[visit+'_'+b+'_mean']=float(np.mean(np.clip(z[i,take],-7,7))) if take.any() else np.nan
        allrows.append(q)
    write_csv(OUT/'neptune_alias_validation.csv',list(allrows[0]),allrows)
    # One representative per 2023 anchor feature and radial alias neighborhood.
    groups=[]
    for i in np.argsort(-scored):
        if any(abs(radii[i]-radii[j])<600 and abs(ph[i]-ph[j])<np.deg2rad(6) for j in groups):continue
        groups.append(i)
    obs=[];stamps=[]
    for j,g in enumerate(fs):
        f=g['f'];b=f['filter'];s=SIGMA[b];data,_=read_product(g['p'],cache,True);scale=1e9*g['pixar_sr']
        for i in groups:
            if not np.isfinite(xy[i,j]).all():continue
            x,y=xy[i,j];flux,err,dx,dy,pk,frac,mad=photometry(data,x,y,s)
            known=[]
            for code,state in g['states'].items():
                if code==899:continue
                dra=(state[0]-g['state'][0]+180)%360-180;l,r,bot,top=g['skybox']
                if not (l-.05<dra<r+.05 and bot-.05<state[1]<top+.05):continue
                xx,yy=g['w'].all_world2pix(*state[:2],0,quiet=True)
                if np.isfinite(xx+yy):known.append((float(np.hypot(x-xx,y-yy)),code))
            near,body=min(known) if known else (np.nan,0)
            obs.append(dict(a_km=radii[i],phase_deg=np.rad2deg(ph[i]),file=f['file'],visit=f['visit'],filter=b,x=x,y=y,good=int(ok[i,j]),
                            score=z[i,j],aperture_mJy=flux*scale,aperture_error_mJy=err*scale,peak_dx=dx,peak_dy=dy,nearest_moon=body,distance_px=near))
            yy,xx=np.mgrid[-15:16,-15:16];grid=[yy+y,xx+x]
            ns,ms=native[j];out=ndi.map_coordinates(ns,grid,order=1);valid=ndi.map_coordinates(ms.astype('f4'),grid,order=1)>.995
            stamps.append((i,f['visit'],b,np.where(valid,out,np.nan)))
    write_csv(OUT/'neptune_alias_representatives.csv',list(obs[0]),obs)
    sets=[max((i for i in groups if radii[i]<180000),key=lambda i:scored[i]),max((i for i in groups if radii[i]>180000),key=lambda i:scored[i])]
    combos=[('01272002001','F2100W'),('01272010001','F2100W'),('01272002001','F2550W'),('01272010001','F2550W'),('01604005001','F1000W')]
    fig,axs=plt.subplots(2,5,figsize=(11,5),layout='constrained',sharex=True,sharey=True)
    for y,i in enumerate(sets):
        for x,(visit_id,b) in enumerate(combos):
            patches=[q for k,vv,bb,q in stamps if k==i and vv==visit_id and bb==b]
            im=np.nanmedian(patches,axis=0);ax=axs[y,x];ax.imshow(im,origin='lower',extent=(-15.5,15.5,-15.5,15.5),vmin=-3,vmax=3,cmap='viridis')
            ax.plot(0,0,'+',color='white',ms=10)
            if y==0:ax.set_title(b+'\n'+{'01272002001':'Dec 2','01272010001':'Dec 5','01604005001':'Nov 2023'}[visit_id])
            if x==0:ax.set_ylabel(f'{radii[i]:,.0f} km\npixels')
    fig.savefig(OUT/'neptune_alias_stamps.png',dpi=140);plt.close(fig)
    print(len(rr),len(groups),[(vname,int(vc[:,j].sum())) for j,vname in enumerate(vs)],int(np.sum(es>=1)),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('headers',type=Path);a=p.parse_args();run(a.cache,a.headers)
