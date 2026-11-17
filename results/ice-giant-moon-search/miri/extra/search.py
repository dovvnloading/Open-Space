"""Native MIRI multi-filter orbit screening with propagated, observer-specific astrometry."""
import argparse,csv,json,warnings
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from fetch import read_product,products
from coverage import OUT,frame,get_states,pixels
from geometry import AU,ARCSEC
from miri_audit import write_csv

SIGMA={'F560W':.75,'F1000W':1.3,'F1280W':1.75,'F1800W':2.4,'F2100W':2.8,'F2550W':3.3}


def maps(g,p,cache):
    dat,_=read_product(p,cache,True);val=dat['SCI'].astype('f4');err=dat['ERR'];dq=dat['DQ'];sig=SIGMA[g['f']['filter']]
    ok=np.isfinite(val)&np.isfinite(err)&(err>0)&((dq==0)|(dq==4))
    yy,xx=np.indices(val.shape); pix=.11
    for code,r in g['states'].items():
        dra=(r[0]-g['state'][0]+180)%360-180;l,rr,b,t=g['skybox']
        if not (l-.003<dra<rr+.003 and b-.003<r[1]<t+.003):continue
        x,y=g['w'].all_world2pix(r[0],r[1],0,quiet=True)
        if not np.isfinite(x+y):continue
        if code in (799,899):
            # Conservative circular sky exclusion encloses the disk and entire regular ring/arc system.
            rad=(110000 if code==799 else 80000)/(r[4]*AU)*ARCSEC/pix+4*sig
        else:
            moon_rad=r[2]/(2*pix) if np.isfinite(r[2]) else 0
            rad=max(7,5*sig+2,moon_rad+(9 if code in (701,702,703,704,705,801) else 5)*sig)
        ok &= np.hypot(xx-x,yy-y)>rad
    val=np.where(ok,val,0);w=ok.astype('f4')
    def blur(sigma):
        den=ndi.gaussian_filter(w,sigma)
        return ndi.gaussian_filter(val,sigma)/np.maximum(den,.001),den
    core,den=blur(sig);back,outer=blur(max(5,4.5*sig))
    valid=ok&(den>.75)&(outer>.7)&(ndi.distance_transform_edt(np.pad(np.ones(val.shape,bool),1))[1:-1,1:-1]>3*sig)
    res=core-back;ny,nx=val.shape
    # Estimate local background and matched-filter scatter on 32x32 blocks; a scene-wide floor
    # prevents anomalously flat blocks from generating large scores.
    block=32;nbx=(nx+block-1)//block;nby=(ny+block-1)//block
    bg=np.full((nby,nbx),np.nan);noise=bg.copy()
    for by,bx in np.ndindex(bg.shape):
        q=np.s_[by*block:(by+1)*block,bx*block:(bx+1)*block];v=res[q][valid[q]]
        if len(v)<80:continue
        med=np.median(v);bg[by,bx]=med;noise[by,bx]=1.4826*np.median(abs(v-med))
    finite=np.isfinite(bg)&(noise>0)
    if not finite.any():return np.zeros(val.shape,'f4'),np.zeros(val.shape,bool)
    nearest=ndi.distance_transform_edt(~finite,return_distances=False,return_indices=True);bg=bg[tuple(nearest)];noise=noise[tuple(nearest)]
    rc=np.array([(yy+.5)/block-.5,(xx+.5)/block-.5]);surface=ndi.map_coordinates(bg,rc,order=1,mode='nearest')
    scale=np.maximum(ndi.map_coordinates(noise,rc,order=1,mode='nearest'),np.percentile(noise[finite],10))
    return np.where(valid,(res-surface)/scale,0).astype('f4'),valid


def visit_scores(raw,support,frames,planet):
    visits=sorted(set(g['f']['visit'] for g in frames));ranks=[];covers=[]
    for v in visits:
        select=np.flatnonzero([g['f']['visit']==v for g in frames]);z=np.clip(raw[:,select],-7,7);ok=support[:,select];n=ok.sum(1)
        z=np.where(ok,z,0);total=z.sum(1)/np.sqrt(np.maximum(n,1))
        if len(select)>=6:
            # Exposure order alternation spreads each half across the association.
            ev=ok[:,::2].sum(1);od=ok[:,1::2].sum(1);cover=(ev>=2)&(od>=2)
            rank=np.minimum.reduce([total/6,z[:,::2].sum(1)/np.sqrt(np.maximum(1,ev))/3.2,z[:,1::2].sum(1)/np.sqrt(np.maximum(1,od))/3.2])
        elif len(select)>=3:
            cover=n>=2
            trim=(z.sum(1)-np.max(np.where(ok,z,-100),axis=1))/np.sqrt(np.maximum(1,n-1))
            rank=np.minimum(total/5.,trim/2.5)
        else:raise ValueError('two exposures need visit-matched combined-filter data')
        ranks.append(np.where(cover,rank,0));covers.append(cover)
    ranks=np.column_stack(ranks);covers=np.column_stack(covers)
    if planet=='Uranus':
        rank=np.sort(ranks,axis=1)[:,-2];good=covers.sum(1)>=2
    else:
        anchor=visits.index('01604005001');others=[j for j in range(len(visits)) if j!=anchor]
        rank=np.minimum(ranks[:,anchor],np.max(ranks[:,others],axis=1));good=covers[:,anchor]&covers[:,others].any(1)
    return rank,good,ranks,covers,visits


def run(planet,cache,heads,grid):
    start,stop,step,dphase=grid
    states=get_states();frames=[]
    for f,p in products():
        if f['planet']!=planet or (planet=='Neptune' and f['planet_on_chip']!='1'):continue
        g=frame(f,p,heads,states[f['file']],{})
        g['score'],g['ok']=maps(g,p,cache);frames.append(g)
    frames.sort(key=lambda g:g['t'])
    # The cross-year search freezes phase at the single 1604 Neptune visit.
    origin=min(g['t'] for g in frames if planet=='Uranus' or g['f']['visit']=='01604005001')
    qc=[dict(file=g['f']['file'],visit=g['f']['visit'],filter=g['f']['filter'],detector_pixels=g['ok'].size,
             usable_pixels=int(g['ok'].sum()),fraction=float(g['ok'].mean())) for g in frames]
    prefix=planet.lower()+f'_{start:g}_{stop:g}';write_csv(OUT/(prefix+'_images.csv'),list(qc[0]),qc)
    az=np.arange(0,360,dphase);rad=np.arange(start,stop+step*.1,step)
    rows=[];leads=[];top=[]
    for a in np.array_split(rad,max(1,(len(rad)+11)//12)):
        rr=np.repeat(a,len(az));ph=np.tile(np.deg2rad(az),len(a));raw=np.zeros((len(rr),len(frames)),'f4');support=np.zeros_like(raw,bool)
        for j,g in enumerate(frames):
            ix,x,y=pixels(g,rr,ph,origin)
            if not len(ix):continue
            xy=[y,x];support[ix,j]=ndi.map_coordinates(g['ok'].astype('f4'),xy,order=1,mode='constant')>.995
            raw[ix,j]=ndi.map_coordinates(g['score'],xy,order=1,mode='constant')
        score,valid,z,cover,visits=visit_scores(raw,support,frames,planet)
        neg,_,_,_,_=visit_scores(-raw,support,frames,planet)
        for r in a:
            s=rr==r
            rows.append(dict(a_km=r,templates=int(s.sum()),eligible=int(valid[s].sum()),passing=int((score[s]>=1).sum()),negative_passing=int((neg[s]>=1).sum()),
                             max_rank=float(score[s].max()),max_negative_rank=float(neg[s].max())))
        for i in np.flatnonzero(score>=1):
            leads.append(dict(a_km=rr[i],phase_deg=np.rad2deg(ph[i]),rank=score[i],**{f'visit_{v}':m for v,m in zip(visits,z[i])}))
        for i in np.argsort(score)[-2:][::-1]:
            top.append(dict(a_km=rr[i],phase_deg=np.rad2deg(ph[i]),rank=score[i],**{f'visit_{v}':m for v,m in zip(visits,z[i])}))
    write_csv(OUT/(prefix+'_radii.csv'),list(rows[0]),rows);write_csv(OUT/(prefix+'_leads.csv'),list(top[0]),leads)
    write_csv(OUT/(prefix+'_top.csv'),list(top[0]),sorted(top,key=lambda r:-r['rank'])[:12])
    summary=dict(planet=planet,grid=grid,epoch_jd_tdb=origin,visits=visits,templates=sum(r['templates'] for r in rows),eligible=sum(r['eligible'] for r in rows),
                 passing=sum(r['passing'] for r in rows),inverse_passing=sum(r['negative_passing'] for r in rows),max_rank=max(r['max_rank'] for r in rows))
    (OUT/(prefix+'_summary.json')).write_text(json.dumps(summary,indent=2)+'\n');print(summary,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('planet',choices=['Uranus','Neptune']);p.add_argument('cache',type=Path)
    p.add_argument('headers',type=Path);p.add_argument('--grid',nargs=4,type=float,required=True);a=p.parse_args();run(a.planet,a.cache,a.headers,a.grid)
