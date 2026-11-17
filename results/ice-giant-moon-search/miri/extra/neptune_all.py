"""Use the 1604 visit to select fine annual-alias trials for both independent 1272 visits."""
import argparse,csv,json
from pathlib import Path
import numpy as np
from scipy import ndimage as ndi
from coverage import OUT,frame,get_states,pixels
from fetch import products
from native_search import maps,visit_scores
from depth import get_xy,sample
from miri_audit import write_csv


def run(cache,head):
    states=get_states();groups=[]
    for f,p in products():
        if f['planet']!='Neptune' or f['planet_on_chip']!='1':continue
        g=frame(f,p,head,states[f['file']],{});g['score'],g['ok']=maps(g,p,cache);groups.append(g)
    groups.sort(key=lambda g:g['t']);current=[g for g in groups if g['f']['visit']=='01604005001'];t0=min(g['t'] for g in current)
    def calc(frames,r,ph):
        xy=get_xy(frames,r,np.deg2rad(ph),t0)
        return sample(frames,[(g['score'],g['ok']) for g in frames],xy)
    def anchor_rank(raw,sup):
        z=np.where(sup,np.clip(raw,-7,7),0);n=sup.sum(1);tot=z.sum(1)
        return np.where(n>=2,np.minimum(tot/np.sqrt(np.maximum(n,1))/5,(tot-np.max(np.where(sup,z,-100),axis=1))/np.sqrt(np.maximum(n-1,1))/2.5),0)
    coarse_radii=np.arange(50000,300001,200);phi=np.arange(0,360,.25);cells=[];summary=[]
    for radius in np.array_split(coarse_radii,64):
        r=np.repeat(radius,len(phi));p=np.tile(phi,len(radius));raw,sup=calc(current,r,p)
        pos=anchor_rank(raw,sup);neg=anchor_rank(-raw,sup)
        for sign,score in [(1,pos),(-1,neg)]:
            for rr,pp,rank in zip(r[score>=.70],p[score>=.70],score[score>=.70]):cells.append((sign,rr,pp,rank))
        for rr in radius:
            keep=r==rr
            summary.append(dict(a_km=rr,phases=len(phi),eligible_anchor=int((sup[keep].sum(1)>=2).sum()),
                                anchor_070=int((pos[keep]>=.70).sum()),negative_anchor_070=int((neg[keep]>=.70).sum())))
    write_csv(OUT/'neptune_anchor_radii.csv',list(summary[0]),summary)
    print('seeds',len(cells),flush=True)
    # Subdivide every selected 200-km x 0.25-deg cell and include both adjacent phases.
    leads=[];counts=[];maxima=[]
    for sign in [1,-1]:
        keys={}
        for s,r,p,z in cells:
            if s==sign:
                for q in ((p-.125)%360,p,(p+.125)%360):keys.setdefault(float(q),set()).add(int(r))
        for p,parents in sorted(keys.items()):
            rs=set()
            for a in parents:
                spacing=.5 if a<=130000 else 1. if a<=190000 else 2.
                rs.update(np.arange(max(50000,a-100),min(300000,a+100)+.01,spacing).tolist())
            rs=np.array(sorted(rs));nvalid=npass=0;top=0
            for r in np.array_split(rs,max(1,len(rs)//12000+1)):
                raw,sup=calc(groups,r,np.full(len(r),p));scores,valid,components,_,vs=visit_scores(sign*raw,sup,groups,'Neptune')
                nvalid+=int(valid.sum());npass+=int(np.sum(scores>=1));top=max(top,float(scores.max()))
                for i in np.flatnonzero(scores>=1):
                    leads.append(dict(sign=sign,a_km=r[i],phase_deg=p,rank=scores[i],**{f'visit_{v}':m for v,m in zip(vs,components[i])}))
            counts.append(dict(sign=sign,phase_deg=p,fine_trials=len(rs),eligible=nvalid,passing=npass,highest_rank=top))
    head=['sign','a_km','phase_deg','rank',*[f'visit_{v}' for v in sorted({g['f']['visit'] for g in groups})]]
    write_csv(OUT/'neptune_all_leads.csv',head,leads);write_csv(OUT/'neptune_all_fine.csv',list(counts[0]),counts)
    report=dict(coarse_templates=len(coarse_radii)*len(phi),anchor_supported=sum(r['eligible_anchor'] for r in summary),epoch_jd_tdb=t0,
                branches={s:dict(fine_trials=sum(q['fine_trials'] for q in counts if q['sign']==s),eligible=sum(q['eligible'] for q in counts if q['sign']==s),
                                 passing=sum(q['passing'] for q in counts if q['sign']==s),max_rank=max(q['highest_rank'] for q in counts if q['sign']==s)) for s in [1,-1]})
    (OUT/'neptune_all_summary.json').write_text(json.dumps(report,indent=2)+'\n');print(report,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('headers',type=Path);a=p.parse_args();run(a.cache,a.headers)
