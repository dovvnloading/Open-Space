"""Screen disjoint repeated visits that have too few exposures for within-visit dropout."""
import argparse,csv,gzip,json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from screen import ROOT,census,reduce
from stationary import clean


def rank(raw,ok,visits):
    z=[];mask=[]
    vset=np.unique(visits)
    for v in vset:
        m=visits==v;cnt=ok[:,m].sum(1);mask.append(cnt>0)
        z.append(np.where(ok[:,m],np.clip(raw[:,m],-7,7),0).sum(1)/np.sqrt(np.maximum(1,cnt)))
    z=np.asarray(z);mask=np.asarray(mask);n=mask.sum(0)
    score=np.full(len(raw),np.inf)
    valid=n>=4
    for side in [0,1]:
        part=z[side::2];available=mask[side::2].sum(0);valid &= available>=2
        score=np.minimum(score,part.sum(0)/np.sqrt(np.maximum(1,available))/5.)
    # Three separate visits and a stack after dropping the highest-scoring whole visit.
    score=np.minimum(score,np.partition(z,-3,axis=0)[-3]/2.5)
    score=np.minimum(score,(z.sum(0)-z.max(0))/np.sqrt(np.maximum(1,n-1))/5.)
    return np.where(valid,score,0.),valid


def run(args):
    inventory=list(csv.DictReader((ROOT.parent/'manifests/observations.csv').open()))
    prior=list(csv.DictReader((ROOT/(args.planet.lower()+'_coverage.csv')).open()))
    fn=ROOT/(args.planet.lower()+'_sparse.csv')
    done={r['epoch']:r for r in csv.DictReader(fn.open())} if fn.exists() else {}
    with gzip.open(ROOT/(args.planet.lower()+'_positions.csv.gz'),'rt') as f:ep=list(csv.DictReader(f))
    outdir=ROOT/'epochs'
    for old in prior:
        if old['status']=='screened' or old['epoch'] in done:continue
        tag=old['epoch'];band=old['filter'];prop=old['proposal'];mjd=int(tag.split('-')[-1])
        rows=sorted([r for r in inventory if r['planet']==args.planet and r['filters']==band and r['proposal_id']==prop
                     and mjd<=float(r['t_min'])<mjd+4],key=lambda r:float(r['t_min']))
        if len(rows)!=int(old['archive_exposures']):raise ValueError('epoch changed '+tag)
        part=[dict(rootname=r['obs_id'],target=args.planet.upper(),visit=r['obs_id'][4:6],filter=band) for r in rows]
        rec=dict(epoch=tag,archive_exposures=len(part),used_exposures=0,used_visits=0,grid=0,supported=0,
                 screen_hits=0,after_static_subtraction=0,highest_after_static='',status='insufficient_registered_visits')
        out=outdir/tag;out.mkdir(parents=True,exist_ok=True);args.cache.mkdir(parents=True,exist_ok=True)
        with ThreadPoolExecutor(4) as ex:list(ex.map(lambda r:census.download(r,args.cache),part))
        from hashlib import sha256
        hashes=[]
        for r in part:
            name=r['rootname']+'_flt.fits';path=args.cache/name
            hashes.append(dict(file=name,bytes=path.stat().st_size,sha256=sha256(path.read_bytes()).hexdigest()))
        reduce.write_csv(out/'inputs.csv',hashes)
        names={r['rootname'] for r in part}
        mapping={'selection.csv':part,'inputs.csv':hashes,'ephemerides.csv':[r for r in ep if r['rootname'] in names]}
        original=reduce.read_csv;reduce.read_csv=lambda n:mapping[n]
        try:fs,reg,_,_=reduce.prep(args.cache,args.planet,allow_missing_control=True)
        finally:reduce.read_csv=original
        reduce.write_csv(out/'registration.csv',reg)
        vs=np.array([f['visit'] for f in fs]);rec.update(used_exposures=len(fs),used_visits=len(np.unique(vs)))
        if rec['used_visits']>=4:
            bank=reduce.grid(args.planet,reduce.DOMAINS[args.planet]);cleaned=None;visible=triggered=kept=0;max_after=0.
            survivors=[]
            for start in range(0,len(bank),24000):
                sub=bank[start:start+24000];raw,ok=reduce.samples(fs,args.planet,sub);z,coverage=rank(raw,ok,vs)
                visible+=int(coverage.sum());ix=np.flatnonzero(z>=1);triggered+=len(ix)
                if len(ix):
                    if cleaned is None:cleaned=clean(fs)
                    clean_raw,clean_ok=reduce.samples(cleaned,args.planet,sub[ix]);score,_=rank(clean_raw,clean_ok,vs)
                    kept+=int((score>=1).sum());max_after=max(max_after,float(score.max()))
                    for k in np.flatnonzero(score>=1):
                        survivors.append(dict(a_km=sub[ix[k],0],phase_rad=sub[ix[k],1],initial=z[ix[k]],after=score[k]))
            rec.update(grid=len(bank),supported=visible,screen_hits=triggered,after_static_subtraction=kept,
                       highest_after_static=round(max_after,4),status='screened')
            if survivors:reduce.write_csv(out/'sparse_survivors.csv',survivors)
        done[tag]=rec;reduce.write_csv(fn,list(done.values()))
        print(tag,rec['status'],rec['after_static_subtraction'],flush=True)
        for r in hashes:(args.cache/r['file']).unlink(missing_ok=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('planet',choices=list(reduce.DOMAINS))
    p.add_argument('cache',type=Path);run(p.parse_args())
