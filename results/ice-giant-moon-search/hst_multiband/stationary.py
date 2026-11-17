"""Recheck every triggered multi-epoch template after subtracting planet-fixed structure."""
import argparse
import csv
import gzip
import json
from pathlib import Path
import sys
from concurrent.futures import ThreadPoolExecutor
import warnings
import numpy as np
from scipy import ndimage as ndi
from screen import ROOT,census,reduce


def clean(fs):
    # For a given exposure, all images in that visit are held out of the static reference.
    # The projected equatorial basis maps the detector offsets between the images.
    result=[]
    for f in fs:
        yy,xx=np.indices(f['score'].shape)
        uv=np.linalg.solve(f['project'],np.stack([(xx-f['c'][0]).ravel(),(yy-f['c'][1]).ravel()]))
        references=[]
        for v in np.unique([g['visit'] for g in fs]):
            if v==f['visit']: continue
            images=[]
            for g in fs:
                if g['visit']!=v:continue
                xy=(g['project']@uv+g['c'][:,None]).reshape(2,*xx.shape)
                support=ndi.map_coordinates(g['support'].astype('f4'),xy[::-1],order=1,mode='constant')>.98
                val=ndi.map_coordinates(g['score'],xy[::-1],order=1,mode='constant')
                images.append(np.where(support,val,np.nan))
            references.append(np.nanmedian(images,axis=0))
        count=np.isfinite(references).sum(axis=0)
        ref=np.nanmedian(references,axis=0);valid=f['support']&(count>=2)
        res=f['score']-ref
        # Correct the extra noise from estimating the template, preserving the original thresholds.
        bins=(np.hypot(xx-f['c'][0],yy-f['c'][1])/12).astype(int)
        for b in np.unique(bins[valid]):
            use=(bins==b)&valid
            if use.sum()<80:valid[use]=False;continue
            med=np.median(res[use]);scale=max(1.,1.4826*np.median(abs(res[use]-med)))
            res[use]=(res[use]-med)/scale
        result.append(dict(f,score=np.where(valid,res,0),support=valid))
    return result


def audit(a):
    inventory=list(csv.DictReader((ROOT.parent/'manifests/observations.csv').open()))
    path=ROOT/(a.planet.lower()+'_vet.csv')
    audited={r['epoch']:r for r in csv.DictReader(path.open())} if path.exists() else {}
    src=list(csv.DictReader((ROOT/(a.planet.lower()+'_coverage.csv')).open()))
    with gzip.open(ROOT/(a.planet.lower()+'_positions.csv.gz'),'rt') as fi: eph=list(csv.DictReader(fi))
    for row in src:
        tag=row['epoch']
        if int(row['dropout_hits'])==0 or tag in audited:continue
        out=ROOT/'epochs'/tag
        hashes=list(csv.DictReader((out/'inputs.csv').open()))
        names={r['file'].split('_')[0] for r in hashes}
        rows=sorted([r for r in inventory if r['obs_id'] in names],key=lambda r:float(r['t_min']))
        part=[dict(rootname=r['obs_id'],target=a.planet.upper(),visit=r['obs_id'][4:6],filter=row['filter']) for r in rows]
        cache=a.shared if row['filter']=='F845M' else a.cache;cache.mkdir(parents=True,exist_ok=True)
        with ThreadPoolExecutor(4) as ex:list(ex.map(lambda r:census.download(r,cache),part))
        mapping={'selection.csv':part,'inputs.csv':hashes,'ephemerides.csv':[e for e in eph if e['rootname'] in names]}
        original=reduce.read_csv;reduce.read_csv=lambda n:mapping[n]
        try:fs,reg,_,_=reduce.prep(cache,a.planet)
        finally:reduce.read_csv=original
        other=clean(fs);visits=np.array([f['visit'] for f in fs]);bank=reduce.grid(a.planet,reduce.DOMAINS[a.planet])
        survivors=[];nraw=0;highest=0.;best_records=[]
        for start in range(0,len(bank),24000):
            sub=bank[start:start+24000];raw,ok=reduce.samples(fs,a.planet,sub);score=reduce.score(raw,ok,visits,True)[0]
            hits=np.flatnonzero(score>=1)
            if len(hits)==0:continue
            nraw+=len(hits)
            raw2,ok2=reduce.samples(other,a.planet,sub[hits]);again=reduce.score(raw2,ok2,visits,True)[0]
            highest=max(highest,float(again.max()))
            for i in np.argsort(again)[-3:]:
                best_records.append(dict(a_km=sub[hits[i],0],phase_rad=sub[hits[i],1],original_rank=score[hits[i]],clean_rank=again[i]))
            for j in np.flatnonzero(again>=1):
                survivors.append(dict(a_km=sub[hits[j],0],phase_rad=sub[hits[j],1],original_rank=score[hits[j]],clean_rank=again[j]))
        if nraw!=int(row['dropout_hits']):raise AssertionError((tag,nraw,row['dropout_hits']))
        best_records.sort(key=lambda r:-r['clean_rank']);reduce.write_csv(out/'static_vet_top.csv',best_records[:32])
        if survivors:reduce.write_csv(out/'static_vet_survivors.csv',survivors)
        audited[tag]=dict(epoch=tag,raw_hits=nraw,after_static_subtraction=len(survivors),highest_rank=round(highest,4))
        reduce.write_csv(path,list(audited.values()));print(tag,nraw,len(survivors),flush=True)
        if row['filter']!='F845M':
            for r in hashes:(cache/r['file']).unlink(missing_ok=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('planet',choices=list(reduce.DOMAINS))
    p.add_argument('cache',type=Path);p.add_argument('shared',type=Path);audit(p.parse_args())
