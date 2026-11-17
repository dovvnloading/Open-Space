"""Reproduce nonphysical offset-grid triggers and examine each native exposure."""
import argparse,csv,json,sys
from pathlib import Path
import warnings
import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.ndimage import map_coordinates

def mad(x):
    a=np.nanmedian(x)
    return float(a),float(1.4826*np.nanmedian(np.abs(x-a)))

def load_ledger(wise,name):
    rp=wise/'receipts'/(name+'.json')
    if rp.exists(): return json.loads(rp.read_text())
    import gzip
    for q in (wise/'compact'/'acquisition.jsonl.gz',wise/'acquisition.jsonl.gz'):
        if q.exists():
            return next(i for i in map(json.loads,gzip.open(q,'rt')) if i.get('target')==name or i.get('name')==name)
    raise FileNotFoundError(name)

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--wise',type=Path,required=True);pa.add_argument('--cache',type=Path,required=True)
    pa.add_argument('--out',type=Path,default=Path(__file__).resolve().parent)
    args=pa.parse_args()
    sys.path[:0]=[str(args.wise),str(args.wise.parent/'orbits'),str(args.wise.parent/'other')]
    from search import Frame,positions,offsets,tangent,pix
    import search
    name='RR-h270'
    rec=load_ledger(args.wise,name)
    visit_id={0:1,1:0,2:1}
    frames=[Frame(q,args.cache/name) for q in rec['products'] if 'int_sha256' in q and 'msk_sha256' in q]
    ss=positions(name,frames);ctrl=offsets(ss)
    which={v:[i for i,s in enumerate(frames) if s.visit==v] for v in visit_id}
    def native(tracks):
        visits=[]
        for v,idx in which.items():
            stack=[]
            for i in idx:
                refs=np.array([frames[j].sample(tracks[i]) for j in which[visit_id[v]]])
                stat=np.nanmedian(refs,axis=0)
                stat[np.isfinite(refs).sum(axis=0)<2]=np.nan
                stack.append(frames[i].sample(tracks[i])-stat)
            visits.append(np.array(stack))
        return visits
    def combine(st):
        vs=[]
        for s in st:
            v=np.nanmedian(s,axis=0);v[np.isfinite(s).sum(axis=0)<2]=np.nan
            vs.append(v)
        return np.array(vs)
    c,scatter=np.array([mad(x) for x in combine(native(ctrl))]).T
    measured=native(ss)
    z=(combine(measured)-c[:,None])/scatter[:,None]
    hit=np.flatnonzero(np.isfinite(z).all(axis=0)&(np.min(z,axis=0)>=3)&(z.mean(axis=0)>=5))
    print('hits',hit.tolist(),z[:,hit],flush=True)
    cols=['template','visit','frame','mjd','ra','dec','science_filtered','reference_filtered','difference','z_control']
    out=[]
    grid=[]
    for k in hit:
        entry=dict(template=int(k),east_as=(int(k)%31)*4-60,north_as=(int(k)//31)*4-60,score=z[:,k].tolist(),visits=[])
        for v,idx in which.items():
            coords=np.array([(ss[i][0][k],ss[i][1][k]) for i in idx])
            flux=[]
            for i,(ra,de) in zip(idx,coords):
                f=frames[i]
                coord=(np.array([ra]),np.array([de]))
                science=float(f.sample(coord)[0]);ref_s=float(np.nanmedian([frames[j].sample(coord)[0] for j in which[visit_id[v]]]))
                flux.append(science-ref_s)
                out.append([k,v,f.id,f.mjd,ra,de,science,ref_s,science-ref_s,(science-ref_s-c[v])/scatter[v]])
            pos=int(np.argmin(np.abs([frames[i].mjd-np.median([frames[j].mjd for j in idx]) for i in idx])))
            entry['visits'].append(dict(visit=v,mjd=frames[idx[pos]].mjd,ra=coords[pos,0],dec=coords[pos,1],flux_median_mad=mad(flux)))
        grid.append(entry)
    args.out.mkdir(exist_ok=True,parents=True)
    # Compare the two approximately same-season scans a year apart. The two
    # predicted summer locations are separated by hundreds of arcseconds.
    original=visit_id.copy()
    for alternative in ({0:2,1:0,2:0},{0:2,1:2,2:0}):
        visit_id.update(alternative)
        c2,s2=np.array([mad(x) for x in combine(native(ctrl))]).T
        az=(combine(native(ss))-c2[:,None])/s2[:,None]
        for entry in grid:
            entry.setdefault('alternate_references',[]).append(dict(reference=alternative,
                visit_z=[float(x) if np.isfinite(x) else None for x in az[:,entry['template']]]))
    visit_id.update(original)
    (args.out/'positions.json').write_text(json.dumps(grid,indent=2)+'\n')
    with (args.out/'individual.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(cols);w.writerows(out)

if __name__=='__main__':
    warnings.filterwarnings('ignore',message='All-NaN slice encountered')
    main()
