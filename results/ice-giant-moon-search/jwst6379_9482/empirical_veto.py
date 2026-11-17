"""Apply an independent known-moon proximity veto to completed links and planted orbits."""
import argparse
import csv,json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
from scipy.spatial import cKDTree
import numpy as np
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base/'jwst9482'))
from link_visits import geometry,search,trials


def proximity(cube_paths,ephemeris,elems,epoch,tmp):
    best=np.full(len(elems),np.inf);body=np.zeros(len(elems),int)
    for v,path in enumerate(cube_paths):
        c=search.load_npz(path)
        # Nominal JPL moon trails plus measured two-visit Perdita track.
        cp=trials(tmp/'e.csv',path,tmp/'perdita.npz',np.array([[ephemeris['a_km'],0,0,0,0,ephemeris['phase_late_deg']]]),2460955.765061979)
        empirical=np.stack([cp['x_path'][0],cp['y_path'][0]],axis=-1)
        trees=[]
        for f in range(len(c['image'])):
            s=np.concatenate([c['known_path'][:,f],empirical[None,f]],axis=0)
            ids=np.repeat(np.r_[c['known_id'],-725],s.shape[1])
            trees.append((cKDTree(s.reshape(-1,2)),ids))
        for off in range(0,len(elems),300):
            selected=np.arange(off,min(off+300,len(elems)))
            b=trials(tmp/'e.csv',path,tmp/'g.npz',elems[selected],epoch)
            for f,(tree,ids) in enumerate(trees):
                xy=np.stack([b['x_path'][:,f],b['y_path'][:,f]],axis=-1)
                dist,index=tree.query(xy)
                arg=np.argmin(dist,axis=1);di=dist[np.arange(len(arg)),arg];ind=index[np.arange(len(arg)),arg]
                take=di<best[selected]
                best[selected[take]]=di[take];body[selected[take]]=ids[ind[take]]
    return best,body


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cube1','cube2','links','plants','empirical','out']:p.add_argument(k,type=Path)
    p.add_argument('--pad',type=float,default=10)
    a=p.parse_args();a.out.mkdir(exist_ok=False,parents=True)
    e=json.loads(a.empirical.read_text())
    links=list(csv.DictReader((a.links/'passing_pairs.csv').open()))
    plants=list(csv.DictReader((a.plants/'plants.csv').open()))
    ip=json.loads((a.plants/'report.json').read_text())
    epoch=json.loads((a.links/'parameters.json').read_text())['epoch1_tdb_jd']
    le=np.zeros((len(links),6));le[:,0]=[float(r['a_km']) for r in links];le[:,-1]=[float(r['phase_first_deg']) for r in links]
    rng=np.random.default_rng(ip['seed']);pe=np.zeros((120,6));pe[:,0]=rng.uniform(79500,80500,120);pe[:,-1]=rng.uniform(0,360,120)
    with TemporaryDirectory() as td:
        pt=Path(td)
        ld,lb=proximity([a.cube1,a.cube2],e,le,epoch,pt)
        pd,pb=proximity([a.cube1,a.cube2],e,pe,ip['epoch_jd_tdb'],pt)
    counts={}
    with (a.out/'link_distances.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['input_row','early_peak','late_peak','early_sign','late_sign','min_trail_sep_px','body','veto'])
        for i,(r,d,b) in enumerate(zip(links,ld,lb,strict=True)):
            v=bool(d<=a.pad)
            w.writerow([i,r['early_peak'],r['late_peak'],r['early_sign'],r['late_sign'],d,b,v])
            key=f"{r['early_sign']}/{r['late_sign']}";g=counts.setdefault(key,[0,0]);g[v]+=1
    kept=[];retained=[]
    with (a.out/'plant_dispositions.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['flux','plant_index','original_recovered','min_trail_sep_px','body','guaranteed_clear','retained'])
        for r in plants:
            i=int(r['plant_index']);clear=pd[i]>a.pad+ip['matched_tolerance_px']
            passed=r['recovered']=='True' and clear
            w.writerow([r['flux'],i,r['recovered'],pd[i],pb[i],clear,passed])
            row=dict(flux=float(r['flux']),status='passed' if passed else 'missed',recovered=bool(passed))
            retained.append(row)
            if clear:kept.append(row)
    rpt=dict(padding_px=a.pad,truth_matching_reserve_px=ip['matched_tolerance_px'],
        definition='Any encounter during any CAL exposure vetoes the complete orbit. ID -725 is empirical Perdita. Other IDs are nominal Horizons/SPK moon trails.',
        links_by_sign_unvetoed_vetoed=counts,random_population=len(pe),
        random_population_guaranteed_clear=int((pd>a.pad+ip['matched_tolerance_px']).sum()),
        completeness_on_original_observable_population=search.summarize_plants(retained),
        completeness_on_observable_and_guaranteed_clear=search.summarize_plants(kept),
        inputs_sha256={str(f):search.file_hash(f) for f in [a.cube1,a.cube2,a.links/'passing_pairs.csv',a.plants/'plants.csv',a.empirical,Path(__file__)]})
    (a.out/'summary.json').write_text(json.dumps(rpt,indent=2)+'\n')


if __name__=='__main__':main()
