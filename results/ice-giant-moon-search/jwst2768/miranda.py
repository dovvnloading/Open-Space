"""Define and apply a data-derived mask for Miranda's co-moving diffraction pattern."""
import argparse,csv,json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent))
from search import Cuts,file_hash,load_npz,prepare,sample,snr,validate_cube


def run(paths,output,report):
    report=Path(report);report.mkdir(parents=True,exist_ok=True)
    cubes=list(map(load_npz,paths));rows=[];last=0.
    for c in cubes:
        m=prepare(c,validate_cube(c),Cuts(psf_support=.8))
        p=c['known_xy'][:,list(c['known_id']).index(705)]
        theta=np.linspace(0,2*np.pi,72,endpoint=False)[:,None]
        for radius in np.arange(7.,101.):
            n,w=sample(m,p[None,:,0]+radius*np.cos(theta),p[None,:,1]+radius*np.sin(theta))
            all_visits=[]
            for visit in np.unique(c['visit']):
                z=c['visit']==visit
                all_visits.append(np.where((w[:,z]>0).sum(1)>=2,snr(n[:,z].sum(1),w[:,z].sum(1)),0))
            second=np.partition(all_visits,-2,axis=0)[-2]
            if second.max()>=5.:last=max(last,radius)
            rows.append([str(c['band']),radius,*np.percentile(second,[50,90,100])])
    guard=3.;radius=float(guard+last)
    with (report/'radial_profile.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['band','radius_px','median_second_visit_score','p90_second_visit_score','max_second_visit_score']);w.writerows(rows)
    for c,dst in zip(cubes,output):
        xy=c['known_path'][list(c['known_id']).index(705)]
        yy,xx=np.indices(c['image'].shape[1:])
        for j,p in enumerate(xy):
            excluded=np.logical_or.reduce([(xx-x)**2+(yy-y)**2<=radius**2 for x,y in p])
            c['valid'][j,excluded]=False;c['geometry_valid'][j,excluded]=False
        c['miranda_mask_px']=radius
        np.savez_compressed(dst,**c)
    (report/'receipt.json').write_text(json.dumps(dict(input_sha256=list(map(file_hash,paths)),
              output_sha256=list(map(file_hash,output)),mask_radius_px=radius,last_reproducible_5score_px=float(last),
              guard_px=guard,PA_samples=72,search_radius_px=[7,100],visits_required=2,exposures_per_visit=2),indent=2)+'\n')
    print('mask radius',radius)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ['f182','f210','out182','out210','report']:p.add_argument(arg)
    a=p.parse_args();run([a.f182,a.f210],[a.out182,a.out210],a.report)
