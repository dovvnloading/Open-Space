"""Off-grid source injections into the registered Neptune MIRI images."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from geometry import generate
from search import (Cuts, TemplateCache, file_hash, load_npz, plants, prepare,
                    scan, summarize_plants, validate_cube, validate_tracks, write_rows)


def run(cube_file, band, scratch, dest):
    dest.mkdir(exist_ok=True,parents=True)
    bank=scratch/f'{band}-annulus.npz';fake=scratch/f'{band}-offgrid.npz'
    generate(cube_file,bank,grid=(180000.,182000.,1000.,2.))
    search_tracks=load_npz(bank)
    rng=np.random.default_rng(2131);n=32
    radii=rng.uniform(180000.,182000.,n);phases=rng.uniform(0,360,n)
    params=dest/'injected_elements.csv'
    with params.open('w') as f:
        w=csv.writer(f,lineterminator='\n');w.writerow(['a_km','e','i_deg','node_deg','peri_deg','mean_deg'])
        for a,p in zip(radii,phases):w.writerow([a,0.,0.,0.,0.,p])
    generate(cube_file,fake,elements=params,model='circular',epoch=float(search_tracks['epoch_jd_tdb']))
    injected=load_npz(fake)
    injected["orbit_id"]=np.array([f"plant-{i:02d}" for i in range(n)])
    c=load_npz(cube_file);validate_tracks(search_tracks,c);validate_tracks(injected,c)
    k=validate_cube(c);cuts=Cuts();weights=prepare(c,k,cuts)
    cache=TemplateCache(c,search_tracks,k)
    base=scan(weights,search_tracks,c['visit'],cuts,cache=cache)
    levels=(.0002,.0006,.002,.006) if band=='F2100W' else (.002,.006,.02,.06)
    records=list(plants(c,k,search_tracks,injected,cuts,base,levels,n,8128,4.,cache=cache))
    write_rows(dest/'plants.csv',records,list(records[0]))
    summary=dict(band=band,trial_domain_km=[180000,182000],phase_domain_deg=[0,360],
                 cube_sha256=file_hash(cube_file),tracks_sha256=file_hash(bank),plants_sha256=file_hash(fake),
                 source='one positive Gaussian PSF per exposure, integrated along off-grid orbit',
                 flux_unit='Jy',matching_radius_pixels=4,results=summarize_plants(records))
    (dest/'depth.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(band,[(r['flux'],r['recovered']) for r in summary['results']],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('band',choices=('F2100W','F2550W'))
    p.add_argument('cache',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();run(a.cache/f'{a.band}.npz',a.band,a.cache,a.output)
