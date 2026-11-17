"""Create the deposited compact strip receipts, trigger table and recovery records."""
import argparse
import csv
import json
from pathlib import Path
from search import HERE
from photometry import save_csv


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--profiles',type=Path,default=HERE/'profiles')
    p.add_argument('--ray-tests',type=Path,required=True)
    args=p.parse_args()
    strips=[];triggers=[];inj=[];ray=[];seeds=[]
    for file in sorted((HERE/'runs').glob('*.json')):
        a=json.loads(file.read_text())
        b=json.loads((args.profiles/file.name).read_text())
        assert len(a['triggers'])==len(b['profiles'])
        info={k:v for k,v in a.items() if k not in ['triggers','injections']}
        info['n_triggers']=len(b['profiles'])
        info['n_profiles_passing']=sum(v['pass_vis'] for v in b['profiles'])
        strips.append(info)
        for row in b['profiles']:
            triggers.append(row)
            if row['pass_vis']:seeds.append(row['id'])
        for row in b['injections']:
            inj.append(dict(tile=a['tile'],row=a['row'],**row))
    for file in sorted(args.ray_tests.glob('*.csv')):
        with file.open() as f:ray.extend(csv.DictReader(f))
    if len(strips)!=32 or not ray:
        raise RuntimeError('run all 32 strips and population injections before exporting')
    assert set(v['tile'] for v in ray)==set(v['tile'] for v in strips)
    save_csv(HERE/'strips.csv',strips)
    save_csv(HERE/'vis_triggers.csv',triggers)
    save_csv(HERE/'ring_injections.csv',inj)
    save_csv(HERE/'lens_injections.csv',ray)
    print('eligible center area [deg2]',sum(v['eligible_binned_centers'] for v in strips)*.2**2/3600**2)
    print('preliminary, passing, injections',len(triggers),len(seeds),len(inj),len(ray))
if __name__=='__main__':main()
