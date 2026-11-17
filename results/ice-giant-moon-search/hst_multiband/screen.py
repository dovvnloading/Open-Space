"""Screen single-epoch, single-filter WFC3 observations, with bounded local cache."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import json
import sys
from pathlib import Path
import time
import numpy as np
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'hst16995'))
sys.path.insert(0,str(ROOT.parent/'hst_f845m'))
import census
import reduce


def run(args):
    cache=args.cache;cache.mkdir(parents=True,exist_ok=True)
    inventory=list(csv.DictReader((ROOT.parent/'manifests/observations.csv').open()))
    rows=[x for x in inventory if x['obs_collection']=='HST' and x['planet']==args.planet]
    with gzip.open(ROOT/(args.planet.lower()+'_positions.csv.gz'),'rt') as f: eph=list(csv.DictReader(f))
    ep_by_name={}
    for r in eph: ep_by_name.setdefault(r['rootname'],[]).append(r)
    summary_file=ROOT/(args.planet.lower()+'_coverage.csv')
    summary={r['epoch']:r for r in csv.DictReader(summary_file.open())} if summary_file.exists() else {}
    names=['F845M','F467M','F763M','F547M','F657N','F390M','F487N','FQ727N','FQ619N']
    for band in names:
        for prop in sorted({r['proposal_id'] for r in rows if r['filters']==band}):
            chosen=sorted([r for r in rows if r['filters']==band and r['proposal_id']==prop],key=lambda r:float(r['t_min']))
            chunks=np.split(np.arange(len(chosen)),np.flatnonzero(np.diff([float(r['t_min']) for r in chosen])>3.)+1)
            for inds in chunks:
                obs=[chosen[k] for k in inds];start=obs[0]['t_min'].split('.')[0]
                tag=f'{args.planet.lower()}-{prop}-{band.lower()}-{start}'
                if tag in summary:continue
                part=[dict(rootname=r['obs_id'],target=args.planet.upper(),visit=r['obs_id'][4:6],filter=band) for r in obs]
                counts={v:sum(r['visit']==v for r in part) for v in {r['visit'] for r in part}}
                rec=dict(epoch=tag,filter=band,proposal=prop,archive_exposures=len(obs),archive_visits=len(counts),used_exposures=0,
                         used_visits=0,grid=0,supported=0,screen_hits=0,dropout_hits=0,top_rank='',status='insufficient_repeated_visits')
                if sum(n>=2 for n in counts.values())>=4:
                    out=ROOT/'epochs'/tag;out.mkdir(exist_ok=True,parents=True)
                    existing=args.shared if band=='F845M' else cache
                    try:
                        with ThreadPoolExecutor(4) as ex:list(ex.map(lambda r:census.download(r,existing),part))
                        hashes=[]
                        for r in part:
                            name=r['rootname']+'_flt.fits';path=existing/name
                            hashes.append(dict(file=name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size))
                        reduce.write_csv(out/'inputs.csv',hashes)
                        mapping={'selection.csv':part,'inputs.csv':hashes,'ephemerides.csv':[e for r in part for e in ep_by_name[r['rootname']]]}
                        original=reduce.read_csv
                        reduce.read_csv=lambda n:mapping[n]
                        try:
                            reduce.run(argparse.Namespace(out=out,cache=existing,planet=args.planet,band=band,no_plants=True))
                        finally:reduce.read_csv=original
                        s=json.loads((out/'summary.json').read_text())
                        rec.update(used_exposures=s['exposures'],used_visits=s['visits'],grid=s['grid_tracks'],supported=s['supported_tracks'],
                                   screen_hits=s['threshold_crossings'],dropout_hits=s['dropout_crossings'],top_rank=round(s['top_dropout_rank'],4),status='screened')
                    except ValueError as e:
                        if str(e)!='No isolated control PSF':raise
                        rec['status']='no_registered_control_PSF'
                    if not args.retain and existing==cache:
                        for r in part:(cache/(r['rootname']+'_flt.fits')).unlink(missing_ok=True)
                summary[tag]=rec
                reduce.write_csv(summary_file,list(summary.values()))
                print(tag,rec['status'],rec['dropout_hits'],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('planet',choices=list(reduce.DOMAINS))
    ap.add_argument('cache',type=Path);ap.add_argument('shared',type=Path,help='directory of reusable F845M FITS')
    ap.add_argument('--retain',action='store_true');run(ap.parse_args())
