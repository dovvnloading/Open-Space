"""Select, acquire and search one latest public FFI visit per matched non-S109 bright source."""
import argparse
import csv
import hashlib
import gzip
import json
import math
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from astropy.time import Time
from tess_stars2px import TESS_Spacecraft_Pointing_Data as Pointings

from ffi_acquire import (GAIA_FIELDS,TIC_FIELDS,TESSCUT,download,propagated,tsv_rows,vizier)
from ffi_photometry import load_pixel_tools,run
from wd_search import search

VISIT_FIELDS = ['target','path','cadence_s','light_fraction','mass_solar','radius_solar']


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows, fields):
    with path.open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fields,lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(args):
    cat=tsv_rows(args.parent)
    matches={r['gaia']:r for r in read_csv(args.matches)}
    primary={r['GaiaEDR3'] for r in read_csv(args.primary)}
    rows=[]
    for raw in cat:
        source=raw['GaiaEDR3']
        m=matches[source]
        star=dict(raw,tic_requested=m['tic'],selection='')
        try:
            mass=float(raw['MassH'])
            radius=math.sqrt(1.3271244e20*mass/(10**float(raw['loggH'])*.01))/6.957e8
        except (ValueError,OverflowError):
            mass,radius=None,None
        star.update(mass_solar=mass,radius_solar=radius)
        if source in primary:
            status='primary geometric roster'
        elif m['status']!='source_id_and_position':
            status='no unique source ID and position TIC'
        elif m['tic'] in {'267574918','60040774','219244444'}:
            status='positive control'
        elif (mass is None or not .1<=mass<=1.44 or not .002<=radius<=.2):
            status='stellar fit unavailable or outside cut'
        else:
            status='request latest visit'
        star['selection']=status
        rows.append(star)
    out=args.record/'extension_parent.csv'
    write_csv(out,rows,rows[0].keys())
    report=dict(catalogue_sha256=sha(args.parent),tic_matches_sha256=sha(args.matches),
                s109_geometric_sha256=sha(args.primary),roster_sha256=sha(out),
                decisions=Counter(r['selection'] for r in rows))
    (args.record/'extension_selection.json').write_text(json.dumps(report,indent=2)+'\n')
    print(report)


def obtain(row, cache, curves, px):
    source=row['GaiaEDR3']
    output=curves/(source+'.json')
    if output.exists():
        return json.loads(output.read_text())
    ra,dec=propagated(row,2025.)
    coords=cache/f'{source}-sectors.json'
    url='https://mast.stsci.edu/tesscut/api/v0.1/sector?'+urlencode(dict(ra=f'{ra:.8f}',dec=f'{dec:.8f}'))
    download(coords,url)
    hits=json.loads(coords.read_text())['results']
    trials=[]
    result=dict(source_id=source,wd=row['WDJname'],status='no available 200-second visit',attempts=trials)
    for hit in sorted(hits,key=lambda h:-int(h['sector'])):
        sector=int(hit['sector'])
        if not 56<=sector<=108:
            continue
        epoch=float(Time(Pointings.midtimes[Pointings.sectors==sector][0],format='jd').jyear)
        ra,dec=propagated(row,epoch)
        star=dict(row,sector=sector,camera=hit['camera'],ccd=hit['ccd'],coordinate_epoch=epoch)
        path=cache/f'{source}-s{sector}.fits'
        download(path,TESSCUT+urlencode(dict(ra=f'{ra:.8f}',dec=f'{dec:.8f}',sector=sector,x=9,y=9,units='px')))
        trials.append(dict(sector=sector,public_fits=path.exists(),receipt=path.name+'.json'))
        if not path.exists():
            continue
        for tag,tab,fields,center,rad in [
            ('gaia','I/355/gaiadr3',GAIA_FIELDS,(float(row['RA_ICRS']),float(row['DE_ICRS'])),180),
            ('tic82','IV/39/tic82',TIC_FIELDS,propagated(row,2000.),30)]:
            a,d=center
            download(cache/f'{source}-{tag}.tsv',vizier(tab,fields,**{'-c':f'{a:.9f} {d:+.9f}','-c.rs':str(rad)}))
        try:
            measurement=run(star,cache,curves,px)
        except ValueError as exc:
            if 'not 200-second FFI' in str(exc):
                continue
            raise
        result.update(measurement)
        result['coordinate_epoch']=epoch
        if result['status']=='included':
            if str(result['tic'])!=row['tic_requested']:
                raise ValueError('TIC crossmatches disagree')
            result.update(mass_solar=float(row['mass_solar']),radius_solar=float(row['radius_solar']))
            with np.load(curves/result['lc_file']) as f:
                result['search']=search(f['time'],f['flux'],float(f['cadence_s']))
        break
    output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('COMPLETE',source,result['status'],len(result.get('search',{}).get('periodic',[])),len(result.get('search',{}).get('isolated',[])),flush=True)
    return result


def acquire(args):
    roster=read_csv(args.record/'extension_parent.csv')
    rows=[r for r in roster if r['selection']=='request latest visit']
    px=load_pixel_tools(args.pixel_tools)
    args.cache.mkdir(parents=True,exist_ok=True)
    args.curves.mkdir(parents=True,exist_ok=True)
    def process(row):
        try:
            return obtain(row,args.cache,args.curves,px)
        except Exception as exc:
            print('FAILED',row['GaiaEDR3'],str(exc),flush=True)
            return dict(source_id=row['GaiaEDR3'],status='acquisition/measurement error',message=str(exc))
    with ThreadPoolExecutor(args.workers) as pool:
        results=list(pool.map(process,rows))
    (args.record/'extension_measurements.json').write_text(json.dumps(results,indent=1,allow_nan=False)+'\n')
    visits=[dict(target=r['source_id'],path=r['lc_file'],cadence_s=r['cadence_s'],
                 light_fraction=r['target_fraction'],mass_solar=r['mass_solar'],radius_solar=r['radius_solar'])
            for r in results if r['status']=='included']
    write_csv(args.curves/'extension_visits.csv',visits,VISIT_FIELDS)
    write_csv(args.record/'extension_visits.csv',visits,VISIT_FIELDS)
    receipts=[dict(file=p.name,**json.loads(p.read_text())) for p in sorted(args.cache.glob('*'))
              if p.name.endswith(('.fits.json','.tsv.json','.json.json'))]
    (args.record/'extension_receipts.json').write_text(json.dumps(receipts,indent=1)+'\n')
    # Freeze dynamic catalogue and coverage responses including their timestamps.
    # FITS URLs have independent hashes of the measured arrays in the visit records.
    text_inputs={}
    for rec in receipts:
        original=args.cache/rec['file'][:-5]
        if original.is_file() and original.suffix in {'.tsv','.json'} and rec.get('status')=='ok':
            if sha(original)!=rec['sha256']:
                raise ValueError('Catalogue or sector receipt mismatch: '+original.name)
            text_inputs[original.name]=original.read_bytes().decode('utf-8')
    (args.record/'extension_catalogue_responses.json.gz').write_bytes(
        gzip.compress(json.dumps(text_inputs,ensure_ascii=False).encode(),mtime=0))
    print(Counter(r['status'] for r in results),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('operation',choices=['prepare','acquire'])
    parser.add_argument('--parent',type=Path)
    parser.add_argument('--matches',type=Path)
    parser.add_argument('--primary',type=Path)
    parser.add_argument('--record',type=Path,required=True)
    parser.add_argument('--cache',type=Path)
    parser.add_argument('--curves',type=Path)
    parser.add_argument('--workers',type=int,default=5)
    parser.add_argument('--pixel-tools',type=Path)
    args=parser.parse_args()
    args.record.mkdir(parents=True,exist_ok=True)
    if args.operation=='prepare':
        prepare(args)
    else:
        acquire(args)
