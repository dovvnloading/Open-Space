"""Check released optical alert catalogs, primary-paper roster and public image footprints."""
import argparse
import csv
import hashlib
import json
import math
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

from sightlines import EVENTS, HERE, KMT, CFA, CATS, catalog, fetch_catalog, urlencode
from neighbours import angle
from support import read
from counterparts import make_queries,run_one


def separation(ra,dec,aa,dd):
    a,d,b,e=map(math.radians,(ra,dec,aa,dd))
    s=math.sin((e-d)/2)**2+math.cos(d)*math.cos(e)*math.sin((b-a)/2)**2
    return math.degrees(2*math.asin(math.sqrt(min(1,s))))*3600


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--ogle24-par',type=Path,required=True);args=p.parse_args()
    all_catalogs=[]
    for name,file,fmt in [('OGLE_2024',args.ogle24_par,'par'),('OGLE_2025_26',KMT.parent/'ogle/catalog.csv','ogle'),
                          ('MOA_PRIME',KMT.parent/'other/index.csv','other'),('known_literature',KMT.parent/'other/known_roster.csv','lit')]:
        buf=file.read_bytes()
        source=[]
        if fmt=='par':
            for l in buf.decode().splitlines():
                vals=l.split()
                if vals and vals[0].startswith('2024-BLG-'):
                    source.append(('OGLE-'+vals[0],angle(vals[3],True),angle(vals[4])))
        else:
            for r in csv.DictReader(buf.decode().splitlines()):
                if fmt=='ogle':source.append(('OGLE-'+r['event'],angle(r['ra_hms'],True),angle(r['dec_dms'])))
                if fmt=='other':source.append((r['event'],angle(r['ra'],True),angle(r['dec'])))
                if fmt=='lit' and r['ra_deg'] and r['dec_deg']:source.append((r['known_id'],float(r['ra_deg']),float(r['dec_deg'])))
        all_catalogs.append((name,hashlib.sha256(buf).hexdigest(),source))
    catalog_output=[]
    jobs=[]
    for event in EVENTS+['KMT-2024-BLG-1762','KMT-2025-BLG-1797','KMT-2025-BLG-2489','KMT-2024-BLG-1152','KMT-2024-BLG-3413']:
        row=catalog(event)
        ra,dec=angle(row['RA'],True),angle(row['Dec'])
        matches=[]
        for name,sha,source in all_catalogs:
            dist=[(separation(ra,dec,a,d),ident) for ident,a,d in source]
            d0,n0=min(dist)
            matches.append(dict(catalog=name,sha256=sha,n_with_coordinates=len(dist),closest=n0,sep_arcsec=round(d0,3),
                                matches_3arcsec=[dict(id=n,sep_arcsec=round(d,3)) for d,n in sorted(dist) if d<=3]))
        catalog_output.append(dict(event=event,catalogs=matches))
        x=dict(event=event[4:],ra_hms=row['RA'],dec_dms=row['Dec'])
        for name,url in make_queries(x).items():
            if 'footprint' in name:
                jobs.append((args.cache,x,name,url))
    with ThreadPoolExecutor(4) as pool:images=list(pool.map(run_one,jobs))
    for x in images:x['event']='KMT-'+x['event']
    with (KMT.parent/'ogle/catalog.csv').open() as f:
        og=next(r for r in csv.DictReader(f) if r['event']=='2026-BLG-0670')
    with ThreadPoolExecutor(4) as pool:
        og_data=list(pool.map(run_one,[(args.cache,og,k,v) for k,v in make_queries(og).items()]))
    a,d=angle(og['ra_hms'],True),angle(og['dec_dms'])
    query=CFA+urlencode({'-source':CATS[0][1],'-c':f'{a:.9f} {d:.9f}','-c.rs':3,
                         '-out.add':'_r','-sort':'_r','-out.max':50,'-out':CATS[0][2]})
    v=fetch_catalog((args.cache,'0670virac',query))
    og_data.append(dict(event='2026-BLG-0670',product='virac_3arcsec',**v))
    (HERE/'OGLE-2026-BLG-0670.json').write_text(json.dumps(og_data,indent=2)+'\n')
    finders=[]
    for event in EVENTS[:2]:
        key=event.split('-')[-1]
        url=f'https://kmtnet.kasi.re.kr/ulens/event/2024/data/KB24{key}/fchart/fchart.jpg'
        _,r=read(args.cache/event,'finder',url)
        finders.append(dict(event=event,receipt=r))
    (HERE/'associations.json').write_text(json.dumps(dict(catalogs=catalog_output,image_footprints=images,
        kmt_finding_charts=finders),indent=2)+'\n')


if __name__=='__main__':main()
