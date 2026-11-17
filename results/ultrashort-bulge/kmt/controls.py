"""Check earlier KMT event catalogs, the KMT planet archive, and ancillary cone services."""
import argparse
import concurrent.futures
import csv
import json
import math
import re
from pathlib import Path
from urllib.parse import urlencode
from lxml import html

from freeze import HERE, BASE, get
from neighbours import angle

LEADS={'KMT-2024-BLG-'+i for i in ('0456','0816','1157','0792','3237','3429','0909','3326')} | {
    'KMT-2025-BLG-'+i for i in ('0300','1213','1446','1507','1760','1836','2510','3173','2093')}


def fetch(args):
    key,url,cache=args
    path=cache/(key+'.raw')
    rec=cache/(key+'.json')
    if not path.exists():
        try:
            raw,r=get(url)
            path.write_bytes(raw)
        except OSError as exc:
            r=dict(url=url,status=getattr(exc,'code',None),error=str(exc))
        rec.write_text(json.dumps(r)+'\n')
    return key


def coord_dist(a,d,b,e):
    return math.hypot((a-b)*math.cos(math.radians(d)), d-e)*3600


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args()
    cache=args.cache/'controls';cache.mkdir(exist_ok=True)
    leads=[r for r in map(json.loads,(args.cache/'index.jsonl').read_text().splitlines()) if r['event'] in LEADS]
    urls={f'kmt-{y}':f'{BASE}{y}/' for y in range(2016,2024)}
    urls['exoplanets']='https://kmtnet.kasi.re.kr/exoKMTNet/assets/data/planets.js?v=20261105b'
    for r in leads:
        name=r['event']
        ra,dec=angle(r['ra'],True),angle(r['dec'])
        adql=("SELECT main_id,otype,ra,dec FROM basic WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
              f"CIRCLE('ICRS',{ra:.9f},{dec:.9f},{2/3600:.10f}))")
        urls[name+'-simbad']='https://simbad.cfa.harvard.edu/simbad/sim-tap/sync?'+urlencode(
            dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=adql))
        urls[name+'-vsx']='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urlencode(
            {'-source':'B/vsx/vsx','-c':f'{ra:.9f} {dec:.9f}','-c.rs':'2','-out.max':'100','-out':'Name,Type,RAJ2000,DEJ2000'})
        if name in {'KMT-2024-BLG-0456','KMT-2025-BLG-1507','KMT-2025-BLG-1213','KMT-2024-BLG-3326','KMT-2024-BLG-3429'}:
            # A 5-arcmin cone dwarfs the time-scale offset between HJD and geocentric JD.
            # SkyBoT predicts known objects only; a negative query cannot veto all moving objects.
            urls[name+'-skybot']='https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?'+urlencode(
                {'RA':f'{ra:.9f}','DEC':f'{dec:.9f}','SR':'0.083333','TIME':str(r['t0']+2400000),
                 '-loc':'500','-mime':'text','-output':'object'})
    with concurrent.futures.ThreadPoolExecutor(5) as run:
        list(run.map(fetch,[(key,url,cache) for key,url in urls.items()]))
    records={k:json.loads((cache/(k+'.json')).read_text()) for k in urls}
    years={}
    historical=[]
    for y in range(2016,2026):
        pp=cache/f'kmt-{y}.raw' if y<2024 else args.cache/f'{y}.html'
        if not pp.exists(): continue
        d=html.fromstring(pp.read_bytes())
        rows=d.xpath('//tr')
        if not rows:continue
        cols=[c.text_content().strip() for c in rows[0].xpath('td')]
        if not all(x in cols for x in ['Event','RA','Dec']):continue
        nh=0
        for tr in rows[1:]:
            row=dict(zip(cols,[c.text_content().strip() for c in tr.xpath('td')]))
            try:
                ra,dec=angle(row['RA'],True),angle(row['Dec']);nh+=1
            except (KeyError,ValueError):continue
            for l in leads:
                if l['event']==row['Event']:continue
                s=coord_dist(ra,dec,angle(l['ra'],True),angle(l['dec']))
                if s<2:
                    historical.append(dict(event=l['event'],match=row['Event'],sep_arcsec=round(s,4)))
        years[y]=nh
    catalog=(cache/'exoplanets.raw').read_text()
    planets=json.loads(re.search(r'^window.KMTNET_PLANETS = (.*);$',catalog,re.M).group(1))
    ffp=[dict(event=r['name'],title=r['title'],pub_year=r['pub_year'],ra=r['ra'],dec=r['dec'],ads=r['ads_link'])
         for r in planets if r.get('is_ffp')]
    cone={}
    for r in leads:
        name=r['event'];d={}
        for svc in ('simbad','vsx','skybot'):
            key=name+'-'+svc
            if key not in records:continue
            path=cache/(key+'.raw')
            if not path.exists():
                d[svc]=dict(status='query_failed')
                continue
            raw=path.read_text(errors='replace')
            if svc=='simbad' and raw.startswith('main_id,'):
                import io
                d[svc]=dict(matches=list(csv.DictReader(io.StringIO(raw))))
            elif svc=='vsx' and '#RESOURCE=yCat_1' in raw and '#INFO\tError=' not in raw:
                data=[s for s in raw.splitlines() if s.strip() and not s.startswith('#')]
                d[svc]=dict(matches=data[3:])
            elif svc=='skybot':
                d[svc]=dict(response=raw[:3000])
            else:
                d[svc]=dict(status='unrecognized_response')
        cone[name]=d
    (HERE/'control-checks.json').write_text(json.dumps(dict(predecessor_rows=years,other_event_within_2arcsec=historical,
        exoplanet_list_date='2026-11-05',kmt_exoplanet_rows=len(planets),archive_ffps=ffp,cones=cone),indent=2)+'\n')
    (HERE/'control-receipts.json').write_text(json.dumps(records,indent=2)+'\n')


if __name__=='__main__':
    main()
