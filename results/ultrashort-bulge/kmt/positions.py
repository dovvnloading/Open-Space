"""Crossmatch both released KMT seasons against 2016–26 KMT and 2025–26 OGLE lists."""
import argparse
import csv
import json
import re
from pathlib import Path
from lxml import html
import numpy as np
from scipy.spatial import cKDTree
from freeze import HERE
from neighbours import angle


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args()
    rows=[]
    for year in range(2016,2027):
        path=args.cache/'controls'/f'kmt-{year}.raw' if year<2024 else args.cache/f'{year}.html'
        if year==2026 and not path.exists():
            # The 2026 release is metadata-only here. Preserve its already frozen
            # position associations when replaying without that mutable index.
            continue
        doc=html.fromstring(path.read_bytes());trs=doc.xpath('//tr')
        hdr=[c.text_content().strip() for c in trs[0].xpath('./td')]
        for tr in trs[1:]:
            row=dict(zip(hdr,[c.text_content().strip() for c in tr.xpath('./td')]))
            try:
                ra,dec=angle(row['RA'],True),angle(row['Dec'])
            except (KeyError,ValueError):continue
            try:t=float(row.get('t_0',''))
            except (TypeError,ValueError):t=None
            rows.append((row['Event'],ra,dec,t))
    with (HERE.parent/'ogle/catalog.csv').open() as h:
        for r in csv.DictReader(h):
            rows.append(('OGLE-'+r['event'],angle(r['ra_hms'],True),angle(r['dec_dms']),float(r['t0_hjd'])-2400000.))
    a,d=np.deg2rad([(r[1],r[2]) for r in rows]).T
    unit=np.column_stack((np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)))
    pairs=cKDTree(unit).query_pairs(2*np.sin(np.deg2rad(2/3600)/2))
    targets={r['event'] for r in map(json.loads,(args.cache/'index.jsonl').read_text().splitlines())}
    records=[]
    for i,j in sorted(pairs):
        for k,l in [(i,j),(j,i)]:
            if rows[k][0] in targets:
                dist=np.rad2deg(2*np.arcsin(np.linalg.norm(unit[k]-unit[l])/2))*3600.
                # Do not reproduce current-season proprietary fit values in the match file.
                delta=None if 'KMT-2026' in rows[l][0] or rows[k][3] is None or rows[l][3] is None else round(rows[l][3]-rows[k][3],4)
                records.append(dict(event=rows[k][0],neighbor=rows[l][0],sep_arcsec=round(dist,4),delta_peak_days=delta))
    # Dedicated KMT planet list includes KMT and OGLE object names; use both exact
    # designation and 2-arcsec position to catch aliases. Published follow-up can precede a list update.
    text=(args.cache/'controls/exoplanets.raw').read_text()
    planets=json.loads(re.search(r'^window.KMTNET_PLANETS = (.*);$',text,re.M).group(1))
    pub=[]
    tree=cKDTree(unit)
    for planet in planets:
        ra,dec=planet.get('ra'),planet.get('dec')
        eligible=set()
        if ra is not None and dec is not None:
            a,d=np.deg2rad([ra,dec]);v=np.array([np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)])
            for idx in tree.query_ball_point(v,2*np.sin(np.deg2rad(2/3600)/2)):
                if rows[idx][0] in targets:eligible.add(rows[idx][0])
        exact=re.search(r'KMT-(?:2024|2025)-BLG-\d{4}',planet['name'])
        if exact and exact.group(0) in targets:eligible.add(exact.group(0))
        for event in sorted(eligible):
            pub.append(dict(event=event,archive_name=planet['name'],is_ffp=planet['is_ffp'],paper=planet['title'],year=planet['pub_year'],ads_link=planet['ads_link']))
    with (HERE/'planet-archive-matches.tsv').open('w') as h:
        w=csv.DictWriter(h,['event','archive_name','is_ffp','paper','year','ads_link'],delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(pub)
    if not (args.cache/'2026.html').exists() and (HERE/'historical-alert-neighbours.tsv').exists():
        with (HERE/'historical-alert-neighbours.tsv').open() as h:
            records.extend(r for r in csv.DictReader(h,delimiter='\t') if r['neighbor'].startswith('KMT-2026-'))
    with (HERE/'historical-alert-neighbours.tsv').open('w') as h:
        w=csv.DictWriter(h,['event','neighbor','sep_arcsec','delta_peak_days'],delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(sorted(records,key=lambda r:(r['event'],float(r['sep_arcsec']))))


if __name__=='__main__':main()
