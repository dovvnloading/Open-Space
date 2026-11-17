"""Build same-season cross-survey ID groups with provenance for every edge.

Only explicit survey cross-identifications join groups. Coordinate-only
neighbors remain in other/matches.csv and are never treated as guaranteed IDs.
"""
import argparse
import collections
import csv
import re
from pathlib import Path
from consolidate import inputs, rows, write_gz, HERE


def canonical(name):
    m=re.fullmatch(r'(MOA|PRIME|OGLE|KMT)-(\d{4})-BLG-(\d+)',name,re.I)
    if m:
        n,yr,i=m.groups()
        width=3 if n.upper() in ('MOA','PRIME') and yr=='2024' else 4
        return f'{n.upper()}-{yr}-BLG-{int(i):0{width}d}'
    return name


def associations():
    spec=inputs();edges=set()
    for path in (HERE/'kmt').glob('catalog-*.tsv'):
        for r in rows(path,'\t'):
            for val in r['Related event'].split():
                m=re.fullmatch('([OMP])B(\\d{2})(\\d+)',val)
                if not m:continue
                survey,yr,number=m.groups()
                other=canonical({'O':'OGLE','M':'MOA','P':'PRIME'}[survey]+'-20'+yr+'-BLG-'+number)
                edges.add((r['Event'],other,'KMT Related event'))
    for r in rows(HERE/'other/official_links.csv'):
        edges.add((canonical(r['id_1']),canonical(r['id_2']),'MOA/PRIME '+r['pair']))
    for r in rows(HERE/'other/known_matches.csv'):
        if r['status']=='known by published ID' and r['known_id']!=r['alert_id']:
            edges.add((canonical(r['known_id']),canonical(r['alert_id']),'published identity'))
    parent={n:n for n in spec}
    def rep(n):
        while n!=parent[n]:
            parent[n]=parent[parent[n]];n=parent[n]
        return n
    edge_rows=[]
    for a,b,source in sorted(edges):
        inside=a in spec and b in spec and spec[a]['year']==spec[b]['year']
        if inside:
            u,v=rep(a),rep(b)
            parent[max(u,v)]=min(u,v)
        edge_rows.append({'a':a,'b':b,'source':source,'joined_inside_numeric_inventory':int(inside)})
    groups=collections.defaultdict(list)
    for n in spec:groups[rep(n)].append(n)
    gr=[{'group_id':k,'n_alert_ids':len(v),'alert_ids':'|'.join(sorted(v))} for k,v in sorted(groups.items())]
    write_gz(HERE/'tables/official-cross-identifications.csv.gz',edge_rows)
    write_gz(HERE/'tables/survey-id-groups.csv.gz',gr)
    print(len(spec),len(groups),sum(len(v)>1 for v in groups.values()),'official-ID groups only')


if __name__=='__main__':associations()
