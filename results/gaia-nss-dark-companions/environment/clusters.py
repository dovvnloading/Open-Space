"""Crossmatch both mass screens against Hunt & Reffert (2024) member assignments."""
import argparse
import csv
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
URL = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/A%2bA/686/A42/members&-out=GaiaDR3,Name,Prob&-out.max=2000000"


def main(cache):
    if not cache.exists():
        with urllib.request.urlopen(URL,timeout=240) as r,cache.open('wb') as w:
            shutil.copyfileobj(r,w)
    selected = {}
    for sample in ('gaia','extension'):
        with (HERE.parent/sample/'mass_screen.csv').open() as f:
            for row in csv.DictReader(f):
                sid=row['source_id']
                selected[sid]=max(float(row['m2_screen_msun']),selected.get(sid,0))
    matches, count = [], 0
    with cache.open() as f:
        for line in f:
            items=line.rstrip().split('\t')
            items[0]=items[0].strip()
            if not items[0].isdecimal(): continue
            sid, name, prob = items
            count+=1
            if sid in selected:
                matches.append(dict(source_id=sid,cluster_name=name.strip(),member_probability=float(prob),
                                    max_screen_msun=round(selected[sid],5)))
    if count!=1291929: raise ValueError(f"Incomplete member table: {count}")
    with (HERE/'cluster_matches.csv').open('w') as f:
        w=csv.DictWriter(f,matches[0],lineterminator='\n');w.writeheader();w.writerows(sorted(matches,key=lambda r:-r['max_screen_msun']))
    with cache.open('rb') as f: sha=hashlib.file_digest(f,'sha256').hexdigest()
    (HERE/'depth.json').write_text(json.dumps(dict(url=URL,sha256=sha,input_rows=count,
        selected_unique=len(selected),matched_unique=len({r['source_id'] for r in matches}),
        above_3_matched=len({r['source_id'] for r in matches if r['max_screen_msun']>3}),
        HD264291_catalog_member=any(r['source_id']=='3378588057203660160' for r in matches)),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path)
    main(p.parse_args().cache)
