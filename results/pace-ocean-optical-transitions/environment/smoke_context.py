"""Classify Florida OCI sites against all NOAA HMS smoke polygons on observation days."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parent
BASE='https://satepsanone.nesdis.noaa.gov/pub/FIRE/web/HMS/Smoke_Polygons/Shapefile/'
DAYS=('2025-05-13','2025-05-17','2025-05-21')

def in_ring(x,y,ring):
    inside=False
    for (x1,y1),(x2,y2) in zip(ring,ring[1:]):
        if (y1>y)!=(y2>y) and x < (x2-x1)*(y-y1)/(y2-y1)+x1:
            inside=not inside
    return inside

def contains(x,y,g):
    parts=g['coordinates'] if g['type']=='MultiPolygon' else [g['coordinates']]
    return any(in_ring(x,y,rings[0]) and not any(in_ring(x,y,r) for r in rings[1:])
               for rings in parts)

def main():
    with open(sys.argv[1]) as f:
        optical=[r for r in csv.DictReader(f) if r['region']=='florida_atlantic']
    rows=[]
    receipts=[]
    with tempfile.TemporaryDirectory() as tmp:
        for day in DAYS:
            short=day.replace('-','')
            url=BASE+f'{short[:4]}/{short[4:6]}/hms_smoke{short}.zip'
            with urlopen(url,timeout=70) as r: raw=r.read()
            path=Path(tmp)/f'hms{short}.zip'
            path.write_bytes(raw)
            geo=json.loads(subprocess.check_output(['ogr2ogr','-f','GeoJSON','/vsistdout/',f'/vsizip/{path}']))
            receipts.append(dict(date=day,url=url,bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest(),polygons=len(geo['features'])))
            for p in optical:
                if day not in (p['start'],p['end']): continue
                x,y=float(p['lon']),float(p['lat'])
                overlaps=[s['properties'] for s in geo['features'] if contains(x,y,s['geometry'])]
                rows.append(dict(start=p['start'],end=p['end'],date=day,lat=y,lon=x,
                                 polygon_count=len(overlaps),classes='|'.join(sorted({s['Density'] for s in overlaps})),
                                 polygon_times_UTC='|'.join(s['Start']+'-'+s['End'] for s in overlaps)))
    with (ROOT/'smoke_footprint.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    (ROOT/'smoke_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')

if __name__=='__main__': main()
