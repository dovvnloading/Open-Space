"""Reproduce the same-band ZTF selection from the per-discovery-date replies."""
import collections
import csv
import gzip
import itertools
import json
import math
from pathlib import Path
from ztf_tap import META, RAW, ROOT, WINDOW

def write_csv(path, data, fields):
    with path.open('w') as file:
        writer = csv.DictWriter(file, fieldnames=fields, lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)

def num(row, field):
    try:
        return float(row[field])
    except (ValueError, TypeError):
        return math.nan

def main():
    objects = collections.defaultdict(list)
    counts = collections.Counter()
    max_magerr = 2.5 / math.log(10) / 5
    for date in range(*WINDOW):
        with gzip.open(META / f'{date}_count.csv.gz','rt') as count_file:
            expected = int(next(csv.DictReader(count_file))['nrows'])
        with gzip.open(RAW / f'{date}_detections.csv.gz','rt') as file:
            rows = list(csv.DictReader(file))
        if len(rows) != expected:
            raise ValueError(f'Incomplete day {date}: {len(rows)}/{expected}')
        counts['raw_detection_rows'] += len(rows)
        for r in rows:
            objects[r['oid']].append(r)
    pairs = []
    for oid, raw in sorted(objects.items()):
        # Distinct candid values preserve exact 64-bit identities.
        points = {r['candid']: r for r in raw}.values()
        counts['unique_candidates'] += len(points)
        good = [r for r in points if num(r,'isdiffpos')==1 and 0 < num(r,'sigmapsf') <= max_magerr and math.isfinite(num(r,'magpsf'))]
        counts['positive_snr5_detections'] += len(good)
        if len(good)>=2:
            counts['objects_two_positive_snr5'] += 1
        baseline_seen = False
        for a,b in itertools.combinations(sorted(good,key=lambda r:num(r,'mjd')),2):
            dt = num(b,'mjd') - num(a,'mjd')
            if a['fid']!=b['fid'] or not (0.02<=dt<=1) or a['pid']==b['pid']:
                continue
            baseline_seen=True
            counts['eligible_pairs']+=1
            dm=num(b,'magpsf')-num(a,'magpsf')
            sigma=math.hypot(num(a,'sigmapsf'),num(b,'sigmapsf'))
            if abs(dm)>=0.30 and abs(dm)>dt and abs(dm)>=3*sigma:
                pairs.append(dict(oid=oid,fid=a['fid'],candid1=a['candid'],candid2=b['candid'],
                     mjd1=a['mjd'],mjd2=b['mjd'],mag1=a['magpsf'],mag2=b['magpsf'],
                     err1=a['sigmapsf'],err2=b['sigmapsf'],dt=dt,dm=dm,
                     mag_per_day=dm/dt,sigma=abs(dm)/sigma,meanra=a['meanra'],meandec=a['meandec'],
                     firstmjd=a['firstmjd'],lastmjd=a['lastmjd'],mjdstarthist=a['mjdstarthist'],
                     drb1=a['drb'],drb2=b['drb'],stellar=a['stellar']))
        if baseline_seen:
            counts['objects_with_eligible_pairs'] += 1
    counts['objects'] = len(objects)
    counts['selected_pairs'] = len(pairs)
    counts['selected_objects'] = len({p['oid'] for p in pairs})
    output=ROOT/'derived'
    output.mkdir(exist_ok=True)
    fields='oid fid candid1 candid2 mjd1 mjd2 mag1 mag2 err1 err2 dt dm mag_per_day sigma meanra meandec firstmjd lastmjd mjdstarthist drb1 drb2 stellar'.split()
    write_csv(output/'ztf_pairs.csv',pairs,fields)
    (output/'ztf_counts.json').write_text(json.dumps(dict(sorted(counts.items())),indent=2)+'\n')
    print(dict(counts))
if __name__=='__main__':
    main()
