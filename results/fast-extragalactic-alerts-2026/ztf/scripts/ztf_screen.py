"""Reproduce the same-band ZTF selection from the per-discovery-date replies."""
import collections
import argparse
import hashlib
import io
import tarfile
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

def main(output):
    # The supplement contains the original count bodies and receipts. Verify
    # every response before parsing measurements or touching outputs.
    archive = ROOT / 'ztf_evidence_supplement.tar.gz'
    with tarfile.open(archive, 'r:gz') as tar:
        for day in range(*WINDOW):
            for tag in (f'{day}_count', f'{day}_detections'):
                rec = json.load(tar.extractfile(f'receipts/ztf/{tag}.json'))
                if tag.endswith('_count'):
                    raw = gzip.decompress(tar.extractfile(f'receipts/ztf/{tag}.csv.gz').read())
                else:
                    raw = gzip.decompress((RAW / (tag + '.csv.gz')).read_bytes())
                if rec.get('http_status') != 200 or len(raw) != rec.get('bytes') or hashlib.sha256(raw).hexdigest() != rec.get('sha256'):
                    raise ValueError('Frozen input mismatch: ' + tag)
    if output.resolve().is_relative_to(ROOT.resolve()):
        raise ValueError('Output must be outside the frozen result')
    if output.exists():
        raise FileExistsError(output)
    objects = collections.defaultdict(list)
    counts = collections.Counter()
    max_magerr = 2.5 / math.log(10) / 5
    for date in range(*WINDOW):
        with tarfile.open(archive, 'r:gz') as tar:
            raw_count = gzip.decompress(tar.extractfile(f'receipts/ztf/{date}_count.csv.gz').read()).decode()
        expected = int(next(csv.DictReader(io.StringIO(raw_count)))['nrows'])
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
    output.mkdir(parents=True, exist_ok=False)
    fields='oid fid candid1 candid2 mjd1 mjd2 mag1 mag2 err1 err2 dt dm mag_per_day sigma meanra meandec firstmjd lastmjd mjdstarthist drb1 drb2 stellar'.split()
    write_csv(output/'ztf_pairs.csv',pairs,fields)
    (output/'ztf_counts.json').write_text(json.dumps(dict(sorted(counts.items())),indent=2)+'\n')
    print(dict(counts))
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('output', type=Path, help='new output directory')
    main(parser.parse_args().output)
