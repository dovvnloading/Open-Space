"""Check the program of every ZTF exposure used in an injection-eligible pair."""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import io
import itertools
import json
import tarfile
import time
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

from simulate import ROOT, load

HERE = Path(__file__).resolve().parent
SUPPLEMENT = ROOT / 'ztf/ztf_evidence_supplement.tar.gz'
ARCHIVE = HERE / 'ztf_program_replies.zip'
TAP = 'https://tap.alerce.online/tap/sync'


def key(fid, mjd):
    return str(fid), round(float(mjd), 7)


def packet_programs(raw, programs):
    packet = json.loads(raw)
    for row in [packet['candidate'], *(packet['prv_candidates'] or []), *(packet.get('fp_hists') or [])]:
        if 'jd' in row and row.get('programid') is not None:
            programs[key(row['fid'], float(row['jd']) - 2400000.5)].add(str(row['programid']))


def add_forced(raw, programs):
    reader = csv.DictReader(io.StringIO(raw.decode()))
    if reader.fieldnames != ['mjd', 'fid', 'programid']:
        raise ValueError('Unexpected TAP result schema')
    for row in reader:
        programs[key(row['fid'], row['mjd'])].add(row['programid'])


def prior_evidence():
    programs = defaultdict(set)
    with tarfile.open(SUPPLEMENT) as source:
        for member in source:
            if member.name.startswith('data/ztf/public_visit_batch_') and member.name.endswith('.csv.gz'):
                add_forced(gzip.decompress(source.extractfile(member).read()), programs)
            elif '/program_packets/' in member.name and member.name.endswith('.json.gz'):
                packet_programs(gzip.decompress(source.extractfile(member).read()), programs)
    return programs


def eligible_visits():
    _, data = load('ztf')
    needed = set()
    for visits in data.values():
        bands = defaultdict(list)
        for (band, _), (mjd, _) in visits.items():
            bands[band].append(mjd)
        for band, times in bands.items():
            for a, b in itertools.combinations(sorted(times), 2):
                if 0.02 <= b - a <= 1:
                    needed.update((key(band, a), key(band, b)))
    return needed


def read_targets(missing):
    # Prefer direct packets. The earlier selection can supply historical rows.
    by_time = defaultdict(list)
    for path in sorted((ROOT / 'ztf/data/ztf').glob('*_detections.csv.gz')):
        with gzip.open(path, 'rt') as stream:
            for row in csv.DictReader(stream):
                visit = key(row['fid'], row['mjd'])
                if visit in missing:
                    by_time[visit].append(row)
    return {visit: min(rows, key=lambda r: (r['has_stamp'] != 'True', r['oid'], r['candid']))
            for visit, rows in sorted(by_time.items())}


def tap(query):
    return TAP + '?' + urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', MAXREC=200000, QUERY=query))


def fetch(item):
    name, url = item
    for attempt in range(3):
        try:
            with urlopen(url, timeout=50) as response:
                raw = response.read()
                receipt = dict(url=url, status=response.status, utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                               bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
            return name, (raw, receipt)
        except OSError:
            if attempt == 2:
                raise
            time.sleep(attempt + 1)


def batches(ids):
    values = sorted(set(ids))
    for start in range(0, len(values), 50):
        yield ','.join("'" + oid + "'" for oid in values[start:start + 50])


def refresh(needed, initial):
    replies = {}
    programs = defaultdict(set, {k: v.copy() for k, v in initial.items()})
    targets = read_targets(needed - programs.keys())
    with ThreadPoolExecutor(max_workers=3) as pool:
        requests = [(f'forced_{i}.csv', tap('SELECT DISTINCT mjd,fid,programid FROM ztf.forced_photometry '
                    f'WHERE oid IN ({batch}) AND mjd >= 61345 AND mjd < 61361'))
                    for i, batch in enumerate(batches(r['oid'] for r in targets.values()))]
        for name, (raw, receipt) in pool.map(fetch, requests):
            replies[name] = (raw, receipt)
            add_forced(raw, programs)
        # One parent packet can supply several missing history rows.
        remaining = [r for k, r in targets.items() if k not in programs]
        print('Still missing after forced-history lookups:', len(remaining), flush=True)
        parent = {}
        requests = [(f'parents_{i}.csv', tap(f'SELECT oid,candid,parent_candid FROM ztf.detection WHERE oid IN ({batch})'))
                    for i, batch in enumerate(batches(r['oid'] for r in remaining if r['has_stamp'] != 'True'))]
        for name, (raw, receipt) in pool.map(fetch, requests):
            replies[name] = (raw, receipt)
            for row in csv.DictReader(io.StringIO(raw.decode())):
                parent[row['candid']] = row['parent_candid']
        requests = {}
        for r in remaining:
            cid = r['candid'] if r['has_stamp'] == 'True' else parent[r['candid']]
            requests[f'packet_{cid}.json'] = 'https://avro.alerce.online/get_avro_info?' + urlencode(dict(oid=r['oid'], candid=cid))
        for name, reply in pool.map(fetch, sorted(requests.items())):
            replies[name] = reply
    # Keep exact response bodies alongside their query URLs, dates and hashes.
    with zipfile.ZipFile(ARCHIVE.with_suffix('.tmp'), 'w', compression=zipfile.ZIP_DEFLATED) as output:
        for name, (raw, _) in sorted(replies.items()):
            output.writestr(name, raw)
        output.writestr('receipts.json', json.dumps({n: r for n, (_, r) in replies.items()}, indent=2) + '\n')
    ARCHIVE.with_suffix('.tmp').replace(ARCHIVE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch', action='store_true', help='Replace the supplementary program receipts from live services')
    options = parser.parse_args()
    programs = prior_evidence()
    needed = eligible_visits()
    if options.fetch:
        refresh(needed, programs)
    with zipfile.ZipFile(ARCHIVE) as source:
        for name, meta in json.loads(source.read('receipts.json')).items():
            raw = source.read(name)
            if hashlib.sha256(raw).hexdigest() != meta['sha256']:
                raise ValueError('Changed reply: ' + name)
            if name.startswith('forced_'):
                add_forced(raw, programs)
            if name.startswith('packet_'):
                packet_programs(raw, programs)
    with (HERE / 'ztf_program_visits.csv').open('w') as out:
        writer = csv.writer(out, lineterminator='\n')
        writer.writerow(['fid', 'mjd', 'programids'])
        writer.writerows((fid, mjd, '|'.join(sorted(programs[fid, mjd]))) for fid, mjd in sorted(needed))
    print('Verified solely program 1:', sum(programs[k] == {'1'} for k in needed), 'of', len(needed))


if __name__ == '__main__':
    main()
