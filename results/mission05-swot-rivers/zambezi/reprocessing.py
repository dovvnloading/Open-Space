"""Retrieve or replay same-pass RiverSP v2 vs D comparisons at two station reaches.

python reprocessing.py retrieve DATA_DIRECTORY D_DIRECTORY
python reprocessing.py compare DATA_DIRECTORY D_DIRECTORY
"""
import csv
import io
import json
import math
import statistics
import struct
import sys
import zlib
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from acquire import HYDROCRON, FIELDS, fetch, sha
from coverage import geometry_ok, value

TRANSLATION = 'https://zenodo.org/records/22259077/files/SWORD_v17b_to_v16_Translation.zip'
SITES = {'12279000221': 'nanas_farm', '12297000161': 'chavuma'}
# Central directory of the DOI-pinned ZIP, v17b/AF_ReachIDs_v17b_vs_v16.csv.
OFFSET, COMPRESSED_SIZE, CRC = 94336333, 521785, 1043955410
MEMBER_SHA256 = '507f50d1e93a1f70e84519103001d4f1759c38a21bd25b9f59f848ec708a7211'


def retrieve(out):
    out.mkdir(parents=True, exist_ok=True)
    end = OFFSET + COMPRESSED_SIZE + 255
    with urlopen(Request(TRANSLATION, headers={'Range': f'bytes={OFFSET}-{end}'}), timeout=60) as response:
        assert response.getheader('Content-Range', '').startswith(f'bytes {OFFSET}-')
        data = response.read()
    header = struct.unpack_from('<4s5H3L2H', data)
    assert header[0] == b'PK\x03\x04'
    skip = 30 + sum(header[-2:])
    uncompressed = zlib.decompress(data[skip:skip + COMPRESSED_SIZE], -15)
    assert zlib.crc32(uncompressed) == CRC and sha(uncompressed) == MEMBER_SHA256
    translations = [r for r in csv.DictReader(io.StringIO(uncompressed.decode())) if r['v16_reach_id'] in SITES]
    assert sorted(r['v16_reach_id'] for r in translations) == sorted(SITES)
    metadata = dict(translation_url=TRANSLATION, member='v17b/AF_ReachIDs_v17b_vs_v16.csv',
                    member_sha256=MEMBER_SHA256, translation=translations, hydrocron=[])
    for tr in translations:
        v17 = tr['v17_reach_id']
        query = dict(feature='Reach', feature_id=v17, start_time='2023-03-29T00:00:00Z',
                     end_time='2025-05-03T23:59:59Z', output='csv', collection_name='SWOT_L2_HR_RiverSP_D',
                     fields=FIELDS.split(',sos_')[0])
        url = HYDROCRON + '?' + urlencode(query)
        raw = fetch(url)
        payload = json.loads(raw)['results']['csv'].encode()
        name = 'vD_' + v17 + '.csv'
        (out / name).write_bytes(payload)
        metadata['hydrocron'].append(dict(file=name, url=url, v16=tr['v16_reach_id'], v17b=v17,
                                          csv_sha256=sha(payload), response_sha256=sha(raw)))
    (out / 'vD_manifest.json').write_text(json.dumps(metadata, indent=2)+'\n')


def loadrows(path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    # For an identical observation appearing in two deliveries, keep the greatest
    # numerical granule revision. Distinct observations in a cycle/pass are ambiguous.
    indexed = {}
    for row in rows:
        key = row['cycle_id'], row['pass_id']
        if key in indexed:
            old = indexed[key]
            if (old['time_str'], old['crid']) != (row['time_str'], row['crid']):
                raise ValueError(f'ambiguous cycle/pass in {path}: {key}')
            revision = lambda r: int(r['granuleUR'].split('_')[-1].split('.')[0])
            if revision(old) >= revision(row):
                continue
        indexed[key] = row
    return indexed


def write(path, rows):
    with path.open('w') as f:
        writer = csv.DictWriter(f, fieldnames=rows[0], lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def compare(data, out):
    manifest = json.loads((out / 'vD_manifest.json').read_text())
    comparisons, summary = [], []
    for src in manifest['hydrocron']:
        a = loadrows(data / f"reach_{src['v16']}.csv")
        b = loadrows(out / src['file'])
        keys = sorted(set(a) & set(b))
        for cycle, pid in keys:
            x, y = a[cycle, pid], b[cycle, pid]
            rec = dict(site=SITES[src['v16']], v16=src['v16'], v17b=src['v17b'], cycle=cycle, pass_id=pid,
                       date_v2=x['time_str'], date_d=y['time_str'], qa_v2=int(geometry_ok(x)), qa_d=int(geometry_ok(y)),
                       reach_q_v2=x['reach_q'], reach_q_d=y['reach_q'], crid_v2=x['crid'], crid_d=y['crid'])
            for field in ('p_lat', 'p_lon', 'p_length', 'p_n_nodes', 'wse', 'width', 'slope', 'loc_offset',
                          'node_dist', 'dark_frac', 'wse_u', 'width_u', 'slope_u'):
                for suffix, row in [('v2', x), ('d', y)]:
                    val = value(row, field)
                    rec[f'{field}_{suffix}'] = val if math.isfinite(val) else ''
            comparisons.append(rec)
        for pid in sorted({key[1] for key in keys}):
            for mode in ('v2_good', 'both_good'):
                common = [r for r in comparisons if r['v16'] == src['v16'] and r['pass_id'] == pid and r['qa_v2']
                          and (mode != 'both_good' or r['qa_d'])]
                rec = dict(site=SITES[src['v16']], pass_id=pid, subset=mode, pairs=len(common))
                for measure in ('wse', 'width', 'slope'):
                    diffs = [r[measure+'_d']-r[measure+'_v2'] for r in common if '' not in (r[measure+'_d'],r[measure+'_v2'])]
                    rec['n_'+measure] = len(diffs)
                    rec['median_d_minus_v2_'+measure] = statistics.median(diffs) if diffs else ''
                summary.append(rec)
    write(out / 'same_pass.csv', comparisons)
    write(out / 'same_pass_summary.csv', summary)


if __name__ == '__main__':
    mode, data, out = sys.argv[1:]
    if mode == 'retrieve':
        retrieve(Path(out))
    elif mode != 'compare':
        raise ValueError(mode)
    compare(Path(data), Path(out))
