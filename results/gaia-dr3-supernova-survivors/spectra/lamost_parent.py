"""Exact-ID crossmatch of every LAMOST DR12 v2 low-resolution WD-table entry to Gaia DR3."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import time
import urllib.parse
import urllib.request

URL = 'https://www.lamost.org/dr12/v2.0/catdl?name=dr12_v2.0_LRS_wd.csv.gz'
SHA = '5ee713b05a86533f528a9c29086387ff5aee9d89ab90543edafa3a49540d425a'
TAP = 'https://gea.esac.esa.int/tap-server/tap/sync'
FIELDS = ('source_id,ra,dec,parallax,parallax_error,pm,pmra,pmdec,ruwe,'
          'visibility_periods_used,duplicated_source,ipd_frac_multi_peak,'
          'ipd_gof_harmonic_amplitude,astrometric_params_solved,phot_g_mean_mag,bp_rp')
SCALE = 4.740470463


def sha(data):
    return hashlib.sha256(data).hexdigest()


def fetch(chunk, cache):
    query = f"SELECT {FIELDS} FROM gaiadr3.gaia_source WHERE source_id IN ({','.join(chunk)}) ORDER BY source_id"
    path = cache / (sha(query.encode()) + '.csv')
    if not path.exists():
        data = urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',
                                           MAXREC=str(len(chunk)+1),QUERY=query)).encode()
        for trial in range(3):
            try:
                with urllib.request.urlopen(TAP, data=data, timeout=120) as response:
                    raw = response.read()
                if not raw.startswith(b'source_id,'):
                    raise ValueError(f'Archive did not return CSV: {raw[:150]}')
                path.write_bytes(raw)
                break
            except Exception:
                if trial == 2:
                    raise
                time.sleep(1 + trial)
    raw = path.read_bytes()
    rows = list(csv.DictReader(io.StringIO(raw.decode())))
    found = [r['source_id'] for r in rows]
    if set(found) - set(chunk) or len(set(found)) != len(found):
        raise ValueError('Unexpected extra or duplicate Gaia records')
    return rows, dict(query_sha256=sha(query.encode()),response_sha256=sha(raw), bytes=len(raw),
                      requested=len(chunk), returned=len(rows),missing=sorted(set(chunk)-set(found)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=Path('.cache'))
    parser.add_argument('--output', type=Path, default=Path('.'))
    args = parser.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    cat = args.cache/'lamost_wd.csv.gz'
    if not cat.exists():
        with urllib.request.urlopen(URL,timeout=120) as response:
            cat.write_bytes(response.read())
    raw=cat.read_bytes()
    if sha(raw)!=SHA:
        raise ValueError('LAMOST archive changed. Inspect the new release before using it.')
    obs = list(csv.DictReader(io.StringIO(gzip.decompress(raw).decode())))
    ids = sorted({r['gaia_source_id'] for r in obs if int(r['gaia_source_id'])>0},key=int)
    chunks=[ids[i:i+350] for i in range(0,len(ids),350)]
    receipts=[]
    rows=[]
    with ThreadPoolExecutor(4) as pool:
        for batch,receipt in pool.map(lambda c:fetch(c,args.cache),chunks):
            rows.extend(batch)
            receipts.append(receipt)
    counts = Counter()
    high=[]
    eligible=[]
    upper=[]
    for row in rows:
        if not all(row[k] for k in ('pm','parallax','parallax_error','ruwe')):
            continue
        if float(row['parallax_error'])<=0:
            continue
        counts['astrometric_solution']+=1
        p,e,pm,ruwe=map(float,(row['parallax'],row['parallax_error'],row['pm'],row['ruwe']))
        if p<=0:
            continue
        counts['positive_parallax']+=1
        sn=p/e
        vt=SCALE*pm/p
        if vt>=500 and sn>=3:
            high.append(dict(**row,snr=sn,raw_vtan=vt))
        if sn<5 or ruwe>=1.4:
            continue
        eligible.append((vt,row['source_id']))
        upper.append((SCALE*pm/(p-1.645*e),row['source_id']))
        if vt>=500: counts['raw_vtan_ge500']+=1
        if vt>=1000: counts['raw_vtan_ge1000']+=1
    counts['raw_vtan_ge1000'] += 0
    counts['raw_vtan_ge500'] += 0
    counts['snr_ge5_ruwe_lt1p4']=len(eligible)
    counts['lamost_records']=len(obs)
    counts['records_with_positive_gaia_id']=sum(int(r['gaia_source_id'])>0 for r in obs)
    counts['unique_positive_gaia_ids']=len(ids)
    counts['unique_matched']=len(rows)
    counts['snr_ge3_raw_vtan_ge500_including_bad_ruwe']=len(high)
    args.output.mkdir(exist_ok=True,parents=True)
    with (args.output/'lamost_high_velocity.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=FIELDS.split(',')+['snr','raw_vtan'])
        writer.writeheader()
        writer.writerows(sorted(high,key=lambda x:-x['raw_vtan']))
    selected={r['source_id'] for r in high}
    with (args.output/'lamost_high_visits.csv').open('w') as f:
        writer=csv.DictWriter(f,fieldnames=['gaia_source_id','obsid','obsdate','designation','planid','lmjd','spid','fiberid','snrg','wd_subclass','z','z_err'], extrasaction='ignore')
        writer.writeheader()
        writer.writerows(r for r in obs if r['gaia_source_id'] in selected)
    controls=list(csv.DictReader((Path(__file__).parent.parent/'control'/'controls.csv').open()))
    membership=set(ids)
    result=dict(control_count=len(controls),controls_in_lamost_WD_table=[r['name'] for r in controls if r['source_id'] in membership],
                maximum_vtan_at_parallax_minus_1p645_sigma_fixed_pm=max(upper),catalog=URL,archive_sha256=sha(raw),tap=TAP,chunk_size=350,
                source_id_order='sorted numerically ascending, split into consecutive chunks',
                counts=dict(counts),maximum_raw_vtan_in_snr5_ruwe14=sorted(eligible,reverse=True)[:3],
                requests=receipts)
    (args.output/'lamost_receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result['counts'],result['maximum_raw_vtan_in_snr5_ruwe14'])


if __name__=='__main__':
    main()
