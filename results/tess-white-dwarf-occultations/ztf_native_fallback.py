"""Replay the four Gaia-target ZTF DR24 HATS extractions and their native shutter joins."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen

import numpy as np
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

TICS = ('357686299', '446165080', '408000326', '55997020')
COLUMNS = ['oid','mjd','exptime','mag','magerr','catflags','filtercode','ra','dec','coordinate_source']
PALOMAR = EarthLocation.from_geodetic(-116.863, 33.356, 1706)
iers.conf.auto_download = False
iers.conf.auto_max_age = None


def retrieve(url, sha, cache, byte_range=None):
    path = cache / sha
    if not path.exists():
        headers = {} if byte_range is None else {'Range': f'bytes={byte_range[0]}-{byte_range[1]}'}
        with urlopen(Request(url,headers=headers), timeout=120) as response:
            if byte_range is not None and response.status != 206:
                raise ValueError('Server did not honor byte range')
            raw = response.read()
        if hashlib.sha256(raw).hexdigest() != sha:
            raise ValueError('Archive response hash changed: '+url)
        path.write_bytes(raw)
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != sha:
        raise ValueError('Cached input hash mismatch')
    return data


def csv_rows(raw):
    return list(csv.DictReader(io.StringIO(raw.decode())))


def hats_rows(meta, cache):
    ids = pa.array([int(r['oid']) for r in csv_rows(retrieve(meta['table_query'],meta['table_sha256'],cache))])
    block = meta['hats_partition']
    # Sparse local Parquet reconstruction using only the original requested ranges.
    path = cache / ('partition-'+str(block['url'].split('Npix=')[-1].split('/')[0])+'.parquet')
    def get_range(entry):
        start,end,sha = entry
        return start, retrieve(block['url'],sha,cache,(start,end))
    with path.open('wb') as stream, ThreadPoolExecutor(3) as pool:
        stream.truncate(block['bytes'])
        for start,raw in pool.map(get_range,block['ranges']):
            stream.seek(start)
            stream.write(raw)
    reader = pq.ParquetFile(path)
    records = []
    for group in block['groups']:
        obj = reader.read_row_group(group, columns=['objectid'])
        if not pc.any(pc.is_in(obj['objectid'], value_set=ids)).as_py():
            continue
        group_data = reader.read_row_group(group, columns=['objectid','objra','objdec','lightcurve'])
        records.extend(group_data.filter(pc.is_in(group_data['objectid'],value_set=ids)).to_pylist())
    return sorted(records,key=lambda r:r['objectid'])


def image_shutters(frames, star):
    rows = csv_rows(frames)
    start = np.asarray([float(r['obsjd'])-2400000.5 for r in rows])
    exposure = np.asarray([float(r['exptime']) for r in rows])
    mid = Time(start+exposure/172800,format='mjd',scale='utc',location=PALOMAR)
    direction = SkyCoord(star['objra'],star['objdec'],unit='deg')
    hmjd = mid.utc.mjd + mid.light_travel_time(direction,kind='heliocentric').to_value('day')
    return start, exposure, hmjd


def build(meta, cache, out):
    stars = hats_rows(meta,cache)
    observations = []
    unmatched = []
    for star in stars:
        oid = str(star['objectid'])
        series = star['lightcurve']
        inputs = meta['terminal_inputs']
        by_time = {}
        if oid in inputs:
            src = inputs[oid]
            for row in csv_rows(retrieve(src['url'],src['sha256'],cache)):
                by_time.setdefault(round(float(row['hjd'])-2400000.5,5),[]).append(row)
        frames = None
        band = ['','zg','zr','zi'][int(oid[-12])]
        image_key = f'{int(oid[:-12])}-{int(oid[-11:-9])}-{int(oid[-9])}-{band}'
        for i, epoch in enumerate(series['hmjd']):
            match = [r for r in by_time.get(round(epoch,5),[])
                     if abs(float(r['mag'])-series['mag'][i])<.00002
                     and abs(float(r['magerr'])-series['magerr'][i])<.00002
                     and int(r['catflags'])==series['catflags'][i]]
            if len(match)==1:
                record = {k:match[0][k] for k in ('mjd','exptime','filtercode','ra','dec')}
                record['coordinate_source'] = 'IRSA epoch'
            else:
                if frames is None:
                    src = inputs[image_key]
                    raw = retrieve(src['url'],src['sha256'],cache)
                    if image_key in meta.get('supplemental_shutter_queries', {}):
                        extra = meta['supplemental_shutter_queries'][image_key]
                        other_chips = retrieve(extra['url'],extra['sha256'],cache)
                        raw = raw.rstrip()+b'\n'+b'\n'.join(other_chips.splitlines()[1:])+b'\n'
                    image_rows = csv_rows(raw)
                    frames = image_shutters(raw,star)
                starts, exposures, predicted = frames
                compatible = np.flatnonzero(abs(predicted-epoch)*86400 <= 1.)
                distinct = {(image_rows[j]['obsjd'],image_rows[j]['exptime']) for j in compatible}
                if len(distinct)!=1:
                    unmatched.append([oid,epoch,sorted(distinct)])
                    continue
                j = compatible[0]
                record = dict(mjd=f'{starts[j]:.10f}',exptime=image_rows[j]['exptime'],filtercode=band,
                              ra=star['objra'],dec=star['objdec'],coordinate_source='HATS reference')
            observations.append(dict(oid=oid,mag=series['mag'][i],magerr=series['magerr'][i],
                                     catflags=series['catflags'][i],**record))
    buf = io.StringIO()
    writer = csv.DictWriter(buf,COLUMNS,lineterminator='\n')
    writer.writeheader()
    writer.writerows(observations)
    body = buf.getvalue().encode()
    if unmatched!=meta['unmatched_shutter_hmjd'] or hashlib.sha256(body).hexdigest()!=meta['sha256']:
        raise ValueError('Reconstructed epoch table disagrees with frozen output')
    filename = out/meta['cache_name']
    filename.write_bytes(gzip.compress(body,mtime=0))
    filename.with_suffix('.receipt').write_text(json.dumps(meta,indent=1)+'\n')
    print(meta['cache_name'], len(observations), flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipts',type=Path,help='catalogue_checks.json.gz or the directory of four native .receipt files')
    p.add_argument('cache',type=Path)
    p.add_argument('out',type=Path)
    a=p.parse_args()
    for part in (a.cache,a.out):
        part.mkdir(parents=True,exist_ok=True)
    if a.receipts.is_dir():
        metas = {tic:json.loads((a.receipts/(tic+'-ztf-native.receipt')).read_text()) for tic in TICS}
    else:
        with gzip.open(a.receipts,'rt') as f:
            metas = {r['tic']:r['queries']['ztf'] for r in json.load(f) if r['tic'] in TICS}
    for tic in TICS:
        build(metas[tic],a.cache,a.out)
