"""Test actual ATLAS shutter intervals against all isolated TESS detection boxes."""
import argparse
import csv
import gzip
import hashlib
import json
import time
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timedelta, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
from astropy.utils import iers

from atlas_epoch_pixels import BASE, measure
from ffi_acquire import propagated

ROOT = Path(__file__).resolve().parent
# WGS84 east longitude, latitude, height (m), ATLAS Operations and Data Processing, §3.
SITES = {'atlas_haleakela': (-156.257,20.7075,3039), 'atlas_mauna_loa': (-155.5763,19.5362,3397),
         'atlas_rio_hurtado': (-70.76498,-30.47103,1609.6),
         'atlas_sutherland': (20.810317,-32.380461,1764)}
iers.conf.auto_download = False
iers.conf.auto_max_age = None


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_windows(scan):
    windows = []
    def add(tic, sector, sample, product, events):
        for i, e in enumerate(events):
            windows.append(dict(tic=str(tic), sector=int(sector), sample=sample, product=product,
                                box=i, epoch_btjd=e['epoch_btjd'], duration_min=e['duration_min']))
    for r in json.loads((ROOT/'extension/extension_measurements.json').read_text()):
        if r['status'] == 'included':
            add(r['tic'],r['sector'],'extension',r['lc_file'],r['search']['isolated'])
    for r in json.loads((ROOT/'record/historical_images.json').read_text()):
        events = [e for e in r['image_checks'] if e.get('diagnostic')=='individual trigger']
        if len(events) != r['n_isolated']:
            raise ValueError('Incomplete historical trigger reconstruction')
        add(r['tic'],r['sector'],'historical',f"TIC{r['tic']}-S{r['sector']}",events)
    if scan.is_dir():
        recs = [json.loads(f.read_text()) for f in sorted(scan.glob('*.json'))]
    else:
        with gzip.open(scan,'rt') as f:
            recs = json.load(f)
    for r in recs:
        add(r['tic'],r['sector'],'short_cadence',r['filename'],r.get('isolated',[]))
    r = json.loads((ROOT/'sector110/tic55997020/search.json').read_text())
    add(r['tic'],r['sector'],'sector110','9x9 direct FFI strips',r['search']['isolated'])
    return windows


def query_groups(windows):
    by_star = defaultdict(list)
    for idx,w in enumerate(windows):
        t = w['epoch_btjd']
        half = .5*w['duration_min']/1440.
        by_star[w['tic']].append((idx,t-half,t+half))
    for tic, boxes in sorted(by_star.items()):
        group = []
        for idx,left,right in sorted(boxes,key=lambda b:b[1]):
            if group and left > max(end for _,_,end in group)+2.:
                yield tic, group
                group = []
            group.append((idx,left,right))
        yield tic, group


def query(job):
    tic, boxes, stars, cache, previous = job
    # Pad both ends by at least one UTC day, comfortably exceeding annual light-time.
    low = min(left for _,left,_ in boxes)
    high = max(right for _,_,right in boxes)
    tt = Time([low,high],format='jd',scale='tdb',val2=2457000)
    dates = tt.utc.to_datetime()
    start, stop = dates[0].date()-timedelta(days=1), dates[1].date()+timedelta(days=2)
    ra,dec = propagated(stars[tic],float(np.mean(tt.jyear)))
    url = BASE+urlencode(dict(ra=round(ra,8), dec=round(dec,8), sources=','.join(SITES),
                              start_date=start.isoformat(), stop_date=stop.isoformat()))
    fn = cache/f'{tic}-{start}-{stop}'
    if previous is not None:
        metadata, content = previous[url]
        if digest(content.encode())!=metadata['sha256']:
            raise ValueError('Frozen query hash mismatch: ' + url)
        fn.write_bytes(content.encode())
        fn.with_suffix('.utc').write_text(metadata['retrieved_utc'])
    for attempt in range(3):
        try:
            if not fn.exists():
                with urlopen(url,timeout=80) as res:
                    raw=res.read()
                fn.write_bytes(raw)
                fn.with_suffix('.utc').write_text(datetime.now(timezone.utc).isoformat())
            raw=fn.read_bytes()
            catch=json.loads(raw)
            obs=catch['data']
            if set(catch['query']['sources']) != set(SITES) or len(obs)!=catch['count']:
                raise ValueError('Incomplete CATCH response')
            break
        except (OSError, ValueError, KeyError) as exc:
            if attempt==2:
                raise RuntimeError(f'{tic} {start}/{stop}: {exc}') from exc
            time.sleep(attempt+1)
    rec=dict(tic=tic,start_utc=str(start),stop_utc=str(stop), url=url,
             retrieved_utc=fn.with_suffix('.utc').read_text(),sha256=digest(raw),bytes=len(raw),
             counts_by_site=dict(Counter(r['source'] for r in obs)),
             ids=[i for i,_,_ in boxes],footprints=len(obs))
    return rec, obs, raw.decode('utf-8')


def time_bounds(obs, star):
    by_site=defaultdict(list)
    for r in obs:
        by_site[r['source']].append(r)
    shutter={}
    for site, records in by_site.items():
        times=np.asarray([[r['mjd_start'],r['mjd_stop']] for r in records])
        duration = (times[:,1]-times[:,0])*86400.
        if any(abs(seconds-r['exposure'])>.1 or seconds<=0 for seconds,r in zip(duration,records)):
            raise ValueError('Inconsistent ATLAS shutter times')
        t=Time(times,format='mjd',scale='utc',location=EarthLocation.from_geodetic(*SITES[site]))
        ras,decs=propagated(star,t.jyear)
        pos=SkyCoord(ras,decs,unit='deg')
        bjd=(t.tdb+t.light_travel_time(pos)).jd - 2457000
        for row,(start,stop) in zip(records,bjd):
            shutter[row['product_id']]=(float(start),float(stop))
    return shutter


def main(scan,cache,out,workers,frozen=None):
    windows=read_windows(scan)
    roster=list(csv.DictReader((ROOT/'extension/extension_parent.csv').open()))
    stars={r['tic_requested']:r for r in roster}
    cache.mkdir(exist_ok=True,parents=True)
    previous=None
    if frozen:
        meta=json.loads((frozen/'coverage.json').read_text())
        packed=(frozen/'responses.json.gz').read_bytes()
        if digest(packed)!=meta['responses_sha256']:
            raise ValueError('CATCH bundle hash mismatch')
        bodies=json.loads(gzip.decompress(packed))
        previous={r['url']:(r,body) for r,body in zip(meta['queries'],bodies,strict=True)}
    jobs=[(tic,boxes,stars,cache,previous) for tic,boxes in query_groups(windows)]
    with ThreadPoolExecutor(workers) as pool:
        responses=list(pool.map(query,jobs))
    queries=[]; photometry=[]
    for group,(rec,obs,_) in zip(jobs,responses):
        tic,boxes,_,_,_=group
        queries.append(rec)
        shutter=time_bounds(obs,stars[tic])
        for idx,left,right in boxes:
            w=windows[idx]
            common=[]
            nearest=None
            for r in obs:
                start,stop=shutter[r['product_id']]
                delta=max(0.,start-right,left-stop)
                distance=delta*86400.
                if nearest is None or distance<nearest:
                    nearest=distance
                overlap=min(right,stop)-max(left,start)
                if overlap>0:
                    common.append(dict(product_id=r['product_id'],site=r['source'],band=r['filter'],
                                       start_btjd=start,stop_btjd=stop,overlap_s=overlap*86400.))
            w.update(query_index=len(queries)-1,matches=common,nearest_shutter_gap_s=nearest)
        matches={x['product_id'] for idx,_,_ in boxes for x in windows[idx]['matches']}
        if matches:
            def night(r):
                # Local mean solar noon-to-noon observing night.
                return (r['source'],r['filter'],int(np.floor(.5*(r['mjd_start']+r['mjd_stop'])+SITES[r['source']][0]/360.-.5)))
            touched={night(r) for r in obs if r['product_id'] in matches}
            local=[r for r in obs if night(r) in touched]
            sub=cache/('pixels-'+tic);sub.mkdir(exist_ok=True)
            with ThreadPoolExecutor(workers) as pool:
                pixels=list(pool.map(lambda r:measure(r,stars[tic],sub),local))
            for r,p in zip(local,pixels):
                p.update(tic=tic,site=r['source'],local_night=night(r)[-1],in_any_box=r['product_id'] in matches)
                photometry.append(p)
    for pixel in photometry:
        if not pixel['in_any_box'] or pixel['status']!='measured':
            continue
        keys = ('tic', 'site', 'band', 'local_night')
        baseline = [p['f25'] for p in photometry if p['status']=='measured' and not p['in_any_box']
                    and all(p[k]==pixel[k] for k in keys)]
        pixel['n_same_night_comparison'] = len(baseline)
        if baseline and np.median(baseline)>0:
            median = float(np.median(baseline))
            pixel['same_night_median_f25'] = median
            pixel['depth_against_same_night'] = 1. - pixel['f25']/median
            pixel['same_night_relative_mad'] = float(1.4826*np.median(abs(np.array(baseline)/median - 1.)))
    compressed=gzip.compress(json.dumps(windows,separators=(',',':'),allow_nan=False).encode(),mtime=0)
    out.mkdir(exist_ok=True,parents=True)
    (out/'windows.json.gz').write_bytes(compressed)
    packed=gzip.compress(json.dumps([body for _,_,body in responses],separators=(',',':')).encode(),mtime=0)
    (out/'responses.json.gz').write_bytes(packed)
    ids={r['product_id'] for _,observations,_ in responses for r in observations}
    report=dict(responses_sha256=digest(packed),distinct_atlas_exposures=len(ids),window_counts=dict(Counter(w['sample'] for w in windows)),targets=len(set(w['tic'] for w in windows)),
                groups=len(queries),matching_boxes=sum(bool(w['matches']) for w in windows),
                matching_images=len({m['product_id'] for w in windows for m in w['matches']}),
                footprints=sum(r['footprints'] for r in queries), windows_sha256=digest(compressed),
                sites_lon_lat_height=SITES,queries=queries,matched_night_photometry=photometry)
    (out/'coverage.json').write_text(json.dumps(report,indent=1,allow_nan=False)+'\n')
    print({k:v for k,v in report.items() if k not in ('queries','matched_night_photometry','sites_lon_lat_height')})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('scan',type=Path)
    p.add_argument('cache',type=Path)
    p.add_argument('out',type=Path)
    p.add_argument('--workers',type=int,default=5)
    p.add_argument('--frozen',type=Path,help='recompute from the deposited CATCH response bundle')
    a=p.parse_args()
    main(a.scan,a.cache,a.out,a.workers,a.frozen)
