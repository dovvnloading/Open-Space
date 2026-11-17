"""Measure every flagged deep-dip window in a separate, Gaia-centred SPOC pixel aperture."""
import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from ffi_acquire import GAIA_FIELDS, download, propagated, tsv_rows, vizier
from ffi_images import event_masks
from ffi_photometry import deposit, gaia_t
from pixel_evidence import profile_center
from replication import direct
from short_pixels import acquire


def field(star, cache, snapshots):
    name = star['GaiaEDR3'] + '-gaia.tsv'
    path = cache / name
    if name in snapshots:
        path.write_text(snapshots[name])
        source = {'archive': 'extension/extension_catalogue_responses.json.gz', 'member': name}
    else:
        url = vizier('I/355/gaiadr3', GAIA_FIELDS,
                     **{'-c': f"{star['RA_ICRS']} {star['DE_ICRS']}", '-c.rs': '180'})
        download(path, url)
        source = json.loads(path.with_suffix('.tsv.json').read_text())
    source['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    stars = tsv_rows(path)
    if sum(r['Source'] == star['GaiaEDR3'] for r in stars) != 1:
        raise ValueError('Gaia cone does not contain one target row')
    return stars, source


def dilution(star, neighbours, wcs, year, x, y, sx, sy):
    coords = [wcs.world_to_pixel_values(*propagated(r, year)) for r in neighbours]
    fractions, contributors = [], []
    for sigma in (.65, .85, 1.05):
        throughput = deposit(x, y, sx, sy, sigma)
        for adverse in (False, True):
            ratios = []
            for row, (nx, ny) in zip(neighbours, coords):
                if row['Source'] == star['GaiaEDR3'] or not row['Gmag']:
                    continue
                ratio = 10**(-.4*(gaia_t(row, adverse)-float(star['tmag'])-(.3 if adverse else 0)))
                ratio *= deposit(x, y, nx, ny, sigma)/throughput
                ratios.append(ratio)
                if sigma == .85 and not adverse:
                    contributors.append(dict(gaia=row['Source'], x=float(nx), y=float(ny), ratio=ratio))
            fractions.append(1/(1+sum(ratios)))
    return dict(nominal_fraction=fractions[2], minimum_fraction=min(fractions),
                neighbours_without_G=sum(not n['Gmag'] for n in neighbours),
                neighbours=sorted(contributors, key=lambda r: -r['ratio'])[:5])


def measure(rec, t, cube, aperture, x, y, sx, sy):
    period = rec.get('period_days')
    duration = rec['duration_min']/1440.
    inside, left, right = event_masks(t, rec['epoch_btjd'], duration, period)
    result = dict(in_points=int(inside.sum()), left_points=int(left.sum()), right_points=int(right.sum()))
    if not (inside.any() and left.any() and right.any()):
        return dict(result, status='incomplete exposure coverage')
    if period:
        # Equal weight per supported orbital event, and equal weight per shoulder.
        groups = np.rint((t-rec['epoch_btjd'])/period).astype(int)
        cycles = np.intersect1d(np.intersect1d(groups[inside], groups[left]), groups[right])
    else:
        groups = np.zeros(len(t), dtype=int)
        cycles = [0]
    if len(cycles) == 0:
        return dict(result, status='no cycle with both shoulders')
    stack = []
    for k in cycles:
        sel = groups == k
        means = [np.nanmean(cube[sel & place], axis=0) for place in (inside, left, right)]
        stack.append(means)
    center, before, after = np.nanmean(stack, axis=0)
    reference = .5*(before+after)
    delta = reference-center
    baseline = reference[aperture].sum()
    if not np.isfinite(baseline) or baseline <= 0 or not np.isfinite(center[aperture]).all():
        return dict(result, status='invalid aperture')
    fit = profile_center(delta)
    def photocenter(image):
        weight = image[aperture]
        if weight.sum() <= 0:
            return [None, None]
        return [float(np.sum(axis[aperture]*weight)/weight.sum()) for axis in (x,y)]
    context = (abs(t-rec['epoch_btjd'])<max(.12, 5*duration)) & (abs(t-rec['epoch_btjd'])>2*duration)
    if period is None and context.sum()>=5:
        local = float(np.nanmedian(cube[context][:,aperture].sum(axis=1)))
        if local > 0:
            result.update(center_vs_local=float(center[aperture].sum()/local-1),
                          shoulders_vs_local=float(baseline/local-1))
    return dict(result, status='measured', cycles=len(cycles),
                depth=float(delta[aperture].sum()/baseline),
                left_minus_right=float((before-after)[aperture].sum()/baseline),
                in_centroid=photocenter(center), out_centroid=photocenter(reference),
                difference_image_e_s=[[float(v) if np.isfinite(v) else None for v in row] for row in delta],
                difference_xy=fit.tolist() if fit is not None else None,
                offset_pixel=float(np.hypot(*(fit-[sx,sy]))) if fit is not None else None)


def visit(job):
    row, triggers, star, neighbours, meta, out, cache = job
    dest = out / (row['filename']+'.json')
    signature = hashlib.sha256(json.dumps(triggers,sort_keys=True).encode()).hexdigest()
    if dest.exists():
        previous=json.loads(dest.read_text())
        if previous.get('trigger_sha256')==signature and 'error' not in previous:
            return
    url = direct(row['filename'])
    rec = dict(product=row['filename'], tic=row['tic'], sector=int(row['sector']),
               cadence_s=int(row['cadence_s']), source_id=star['GaiaEDR3'], trigger_sha256=signature,
               gaia_query=meta, url=url, retrieved_utc=datetime.now(timezone.utc).isoformat())
    try:
        with TemporaryDirectory(dir=cache) as tmp:
            path=Path(tmp)/row['filename']
            rec.update(sha256=acquire(path,url), bytes=path.stat().st_size)
            with fits.open(path,memmap=True) as h:
                if int(h[0].header['TICID'])!=int(row['tic']) or h[0].header['SECTOR']!=rec['sector']:
                    raise ValueError('pixel identity mismatch')
                data = h[1].data
                if abs(h[1].header['TIMEDEL']*86400-rec['cadence_s'])>.1:
                    raise ValueError('pixel cadence mismatch')
                ok=(data['QUALITY']==0) & np.isfinite(data['TIME'])
                t=data['TIME'][ok].astype(float)
                year=2000+(np.median(t)+2457000.-2451545.)/365.25
                wcs=WCS(h[2].header,naxis=2)
                sx,sy=map(float,wcs.world_to_pixel_values(*propagated(star,year)))
                y,x=np.indices(h[2].data.shape)
                ap=np.hypot(x-sx,y-sy)<=1.5
                if not ap.any() or min(sx,sy,x.max()-sx,y.max()-sy)<1.5:
                    raise ValueError('truncated aperture')
                cube=data['FLUX'][ok]
                ok=np.isfinite(cube[:,ap]).all(axis=1)
                t,cube=t[ok],cube[ok]
                rec.update(target_xy=[sx,sy], aperture_pixels=int(ap.sum()), quality0_finite=len(t),
                           **dilution(star,neighbours,wcs,year,x[ap],y[ap],sx,sy))
                rec['events']=[]
                for ev in triggers:
                    metric=measure(ev,t,cube,ap,x,y,sx,sy)
                    if 'depth' in metric:
                        metric['depth_if_on_wd']=metric['depth']/rec['nominal_fraction']
                    rec['events'].append(dict(trigger=ev,pixels=metric))
    except Exception as exc:
        rec['error']=repr(exc)
    dest.write_text(json.dumps(rec,allow_nan=False)+'\n')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('roster','scan','cache','out'):
        p.add_argument(key,type=Path)
    p.add_argument('--workers',type=int,default=3)
    a=p.parse_args()
    for path in (a.cache,a.out):
        path.mkdir(exist_ok=True,parents=True)
    rows=list(csv.DictReader(a.roster.open()))
    available=defaultdict(list)
    for r in rows:
        available[r['tic'],int(r['sector'])].append(r)
    events=defaultdict(list)
    if a.scan.is_dir():
        scans=[json.loads(path.read_text()) for path in sorted(a.scan.glob('*.json'))]
    else:
        with gzip.open(a.scan,'rt') as stream:
            scans=json.load(stream)
    for r in scans:
        for kind in ('periodic','isolated'):
            events[r['tic'],r['sector']].extend(dict(kind=kind,lc=r['filename'],**ev) for ev in r.get(kind,[]))
    stars={r['GaiaEDR3']:r for r in tsv_rows(Path('parent.tsv'))}
    matches={r['tic']:dict(stars[r['gaia']],tmag=r['tmag']) for r in csv.DictReader(open('tic_match.csv')) if r['status']=='source_id_and_position'}
    with gzip.open('extension/extension_catalogue_responses.json.gz','rt') as stream:
        snapshots=json.load(stream)
    hosts=sorted({tic for (tic,sector),items in events.items() if items})
    fields={tic:field(matches[tic],a.cache,snapshots) for tic in hosts}
    jobs=[]
    for (tic,s),ev in events.items():
        if ev:
            row=max(available[tic,s],key=lambda r:int(r['cadence_s']))
            jobs.append((row,ev,matches[tic],*fields[tic],a.out,a.cache))
    with ThreadPoolExecutor(a.workers) as pool:
        list(pool.map(visit,jobs))

if __name__=='__main__':
    main()
