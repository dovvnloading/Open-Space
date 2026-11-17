"""Remeasure each frozen QLP trigger on a SPOC full-frame 9x9 TESScut stamp."""
import argparse
import csv
import hashlib
import io
import json
import re
import time
import zipfile
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np


def header(payload, offset):
    values = {}
    while offset + 80 <= len(payload):
        line = payload[offset:offset + 80].decode('ascii')
        offset += 80
        if line[:8].strip() == 'END':
            return values, ((offset + 2879) // 2880) * 2880
        if line[8:10] == '= ':
            values[line[:8].strip()] = line[10:].split('/')[0].strip().strip("'").strip()
    raise ValueError('truncated FITS header')


def read_stamp(payload):
    h, ext = header(payload, 0)
    if h['NAXIS'] != '0':
        raise ValueError('expected empty primary HDU')
    e, offset = header(payload, ext)
    formats = {'J': '>i4', 'I': '>i2', 'E': '>f4', 'D': '>f8', 'K': '>i8', 'B': 'u1'}
    cols = []
    for i in range(1, int(e['TFIELDS']) + 1):
        name = e[f'TTYPE{i}']
        repeat, code = re.fullmatch(r'(\d*)([JIEDKBA])', e[f'TFORM{i}']).groups()
        size = int(repeat or '1')
        if code == 'A':
            col = (name, f'S{size}')
        elif size == 1:
            col = (name, formats[code])
        else:
            col = (name, formats[code], (size,))
        cols.append(col)
    dtype = np.dtype(cols)
    if dtype.itemsize != int(e['NAXIS1']):
        raise ValueError('unexpected FITS table width')
    return h, e, np.frombuffer(payload, dtype, int(e['NAXIS2']), offset)


def target_pixel(e, ra, dec):
    # FITS celestial binary-table WCS for column 5, 0-based result (x,y).
    vals = [float(e[f'{k}CRVL5']) for k in (1, 2)]
    a, d, a0, d0 = np.radians([ra, dec, *vals])
    da = a - a0
    cosc = np.sin(d0)*np.sin(d) + np.cos(d0)*np.cos(d)*np.cos(da)
    uv = np.degrees([np.cos(d)*np.sin(da)/cosc,
                     (np.cos(d0)*np.sin(d)-np.sin(d0)*np.cos(d)*np.cos(da))/cosc])
    pc = np.array([[float(e[f'{i}{j}PC5']) for j in (1,2)] for i in (1,2)])
    scales = np.array([float(e[f'{i}CDLT5']) for i in (1,2)])
    ref = np.array([float(e[f'{i}CRPX5'])-1 for i in (1,2)])
    return ref + np.linalg.solve(scales[:,None]*pc, uv)


def flux_fit(t, y, peak, duration, ingress):
    # Same triangle as the selection. Fit a local slope and intercept per series.
    u = (t - peak) / duration
    use = (u >= -ingress-.5) & (u <= 1.5-ingress) & np.isfinite(y).all(axis=1)
    u = u[use]
    if min(sum(u < -ingress), sum(u > 1-ingress), sum((u >= -ingress) & (u <= 1-ingress))) < 5:
        return None
    profile = np.maximum(0, np.minimum(1+u/ingress, 1-u/(1-ingress)))
    design = np.column_stack((np.ones(len(u)), u, -profile))
    coefficients = np.linalg.lstsq(design, y[use], rcond=None)[0]
    return coefficients[0], coefficients[2], len(u)


def measure(stamp, event, xy):
    t = stamp['TIME']
    keep = (stamp['QUALITY']==0) & np.isfinite(t)
    t = t[keep]
    image = stamp['FLUX'][keep].astype(float)
    ny = nx = 9
    yy, xx = (v.ravel() for v in np.indices((ny,nx)))
    edge = (yy==0)|(yy==ny-1)|(xx==0)|(xx==nx-1)
    # Local per-cadence perimeter median removes an additive sky level.
    sky = np.nanmedian(image[:,edge],axis=1)
    sub = image-sky[:,None]
    aperture = (abs(xx-round(xy[0]))<=1) & (abs(yy-round(xy[1]))<=1)
    aperture_series = sub[:,aperture].sum(axis=1)
    ref = np.column_stack([aperture_series, sky])
    peak, dur, frac = [float(event[k]) for k in ('peak_btjd','duration_h','ingress_fraction')]
    dur /= 24
    fit = flux_fit(t, ref, peak, dur, frac)
    if fit is None:
        return dict(status='insufficient pixel coverage')
    baseline, depth, count = fit
    relative = depth[0]/baseline[0] if baseline[0]>0 else None
    weights = []
    for other in np.arange(t.min()+dur,t.max()-dur,.25):
        if abs(other-peak)<dur: continue
        control = flux_fit(t,ref,other,dur,frac)
        if control and control[0][0]>0:
            weights.append(control[1][0]/control[0][0])
    imgfit = flux_fit(t, sub, peak, dur, frac)
    metrics = dict(status='ok',points=count,aperture_fractional_depth=relative,
                   aperture_baseline_e_per_s=baseline[0],aperture_deficit_e_per_s=depth[0],
                   perimeter_deficit_e_per_s_per_pixel=depth[1],control_windows=len(weights),
                   controls_as_deep=sum(v>=relative for v in weights) if relative is not None else None,
                   control_depth_p95=float(np.percentile(weights,95)) if weights else None)
    if imgfit:
        diff = imgfit[1]
        # Heuristic location within radius 3 of target. Positive pixels only.
        core=(xx-xy[0])**2+(yy-xy[1])**2<=9
        weights=np.maximum(0, diff-np.median(diff[edge]))*core
        if weights.sum()>0:
            cx,cy=np.average(xx,weights=weights),np.average(yy,weights=weights)
            metrics.update(positive_deficit_centroid_x=cx,positive_deficit_centroid_y=cy,
                           centroid_offset_px=float(np.hypot(cx-xy[0],cy-xy[1])))
        metrics['deficit_image_e_per_s']=np.round(diff.reshape(ny,nx),3).tolist()
    return metrics


def get_star(item):
    tic, events = item
    r=events[0]
    params=dict(ra=r['ra_obj'],dec=r['dec_obj'],sector=r['sector'],x=9,y=9,units='px')
    url='https://mast.stsci.edu/tesscut/api/v0.1/astrocut?'+urlencode(params)
    receipt=dict(tic=tic,url=url,retrieved_utc=datetime.now(timezone.utc).isoformat())
    for attempt in range(3):
        try:
            with urlopen(url,timeout=150) as resp:
                raw=resp.read()
            with zipfile.ZipFile(io.BytesIO(raw)) as z:
                files=[p for p in z.namelist() if p.endswith('.fits')]
                if len(files)!=1: raise ValueError(f'unexpected cutout members {files}')
                payload=z.read(files[0])
            h,e,stamp=read_stamp(payload)
            if int(h['SECTOR'])!=int(r['sector']): raise ValueError('sector mismatch')
            xy=target_pixel(e,float(r['ra_obj']),float(r['dec_obj']))
            if not all((xy>=0)&(xy<9)): raise ValueError(f'target outside stamp {xy}')
            receipt.update(member=files[0],fits_bytes=len(payload),fits_sha256=hashlib.sha256(payload).hexdigest(),
                           flux_sha256=hashlib.sha256(stamp['FLUX'].tobytes()).hexdigest(),target_xy=xy.tolist(),
                           camera=h.get('CAMERA'),ccd=h.get('CCD'),status='ok',
                           events=[{'peak_btjd':float(v['peak_btjd']),**measure(stamp,v,xy)} for v in events])
            break
        except Exception as err:
            receipt.update(status='error',error=str(err))
            if isinstance(err,(KeyError,ValueError)):break
            time.sleep(2**attempt)
    print(tic,receipt['status'],receipt.get('error',''),flush=True)
    return receipt


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('triggers',type=Path)
    p.add_argument('out',type=Path)
    args=p.parse_args()
    rows={}
    for r in csv.DictReader(args.triggers.open()):
        rows.setdefault(r['tic'],[]).append(r)
    with ThreadPoolExecutor(3) as pool:
        result=list(pool.map(get_star,rows.items()))
    args.out.write_text(json.dumps(result,indent=1,allow_nan=False,default=lambda v:v.item())+'\n')


if __name__=='__main__':
    main()
