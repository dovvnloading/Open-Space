"""Fit known transits to independent JWST white-light products. Cache is kept outside the results."""
import argparse
import csv
import hashlib
import io
import sys
from pathlib import Path
from urllib.request import urlopen
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
DURATIONS = [36.06, 42.03, 48.87, 55.76, 62.85, 68.24, 76.7]
SCALE = 0.0015


def load_product(row, cache):
    path = cache / row['filename']
    if not path.exists():
        url = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:JWST/product/' + row['filename']
        path.write_bytes(urlopen(url,timeout=90).read())
    data = path.read_bytes()
    if hashlib.sha256(data).hexdigest() != row['sha256']:
        raise ValueError(f'Product revision does not match manifest: {path}')
    return np.loadtxt(line for line in io.BytesIO(data) if line[:1].isdigit())


def fit(x, y, duration, degree, bootstrap):
    def model(p):
        baseline = np.polynomial.polynomial.polyval(x / 100, p[:degree+1])
        depth, delta, length, ingress = p[degree+1:]
        h = np.clip((length/2 - abs(x-delta)) / ingress, 0, 1)
        return baseline - depth * h*h*(3-2*h)
    lo = [.8] + [-.6]*degree + [0, -12, duration*.8, duration*.04]
    hi = [1.2] + [.6]*degree + [.03, 12, duration*1.2, duration*.18]
    start = [1.] + [0.]*degree + [.006, 0., duration, duration*.09]
    def solve(z, initial):
        return least_squares(lambda p: model(p)-z, initial, bounds=(lo,hi),loss='soft_l1',f_scale=SCALE)
    sol = solve(y, start)
    residual = y-model(sol.x)
    # Wild block bootstrap. One sign per elapsed-time block preserves within-block covariance.
    # Gaps or masked neighboring planets retain their original elapsed-time labels.
    rng = np.random.default_rng(6422)
    errors = []
    if bootstrap:
        for minutes in [3,9]:
            blocks = np.floor((x-x.min())/minutes).astype(int)
            offsets=[]
            for _ in range(60):
                signs = rng.choice([-1,1],size=blocks.max()+1)
                p = solve(model(sol.x)+residual*signs[blocks],sol.x).x
                offsets.append(60*p[degree+2])
            errors.append(float((np.percentile(offsets,84)-np.percentile(offsets,16))/2))
    return 60*sol.x[degree+2], 1e6*sol.x[degree+1], sol.x[degree+3], errors


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('cache',type=Path)
    args=ap.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    with (HERE/'products.csv').open() as f:
        products = list(csv.DictReader(f))
    with (HERE/'forecast_subset.csv').open() as f:
        forecast = list(csv.DictReader(f))
    output=[]
    for row in products:
        raw = load_product(row,args.cache)
        # Both archive time columns are reduced Julian dates despite the BJD label.
        # First average in fixed 20-second bins to reduce correlated oversampling and bound fit cost.
        raw = raw[np.isfinite(raw).all(axis=1)]
        ids = np.floor((raw[:,1]-raw[0,1])*4320).astype(int)
        starts = np.r_[0,np.flatnonzero(np.diff(ids))+1]
        t = np.add.reduceat(raw[:,1],starts)/np.diff(np.r_[starts,len(ids)]) - 49999.5
        y = np.array([np.median(block) for block in np.split(raw[:,2],starts[1:])]); y /= np.median(y)
        selected = [r for r in forecast if min(t)-.04<float(r['bjd_minus_2450000'])<max(t)+.04]
        for current in selected:
            planet = current['planet']
            dur = DURATIONS[ord(planet)-ord('b')]
            mid=float(current['bjd_minus_2450000']); x=(t-mid)*1440
            good=(abs(x)<dur*1.8) & (abs(y-1)<.1)
            for other in selected:
                if other['planet']!=planet:
                    neighbor_duration=DURATIONS[ord(other['planet'])-ord('b')]
                    good &= abs(t-float(other['bjd_minus_2450000']))*1440>neighbor_duration/2+3
            if np.sum(good)<100:
                continue
            local, depth, length, errors = fit(x[good],y[good],dur,2,True)
            alternative,_,_,_ = fit(x[good],y[good],dur,1,False)
            output.append([row['filename'][:12],planet,current['epoch'],f'{mid:.9f}',
                           f'{2450000+mid+local/86400:.9f}',round(local,2),round(alternative,2),
                           round(errors[0],2),round(errors[1],2),round(float(current['sigma_days'])*86400,2),
                           round(depth),round(length,2)])
    with (HERE/'times.csv').open('w') as f:
        writer=csv.writer(f)
        writer.writerow(['visit','planet','epoch','forecast_bjd_minus_2450000','measured_bjd_tdb',
                         'oc_seconds','linear_baseline_oc_seconds','block_3min_error_seconds','block_9min_error_seconds',
                         'forecast_error_seconds','depth_ppm','duration_minutes'])
        writer.writerows(output)

if __name__=='__main__':
    main()
