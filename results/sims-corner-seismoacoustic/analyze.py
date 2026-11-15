"""Reproduce the two bolide wavefront measurements from EarthScope miniSEED."""
import argparse
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from obspy import UTCDateTime, read
from scipy.optimize import least_squares
from scipy.signal import butter, sosfiltfilt

ROOT = Path(__file__).resolve().parent
CENTERS = np.arange(-1500, 2700, dtype=float) + 0.5
SHIFTS = [s for s in range(-1200, 1201, 20) if abs(s) >= 200]
BANDS = [(0.6, 6), (1, 4), (2, 8)]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def trace_path(cache, row):
    return cache / (row['event'] + '.' + row['nslc'] + '.ms')


def retrieve(cache, row):
    if not row['sha256']:
        return
    path = trace_path(cache, row)
    if not path.exists():
        with urlopen(row['url'], timeout=90) as response:
            payload = response.read()
        if hashlib.sha256(payload).hexdigest() != row['sha256']:
            raise ValueError('Returned miniSEED differs from the frozen bytes: ' + str(path))
        path.write_bytes(payload)
    if digest(path) != row['sha256']:
        raise ValueError('Checksum mismatch: ' + str(path))


def distance(lat, lon, lats, lons):
    a, b, c, d = np.deg2rad(lat), np.deg2rad(lon), np.deg2rad(lats), np.deg2rad(lons)
    return 12742 * np.arcsin(np.sqrt(np.sin((c-a)/2)**2 + np.cos(a)*np.cos(c)*np.sin((d-b)/2)**2))


def curves(path, epoch, band):
    """One-second RMS and centered four-second RMS, split at gaps."""
    sums = np.zeros(len(CENTERS))
    counts = np.zeros(len(CENTERS))
    for tr in read(str(path)).merge(method=-1).split():
        data = np.asarray(tr.data, dtype=float)
        rate = float(tr.stats.sampling_rate)
        if len(data) <= 30 * rate:
            continue
        filt = butter(3, band, btype='bandpass', fs=rate, output='sos')
        y = sosfiltfilt(filt, data)
        seconds = np.arange(len(y)) / rate + float(tr.stats.starttime - epoch)
        idx = np.floor(seconds - CENTERS[0] + .5).astype(int)
        valid = (idx >= 0) & (idx < len(CENTERS))
        # Avoid filter/record boundaries and require full-second sample coverage.
        edge = int(10 * rate)
        valid[:edge] = False
        valid[-edge:] = False
        bins = np.bincount(idx[valid], minlength=len(CENTERS))
        full = bins >= .98 * rate
        sums[full] += np.bincount(idx[valid], weights=y[valid]**2, minlength=len(CENTERS))[full]
        counts[full] += bins[full]
    one = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    four = np.convolve(one, np.ones(4)/4, mode='same')
    # For [i-2,i-1,i,i+1], centroid time is CENTERS[i] - 0.5.
    return np.sqrt(one), np.sqrt(four)


def maximum(y, times, lo, hi):
    indices = np.flatnonzero((times >= lo) & (times <= hi) & np.isfinite(y))
    if not len(indices):
        return np.nan, np.nan
    index = indices[np.argmax(y[indices])]
    return float(times[index]), float(y[index])


def gate(row, alt):
    slant = np.hypot(row['range_km'], alt)
    return slant / .36 - 10, slant / .27 + 10


def model(params, rows):
    lat, lon, alt, speed = params
    ranges = distance(lat, lon, [r['lat'] for r in rows], [r['lon'] for r in rows])
    return np.hypot(ranges, alt - np.array([r['elev_km'] for r in rows])) / speed


def summarize(row, alt, one, four):
    lo, hi = gate(row, alt)
    t = CENTERS - .5
    peak, val = maximum(four, t, lo, hi)
    bg = (t < lo - 90) | (t > hi + 90)
    med = float(np.nanmedian(four[bg]))
    outside = float(np.nanmax(four[bg]))
    i, amp = maximum(one, CENTERS, lo, hi)
    return {k:row[k] for k in ('event','nslc','range_km','lat','lon','elev_km')} | {
        'peak_s':peak, 'peak_to_median':val/med, 'peak_to_offgate_max':val/outside,
        'background_median':med, 'offgate_max':outside, 'one_sec_peak_s':i,
        'one_sec_peak_to_median':amp/float(np.nanmedian(one[bg])),
        'coverage_fraction':float(np.mean(np.isfinite(four))),
    }


def write_csv(path, rows):
    with path.open('w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def wa_holdout(rows, energy):
    # Near-station measurements select calibration. Far stations are all candidates.
    near = [r for r in rows if r['event']=='wa' and r['range_km']<=70 and
            r['peak_to_median'] >= 6 and r['peak_to_offgate_max'] > 1]
    far = [r for r in rows if r['event']=='wa' and r['range_km']>70]
    start = [47.7, -119.4, 30, .32]
    bounds = ([47.3, -119.8, 10, .26], [48.1, -119., 70, .40])
    predicted = np.array([r['peak_s'] for r in near])
    fit = least_squares(lambda p:model(p,near)-predicted, start, bounds=bounds)
    times = model(fit.x, far)
    # Test the near-source fit with a fixed +/-25 s tolerance at all far stations.
    tests = []
    ctrl = []
    for shift in [0] + SHIFTS:
        hits = 0
        for r, pred in zip(far,times):
            _, four = energy['wa',r['nslc']]
            t, v = maximum(four, CENTERS-.5, pred+shift-25, pred+shift+25)
            ratio = v/r['background_median']
            hits += int(ratio>=6)
            if shift == 0:
                tests.append({'nslc':r['nslc'], 'range_km':r['range_km'],
                              'predicted_s':float(pred),'window_peak_s':t,'peak_to_median':ratio,
                              'broad_gate_peak_minus_prediction_s':r['peak_s']-pred})
        ctrl.append({'shift_s':shift,'count':hits})
    # Compare fixed-axis reference positions under the same free celerity range.
    comparator = []
    for name, lat, lon, alt in [('CNEOS peak',47.7,-119.4,30.0),
                               ('Skyfalls optical endpoint',47.85,-119.408,45.5)]:
        dist = model([lat,lon,alt,1.0],near)
        speed = float(np.clip(np.dot(dist,dist)/np.dot(dist,predicted),.26,.40))
        comparator.append({'source':name,'best_speed_kmps':speed,
                           'rms_s':float(np.sqrt(np.mean((dist/speed-predicted)**2)))})
    return {'calibrators':[r['nslc'] for r in near],'fixed_position_models':comparator,
            'lat_lon_heightkm_speedkmps':fit.x.tolist(),
            'fit_rms_s':float(np.sqrt(np.mean(fit.fun**2))), 'test_window_half_width_s':25,
            'holdout':tests,'shift_counts':ctrl}


def alaska(rows, energy):
    # All sites with pressure and vertical ground-motion data inside 200 km.
    ak = {r['nslc']:r for r in rows if r['event']=='ak'}
    pairs = []
    for key,p in ak.items():
        if not key.endswith('BDF'):
            continue
        z = ak.get(key.split('.EP.')[0]+'..BHZ')
        if z:
            pairs.append((p,z))
    results, controls = [], []
    for shift in [0] + SHIFTS:
        hits = 0
        for p,z in pairs:
            lo,hi = gate(p,35.5)
            tp,vp = maximum(energy['ak',p['nslc']][1], CENTERS-.5, lo+shift, hi+shift)
            tz,vz = maximum(energy['ak',z['nslc']][1], CENTERS-.5, lo+shift, hi+shift)
            rp,rz = vp/p['background_median'],vz/z['background_median']
            yes = abs(tp-tz)<=3 and min(rp,rz)>=4
            hits += int(yes)
            if shift == 0:
                local = abs((CENTERS-.5) - tp) <= 15
                ep = energy['ak',p['nslc']][0][local]
                ez = energy['ak',z['nslc']][0][local]
                usable = np.isfinite(ep) & np.isfinite(ez)
                corr = float(np.corrcoef(ep[usable],ez[usable])[0,1])
                results.append({'station':p['nslc'].split('.')[1], 'range_km':p['range_km'],
                                'envelope_correlation_30s':corr,
                                'pressure_peak_s':tp, 'vertical_peak_s':tz,
                                'pressure_rms_ratio':rp,'vertical_rms_ratio':rz,'passes':bool(yes)})
        controls.append({'shift_s':shift,'count':hits})
    return {'colocated':results,'shift_counts':controls}


def calendar_controls(cache, manifest):
    records = json.loads((ROOT/'calendar-controls.json').read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda r:retrieve(cache,r), records))
    result = []
    for event in ('wa','ak'):
        source = manifest['sources'][event]['cneos']
        # On-event measurements use the same short window as the placebo dates.
        template = [r for r in records if r['real_event']==event and r['days']==-1]
        selected = {r['nslc'] for r in template}
        live = [dict(r,days=0,epoch=source['date']) for r in manifest['waveforms']
                if r['event']==event and r['nslc'] in selected]
        for day in (0,-7,-1,1,7):
            cohort = live if day==0 else [r for r in records if r['real_event']==event and r['days']==day]
            peaks = {}
            for r in cohort:
                if not r['sha256']:
                    continue
                one,four = curves(trace_path(cache,r),UTCDateTime(r['epoch']),BANDS[0])
                # Enforce the same 10-s filter-edge exclusion at all epochs.
                invalid = (CENTERS < -590) | (CENTERS > 770)
                one[invalid] = np.nan
                four[invalid] = np.nan
                peaks[r['nslc']]=summarize(r,float(source['alt']),one,four)
            if event=='wa':
                passed=sum(r['peak_to_median']>=6 and r['peak_to_offgate_max']>1 for r in peaks.values())
                tested=len(peaks)
            else:
                passed=tested=0
                for nslc, pressure in peaks.items():
                    ground=peaks.get(nslc.split('.EP.')[0]+'..BHZ')
                    if nslc.endswith('BDF') and ground:
                        tested+=1
                        passed+=int(abs(pressure['peak_s']-ground['peak_s'])<=3 and
                                    min(pressure['peak_to_median'],ground['peak_to_median'])>=4)
            result.append({'event':event,'days':day,'passing_sites':int(passed),'covered_sites':tested})
    return result


def sounding():
    series = {}
    for line in (ROOT/'spokane-igra.txt').read_text().splitlines():
        if line.startswith('#'):
            k = line[:26].strip()
            series[k] = []
            continue
        height, temp, wind = (int(line[a:b]) for a,b in [(16,21),(22,27),(46,51)])
        if min(height,wind)>=0 and temp > -1500:
            c = np.sqrt(1.4*287.05*(273.15+temp/10))
            series[k].append([height, c, wind/10, c+wind/10])
    return {k:{'levels':len(v),'max_height_m':max(r[0] for r in v),
               'max_c_plus_wind_mps':max(r[3] for r in v)} for k,v in series.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True,parents=True)
    args.output.mkdir(exist_ok=True,parents=True)
    manifest = json.loads((ROOT/'manifest.json').read_text())
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda r:retrieve(args.cache,r), manifest['waveforms']))
    output = {}
    for band in BANDS:
        energy,rows = {},[]
        for row in manifest['waveforms']:
            if not row['sha256']:
                continue
            source = manifest['sources'][row['event']]['cneos']
            key = row['event'],row['nslc']
            energy[key] = curves(trace_path(args.cache,row),UTCDateTime(source['date']),band)
            rows.append(summarize(row,float(source['alt']),*energy[key]))
        k = '-'.join(map(str,band))
        write_csv(args.output/('stations-'+k+'.csv'),rows)
        output[k]={'wa':wa_holdout(rows,energy), 'ak':alaska(rows,energy)}
    output['calendar_controls']=calendar_controls(args.cache,manifest)
    output['spokane_soundings']=sounding()
    (args.output/'analysis.json').write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
