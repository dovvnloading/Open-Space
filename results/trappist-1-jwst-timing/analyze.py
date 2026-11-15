"""Audit JWST DDT 12492 for the proposed 81-minute transit using real flux injections."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
DURATION = 81.0
KNOWN = [(2461021.1376, 36.06), (2461021.2406, 55.76)]  # Agol forecast, b and e


def shape(time, center, duration=DURATION, ingress=8.1):
    return np.clip((duration / 2 - np.abs(time - center)) / ingress, 0, 1)


def filter_matrix(time, n, known, centers, degree, baseline, duration=DURATION):
    """Local least squares filters for depth, fitted simultaneously with a smooth baseline."""
    weights = np.zeros((len(centers), len(time)))
    eligible = np.zeros(len(centers), dtype=bool)
    for i, center in enumerate(centers):
        x = time - center
        half = duration / 2
        select = ~known & (np.abs(x) < half + baseline)
        tpl = shape(time[select], center, duration, duration / 10)
        # A full footprint, at least 60% of the in-transit template, and two baseline flanks.
        inside = (time.min() <= center - half and time.max() >= center + half)
        supported = (np.sum(tpl) >= .6 * duration * .9 and
                     np.sum(select & (x < -half)) >= 20 and
                     np.sum(select & (x > half)) >= 20)
        if not inside or not supported:
            continue
        design = np.column_stack([((x[select] / 100) ** deg) for deg in range(degree + 1)])
        w = np.minimum(n[select], 44) / 44
        residual = tpl - design @ np.linalg.lstsq(design * np.sqrt(w[:, None]), tpl * np.sqrt(w), rcond=None)[0]
        weights[i, select] = -residual * w / np.dot(tpl, residual * w)
        eligible[i] = True
    return weights, eligible


def main():
    a = np.genfromtxt(HERE/'one_minute.csv', delimiter=',', names=True)
    bjd = a['bjd_tdb_minus_2400000_5']
    t = (bjd - bjd[0]) * 1440
    y = a['relative_flux']
    n = a['count']
    # A large terminal flare truncates the clean record. The last eight minutes remain in the injection denominator.
    known = t > t.max() - 8.0
    for center, duration in KNOWN:
        known |= np.abs(bjd - (center - 2400000.5)) * 1440 <= duration/2 + 8
    centers = np.arange(0, np.floor(t[-1]) + 1, 2.)
    injected = np.arange(1, np.floor(t[-1]), 2.)  # odd-minute locations, never on search grid
    models = [('primary', 1, 80), ('quadratic', 2, 100)]
    stats = {}
    grid = {}
    for label, degree, baseline in models:
        W, ok = filter_matrix(t, n, known, centers, degree, baseline)
        W = W[ok]
        scanned = centers[ok]
        depth = W @ (y - 1)
        # Inverted-flux scan has the same baseline, masked points and template as the dip scan.
        # A fixed 0.5 ppt floor prevents a single unusually quiet visit from lowering the threshold.
        threshold = max(.0005, -depth.min())
        # Inject at the 1-minute photometry stage before any local baseline is fitted.
        # The linear filtering matrix is algebraically identical to rerunning every fit.
        shapes = shape(t[:,None], injected[None,:], DURATION, DURATION*.14)
        transfer = -W @ (y[:,None]*shapes)
        rank = depth.argmax()
        stats[label] = {'centers_searched':len(scanned), 'max_dip_ppm':float(depth[rank]*1e6),
                        'max_dip_bjd_minus_2400000_5':float(bjd[0]+scanned[rank]/1440),
                        'max_brightening_ppm':float(-depth.min()*1e6), 'threshold_ppm':float(threshold*1e6)}
        loc = np.abs(scanned[:,None]-injected[None,:])
        # Denominator contains equally spaced centers over the actual integration span.
        results = []
        for ppm in [200,350,500,750,1000,1500,2000,5000,10000,20000]:
            amps = depth[:,None] + ppm*1e-6*transfer
            peak = np.argmax(amps,axis=0)
            recovered = ((amps[peak,np.arange(len(peak))] > threshold) &
                         (loc[peak,np.arange(len(peak))] <= 15))
            results.append([ppm, int(recovered.sum()),len(recovered)])
            grid[(label,ppm)] = recovered
        stats[label]['injections'] = results
        eligible = np.min(loc, axis=0) < 1.1
        stats[label]['eligible_injections'] = int(eligible.sum())
        stats[label]['eligible_recoveries'] = [
            [ppm, int(np.sum(grid[(label, ppm)] & eligible))] for ppm, _, _ in results
        ]
    with (HERE/'injections.csv').open('w') as f:
        w = csv.writer(f); w.writerow(['bjd_tdb_center','minute_from_first_bin']+
                [f'{label}_{ppm}ppm_recovered' for label, ppm in grid])
        for i,v in enumerate(injected):
            w.writerow([f'{bjd[0]+2400000.5+v/1440:.8f}',int(v)]+[int(z[i]) for z in grid.values()])
    # Sensitivity to the terminal-cut choice and to the proposed duration.
    # Each refit uses the exact same archived photometry.
    control_rows = []
    planet_mask = np.zeros(len(t), dtype=bool)
    for center, duration in KNOWN:
        planet_mask |= np.abs(bjd - (center - 2400000.5)) * 1440 <= duration/2 + 8
    for terminal_minutes in [0, 6, 8, 12, 20]:
        for duration in [60., 81., 105.]:
            mask = planet_mask | (t > t.max() - terminal_minutes)
            W, ok = filter_matrix(t, n, mask, centers, 1, 80, duration)
            depths = W[ok] @ (y - 1)
            control_rows.append([terminal_minutes, duration, int(ok.sum()),
                                 round(float(depths.max()*1e6),2),
                                 round(float(-depths.min()*1e6),2)])
    with (HERE/'sensitivity.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(['terminal_mask_minutes','duration_minutes','search_centers','max_dip_ppm','max_brightening_ppm'])
        writer.writerows(control_rows)
    stats['input_bins'] = len(t)
    stats['first_bin_bjd_tdb'] = float(2400000.5+bjd[0])
    stats['last_bin_bjd_tdb'] = float(2400000.5+bjd[-1])
    stats['input_sha256'] = hashlib.sha256((HERE/'one_minute.csv').read_bytes()).hexdigest()
    (HERE/'summary.json').write_text(json.dumps(stats,indent=2)+'\n')
    scans = []
    for _, degree, window in models:
        filt, valid = filter_matrix(t, n, known, centers, degree, window)
        scans.append((filt @ (y - 1), valid))
    with (HERE/'scan.csv').open('w') as f:
        w = csv.writer(f)
        w.writerow(['bjd_tdb','primary_depth_ppm','quadratic_depth_ppm'])
        for ix, c in enumerate(centers):
            w.writerow([f'{2400000.5+bjd[0]+c/1440:.8f}']+
                       [round(v[ix]*1e6,2) if m[ix] else '' for v,m in scans])


if __name__=='__main__':
    main()
