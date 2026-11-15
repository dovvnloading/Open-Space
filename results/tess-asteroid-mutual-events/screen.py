"""Blind search for recurring residual modulations in the TSSYS-DR1 archive.

Usage: OPENBLAS_NUM_THREADS=1 python screen.py /path/to/release /path/to/output
Run before the candidate-specific follow-up in measure.py.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.linalg import lstsq
from scipy.optimize import minimize_scalar

COUNTERS = ('release', 'numbered_and_280_rows', 'usable_260_rows', 'baseline_and_spin',
            'adequate_amplitude', 'retained')
HARMONICS = (.25, .333, .5, .667, 1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)


def harmonics(t, f, order):
    phase = 2*np.pi*t[:, None]*f*np.arange(1, order+1)
    return np.column_stack([np.sin(phase), np.cos(phase)])


def spin_residual(t, y, f):
    design = np.column_stack([np.ones(len(t)), t, t*t, harmonics(t, f, 6)])
    coefs = lstsq(design, y, lapack_driver='gelsy')[0]
    return y - design @ coefs


def sine(t, y, freq):
    coefs = np.linalg.lstsq(harmonics(t, freq, 1), y-np.mean(y), rcond=None)[0]
    return np.linalg.norm(coefs), np.arctan2(*coefs)/(2*np.pi)%1


def evaluate(root, row, counts):
    ident = row[0]
    counts['release'] += 1
    if not ident.isdigit() or int(row[5]) < 280:
        return
    counts['numbered_and_280_rows'] += 1
    with (root/'lightcurves_spectra'/f'{ident}.lc').open() as handle:
        data = np.loadtxt((l for l in handle if len(l.split()) == 13),
                          usecols=(1, 2, 3, 4, 5, 6, 7, 9, 10, 11))
    ok = np.all(np.isfinite(data),axis=1) & (data[:,7]==0) & (data[:,4]>0) & (data[:,4]<.15)
    a = data[ok]
    if len(a) < 260:
        return
    counts['usable_260_rows'] += 1
    a = a[a[:,4] < max(.02, 3*np.median(a[:,4]))]
    t = a[:,0] - a[0,0]
    span = np.ptp(t)
    f0 = float(row[1])
    if span < 5 or not .5 < f0 < 11:
        return
    counts['baseline_and_spin'] += 1
    y = a[:,3] - 5*np.log10(a[:,8]*a[:,9])
    spin = minimize_scalar(lambda f: np.sum(spin_residual(t,y,f)**2),
                           bounds=(f0-.45/span, f0+.45/span),method='bounded').x
    res = spin_residual(t,y,spin)
    sigma = np.median(abs(res-np.median(res)))/.6745
    good = abs(res-np.median(res)) < max(6*sigma,.1)
    t = t[good]
    res = np.clip(spin_residual(t,y[good],spin),-4*sigma,4*sigma)
    # The TSSYS-DR1 cadence is half an hour. Gaps contribute zero weighted residual.
    n = np.rint((t-t[0])*48).astype(int)
    signal = np.zeros(max(n)+1)
    signal[n] = res-np.mean(res)
    spectrum = abs(np.fft.rfft(signal,n=8192))**2
    freq = np.fft.rfftfreq(8192,1/48)
    valid = (freq>max(.48,7/span)) & (freq<4)
    for h in HARMONICS:
        for alias in (spin*h,abs(48-spin*h),96-spin*h):
            valid &= abs(freq-alias) > 1.4/span
    if not np.any(valid):
        return
    peak = freq[valid][np.argmax(spectrum[valid])]
    amplitude, _ = sine(t,res,peak)
    if amplitude < .012:
        return
    counts['adequate_amplitude'] += 1
    first = t < np.median(t)
    a1, p1 = sine(t[first],res[first],peak)
    a2, p2 = sine(t[~first],res[~first],peak)
    score = amplitude/sigma*np.sqrt(len(t)/2)
    phase_delta = abs((p1-p2+.5)%1-.5)
    ratio = min(a1,a2)/amplitude
    if score < 9 or phase_delta > .15 or ratio < .45:
        return
    counts['retained'] += 1
    return (ident,round(24/spin,4),round(24/peak,3),round(2*amplitude,4),
            round(sigma,4),round(score,2),round(phase_delta,3),row[2],row[6],row[9],
            len(t),round(span,2),round(ratio,2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('release',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    counts = dict.fromkeys(COUNTERS,0)
    with (args.release/'release.merge').open() as handle:
        # One record per identifier in the catalog; preserve release order for tied scores.
        records = {line.split()[0]:line.split() for line in handle}
    results = [r for row in records.values() if (r:=evaluate(args.release,row,counts)) is not None]
    with (args.output/'screen.csv').open('w',newline='') as handle:
        writer = csv.writer(handle,lineterminator='\n')
        writer.writerow('mpc spin_h secondary_h pk_to_pk_mag mad_mag ranking_score half_phase_delta release_type sector tess_mag n span_d min_half_amp_ratio'.split())
        writer.writerows(sorted(results,key=lambda r:-r[5]))
    (args.output/'coverage.json').write_text(json.dumps(counts,indent=2)+'\n')


if __name__ == '__main__':
    main()
