"""Check the published ephemerides against a separately extracted TESS product."""
import csv
from pathlib import Path
import sys
from zipfile import ZipFile

import numpy as np

root = Path(__file__).resolve().parent
sys.path.insert(0, str(root.parent))
from measure import event_rows


def alignment_times(ident, jd):
    # Both archives used half-hour FFIs. Map to Pál timestamps only where the
    # nearest reported exposure is within half a cadence. No flux enters the join.
    with ZipFile(root.parent/'target-inputs.zip') as archive:
        with archive.open(f'lightcurves_spectra/{ident}.lc') as handle:
            pal = np.loadtxt(handle, usecols=1)
    idx = np.searchsorted(pal, jd)
    left = np.clip(idx-1,0,len(pal)-1)
    right = np.clip(idx,0,len(pal)-1)
    match = np.where(abs(pal[left]-jd)<abs(pal[right]-jd),left,right)
    ok = abs(pal[match]-jd)<15/1440
    return pal[match], ok


def check(a, jd, row, mode, offset):
    ident = int(row['mpc'])
    mag = a[:,3]
    t = jd-jd.min()
    x = t/np.ptp(t)-.5
    spin = 24/float(row['spin_h'])
    p = float(row['recurrence_h'])
    center = float(row['first_center_jd'])
    width = float(row['box_width_h'])/p
    for order in (6,14,28):
        angles = 2*np.pi*t[:,None]*spin*np.arange(1,order+1)
        model = np.column_stack([np.ones(len(t)),x,x*x,np.sin(angles),np.cos(angles)])
        coef = np.linalg.lstsq(model,mag,rcond=None)[0]
        resid = mag-model@coef
        events = event_rows(jd-center,resid,p,0,width)
        depths = {int(k):d for k,_,_,d in events}
        pairs = np.array([v-depths[k+1] for k,v in depths.items() if k%2==0 and k+1 in depths])
        boot = np.median(np.random.default_rng(1842).choice(events[:,3],(10000,len(events)),replace=True),axis=1)
        yield [ident,mode,offset,order,len(jd),float(np.median(a[:,4])),len(events),int(sum(events[:,3]>0)),
               float(np.median(events[:,3])),*np.quantile(boot,[.025,.975]),len(pairs),int(sum(pairs>0)),np.median(pairs)]


def main():
    source = np.loadtxt(root/'photometry.tsv')
    with (root.parent/'measurements.csv').open() as handle:
        targets = list(csv.DictReader(handle))[:2]
    rows = []
    for row in targets:
        ident = int(row['mpc'])
        data = source[source[:,0]==ident]
        jd = data[:,1]+2400000.5
        rows.extend(check(data,jd,row,'as_reported',0))
        matched, ok = alignment_times(ident,jd)
        rows.extend(check(data[ok],matched[ok],row,'matched_cadence',np.median((matched-jd)[ok])*1440))
    with (root/'test.csv').open('w',newline='') as handle:
        writer = csv.writer(handle,lineterminator='\n')
        writer.writerow('mpc times median_time_shift_minutes spin_harmonics n median_formal_err_mag covered_cycles positive median_depth_mag low high pairs even_deeper median_even_minus_odd_mag'.split())
        writer.writerows(rows)


if __name__ == '__main__':
    main()
