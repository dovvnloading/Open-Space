"""Recompute Gaia-to-WISE offsets from the printed (rounded) Ren+ coordinates."""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent

def angle(s):
    h,m,sec = [float(x) for x in s.split(":")]
    return math.copysign(abs(h) + m/60 + sec/3600, -1 if s.startswith("-") else 1)

def sep(ra, dec, ra0, dec0):
    d1, d2 = math.radians(dec), math.radians(dec0)
    h = math.sin((d1-d2)/2)**2+math.cos(d1)*math.cos(d2)*math.sin(math.radians(ra-ra0)/2)**2
    return 7200*math.degrees(math.asin(math.sqrt(h)))

def main():
    stars = {r['label']:r for r in csv.DictReader((ROOT/'candidates.csv').open())}
    results = []
    for r in csv.DictReader((ROOT/'ren_centroids.csv').open()):
        s = stars[r['label']]
        row = dict(label=r['label'],band=r['band'],published_offset_arcsec=r['published_offset_arcsec'])
        for extra in (0,-16):
            dt = float(r['epoch'])+extra-float(s['ref_epoch'])
            dec=float(s['dec'])
            ra = float(s['ra'])+dt*float(s['pmra'])/(3.6e6*math.cos(math.radians(dec)))
            d = sep(15*angle(r['ra_hms']),angle(r['dec_dms']),ra,dec+dt*float(s['pmdec'])/3.6e6)
            row['direct_offset_arcsec' if extra==0 else 'earlier_epoch_diagnostic_arcsec'] = round(d,4)
        results.append(row)
    with (ROOT/'ren_offsets_recomputed.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(results[0]),lineterminator='\n')
        w.writeheader()
        w.writerows(results)

if __name__=='__main__':
    main()
