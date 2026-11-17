"""Reproduce the original optical matching procedure and measure the TXS weekly bin.

Run with a scratch directory; the three optical inputs are from the authors' tagged
software and the public LAT curve is pinned by the accompanying compressed snapshot.
"""
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
from math import cos, sin, radians, degrees, asin, sqrt
from pathlib import Path
import shlex
import sys
from urllib.request import urlopen

from statistics import median

HERE = Path(__file__).resolve().parent
GIT = '91cab012adad9dbc34de6fa5d32a61f457193da3'
BASE = f'https://raw.githubusercontent.com/sjoertvv/HelloLancel/{GIT}/hellolancel/data/'
LCR = ('https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/queryDB.php?'
       'typeOfRequest=lightCurveData&source_name=4FGL%20J0509.4%2B0542&cadence=weekly'
       '&flux_type=photon&index_type=fixed&ts_min=4')
UTC = timezone.utc


def separation(x, y, a, b):
    h = sin(radians(b-y)/2)**2 + cos(radians(y))*cos(radians(b))*sin(radians(a-x)/2)**2
    return degrees(2*asin(sqrt(h)))


def table(raw):
    lines = [line for line in raw.decode().splitlines() if line and not line.startswith('#')]
    keys = lines[0].split()
    return [dict(zip(keys, shlex.split(s, comments=False))) for s in lines[1:]]


def write_csv(path, rows):
    with path.open('w', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=rows[0])
        w.writeheader()
        w.writerows(rows)


def main(cache):
    cache.mkdir(exist_ok=True, parents=True)
    receipts = []

    def acquire(part):
        dest = cache / Path(part).name
        if not dest.exists():
            dest.write_bytes(urlopen(BASE+part, timeout=45).read())
        raw = dest.read_bytes()
        receipts.append(dict(url=BASE+part, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
        return table(raw)
    flares = acquire('ZTF_NEOWISE/ZTF_neoWISE_flares_acflares.dat')
    ice = acquire('IC/IC_alerts_merged.txt')
    notices = acquire('IC/AMON_ICECUBE_GOLD_AND_BRONZE.dat')+acquire('IC/AMON_ICECUBE_EHE.dat')
    dates = {r['Date'] for r in notices if float(r['Signalness']) > 0}
    aliases = {'ZTF19aaejtoy': 'AT2019aalc', 'ZTF19aapreis': 'AT2019dsg', 'ZTF19aatubsj': 'AT2019fdr'}
    used = 0
    pairs = []
    for r in ice:
        date = datetime.strptime(r['Event'][2:8], '%y%m%d').replace(tzinfo=UTC)
        ra, rp, rm, dec, dp, dm = map(
            float, (r[k] for k in ['RA', 'RA_plus', 'RA_min', 'Dec', 'Dec_plus', 'Dec_min']))
        # Reproduce setup_ice.py's original area expression and date-at-midnight convention.
        area = degrees(sin(radians(dec+dm))-sin(radians(dec-dp)))*(rp+rm)
        if date.year >= 2021 or area >= 100 or date.strftime('%y/%m/%d') not in dates:
            continue
        used += 1
        for f in flares:
            delay = date.timestamp()/86400+2440587.5-float(f['flare_peak_jd'])
            a, b = float(f['ra']), float(f['dec'])
            if not (0 < delay < min(365, float(f['echo_duration_w1']), float(f['flare_postpeak_duration_day']))):
                continue
            if ra-rm < a < ra+rp and dec-dm < b < dec+dp:
                pairs.append(dict(source=aliases.get(f['name'], f['name']), ztf=f['name'], alert=r['Event'], ra=a, dec=b, alert_ra=ra, alert_dec=dec, separation_deg=round(
                    separation(ra, dec, a, b), 5), peak_mjd=round(float(f['flare_peak_jd'])-2400000.5, 6), delay_days=round(delay, 4), geometry='original authors GCN rectangle'))
    assert used == 43 and len(flares) == 63 and {p['source']
                                                 for p in pairs} == set(aliases.values()) and len(pairs) == 3
    write_csv(HERE/'historical_optical_pairs.csv', pairs)
    txs = gzip.decompress((HERE/'txs_weekly.json.gz').read_bytes())
    lc = json.loads(txs)
    receipts.append(dict(url=LCR, bytes=len(txs), sha256=hashlib.sha256(txs).hexdigest()))
    ts = dict(lc['ts'])
    converged = dict(lc['fit_convergence'])
    bounds = {t: (lo, hi) for t, lo, hi in lc['flux_error']}
    bins = []
    for met, flux in lc['flux']:
        if met not in bounds:
            continue
        lo, hi = bounds[met]
        sigma = (hi-lo)/2
        if converged.get(met) == 0 and ts[met] >= 9 and sigma > 0 and flux/sigma >= 3 and flux > 0:
            bins.append((met, flux, sigma))
    baseline = median(f for _, f, _ in bins)
    # Mission elapsed time is TT elapsed seconds since 2001-01-01 UTC. There are five leap seconds by 2017.
    event = datetime(2017, 9, 22, 20, 54, 30, 430000, tzinfo=UTC)
    met0 = datetime(2001, 1, 1, tzinfo=UTC).timestamp()
    ar, ad = 77.43, 5.72  # GCN 21916
    sr, sd = (5+9/60+25.9644/3600)*15, 5+41/60+35.333/3600  # radio coordinates, ATel 11419
    met, flux, sigma = min(bins, key=lambda row: abs(met0+row[0]-5-event.timestamp()))
    row = dict(source='TXS 0506+056', alert='IC170922A', ra=sr, dec=sd, alert_ra=ar, alert_dec=ad, separation_deg=round(separation(ar, ad, sr, sd), 5), within_gcn_extrema=ar-.80 <= sr <= ar+1.30 and ad-.40 <= sd <= ad+.70, weekly_midpoint=datetime.fromtimestamp(
        met0+met-5, UTC).isoformat(), delay_days=round((event.timestamp()-met0-met+5)/86400, 5), flux_ph_cm2_s=flux, flux_err=sigma, ts=ts[met], median_good_flux=baseline, flux_over_median=round(flux/baseline, 4), passes_bin_threshold=flux >= 3*baseline)
    write_csv(HERE/'historical_txs.csv', [row])
    (HERE/'receipts.json').write_text(json.dumps(receipts, indent=2)+'\n')
    print(f'{used} alerts x {len(flares)} flares: {len(pairs)} historical optical pairs. TXS weekly bin ratio {flux/baseline:.2f}.')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
