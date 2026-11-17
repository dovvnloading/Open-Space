"""Compare named GCN triggers to every measured optical pair (requires astropy-healpix)."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy_healpix import HEALPix

HERE = Path(__file__).resolve().parent


def read_pairs(path, tai=False):
    with path.open() as handle:
        for r in csv.DictReader(handle):
            # Broker Rubin MJD is TAI. The offset is 37 s throughout this 2026 window.
            offset = 37 / 86400 if tai else 0.0
            yield (r['oid'], float(r.get('ra', r.get('meanra'))), float(r.get('dec', r.get('meandec'))),
                   float(r.get('t1', r.get('mjd1'))) - offset,
                   float(r.get('t2', r.get('mjd2'))) - offset)


def sky_map(filename):
    with fits.open(filename) as hs:
        header = hs[1].header
        pixels = HEALPix(nside=header['NSIDE'], order=header['ORDERING'].lower())
        prob = hs[1].data['PROBABILITY'].flatten().astype(float)
    if np.min(prob) < 0 or not .98 < prob.sum() < 1.02:
        raise ValueError(f'Unexpected probability map in {filename}')
    # Cumulative area ranked by posterior density. Equal-area fixed-order HEALPix maps.
    # Handle equal-probability pixels conservatively (include all ties in the cumulative mass).
    ordered = np.sort(prob)[::-1]
    csum = np.cumsum(ordered) / ordered.sum()
    return pixels, prob, ordered, csum


def run(rubin, ztf, cache):
    cache.mkdir(exist_ok=True, parents=True)
    with (HERE / 'trigger_localizations.csv').open() as handle:
        catalog = list(csv.DictReader(handle))
    pairs = list(read_pairs(rubin, tai=True)) + list(read_pairs(ztf))
    out = []
    map_receipts = []
    for tr in catalog:
        trigger_mjd = Time(tr['utc'], format='isot', scale='utc').mjd
        tloc = SkyCoord(float(tr['ra_deg']) * u.deg, float(tr['dec_deg']) * u.deg)
        if tr['map_id']:
            name = tr['map_id']
            url = f'https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/triggers/2026/bn{name}/quicklook/glg_healpix_all_bn{name}.fit'
            path = cache / (tr['gcn'] + '.fit')
            if not path.exists():
                with urlopen(url,timeout=60) as r:
                    path.write_bytes(r.read())
            raw = path.read_bytes()
            map_receipts.append(dict(gcn=tr['gcn'],url=url,size=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
            pixels, prob, ordered, csum = sky_map(path)
        by_object = {}
        for oid, ra, dec, t1, t2 in pairs:
            if not t1 - 3 <= trigger_mjd <= t2 + 1:
                continue
            sc = SkyCoord(ra * u.deg, dec * u.deg)
            dist = tloc.separation(sc).deg
            if tr['map_id']:
                idx = pixels.lonlat_to_healpix(sc.ra, sc.dec)
                rank = np.searchsorted(-ordered, -prob[idx], side='right') - 1
                credible = float(csum[rank])
                hit = credible <= .90
            else:
                credible = None
                hit = dist <= float(tr['radius_deg'])
            # Keep the temporally compatible pair closest to the trigger in its first epoch.
            if oid not in by_object or abs(t1 - trigger_mjd) < abs(by_object[oid][4]):
                by_object[oid] = [oid, tr['event'], tr['gcn'], tr['method'], t1 - trigger_mjd,
                                  round(dist, 6), '' if credible is None else round(credible, 8), int(hit)]
        out.extend(by_object.values())
    with (HERE / 'sky_time_matches.csv').open('w') as handle:
        w = csv.writer(handle,lineterminator='\n')
        w.writerow(['oid','event','gcn','localization','optical_t1_minus_trigger_d','separation_deg','gbm_credible_level','inside_region'])
        w.writerows(out)
    (HERE / 'map_receipts.json').write_text(json.dumps(map_receipts,indent=2)+'\n')
    hits = [row for row in out if row[-1]]
    print(f'{len(set(p[0] for p in pairs))} optical objects, {len(catalog)} localizations, {len(out)} temporally compatible object-localization comparisons, {len(hits)} overlaps')
    for row in hits:
        print(row)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('rubin', type=Path)
    p.add_argument('ztf', type=Path)
    p.add_argument('--cache', type=Path, default=Path('map-cache'))
    a = p.parse_args()
    run(a.rubin, a.ztf, a.cache)
