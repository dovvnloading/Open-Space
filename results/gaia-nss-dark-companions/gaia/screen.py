"""Screen the full ESA Gaia DR3 two-body release using the public IsocLum primary masses."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import shutil
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

BASE = 'https://cdn.gea.esac.esa.int/Gaia/gdr3/'
INPUTS = {
    'binary': ('Performance_verification/binary_masses/BinaryMasses-001.csv.gz', 'ec2b0ec0ca54d7ee67a232cf07d4befb'),
    'nss': ('Non-single_stars/nss_two_body_orbit/NssTwoBodyOrbit_1.csv.gz', '6ca835765e60066bbbeb353a83b749c7'),
}
BAD_ASTROMETRY = {'4698424845771339520', '5765846127180770432', '522135261462534528', '1712614124767394816'}
SB = {'SB1', 'SB1C'}
# Derived from IAU nominal solar GM (1.3271244e20 m^3/s^2).
SB_FACTOR = 86400 * 1000**3 / (2 * math.pi * 1.3271244e20)
FMT = ('source_id nss_solution_type period period_error eccentricity mass_basis m1 m1_lower m1_upper '
       'f_msun f_rv_msun m2_screen_msun parallax parallax_error significance goodness_of_fit '
       'efficiency conf_spectro_period flags rv_n_good_obs_primary period_near_scan_or_year '
       'esa_withdrawn_astrometry has_sb2_or_eclipse mass_method').split()


def inputs(cache):
    metadata = {}
    for name, (remote, expected) in INPUTS.items():
        dest = cache / (name + '.csv.gz')
        if not dest.exists():
            temp = dest.with_suffix('.part')
            with urllib.request.urlopen(BASE + remote, timeout=180) as src, temp.open('wb') as out:
                shutil.copyfileobj(src, out)
            temp.replace(dest)
        with dest.open('rb') as f:
            actual = hashlib.file_digest(f, 'md5').hexdigest()
        if actual != expected:
            raise ValueError(f'{dest}: MD5 differs from ESA _MD5SUM.txt')
        with dest.open('rb') as f:
            metadata[name] = dict(url=BASE + remote, bytes=dest.stat().st_size, md5=actual,
                                  sha256=hashlib.file_digest(f, 'sha256').hexdigest())
    return metadata


def records(path):
    with gzip.open(path, 'rt', newline='') as file:
        yield from csv.DictReader(line for line in file if not line.startswith('#'))


def number(row, field):
    v = row.get(field)
    return math.nan if v in (None, '', 'null') else float(v)


def mass_ratio_root(f, m1):
    # For positive masses M2 is in [f, f+2*M1] and the equation increases monotonically.
    lo, hi = f, f + 2 * m1
    for _ in range(48):
        m = (lo + hi) / 2
        if m**3 < f * (m1 + m)**2:
            lo = m
        else:
            hi = m
    return (lo + hi) / 2


def mass_functions(r):
    p, e = number(r, 'period'), number(r, 'eccentricity')
    typ = r['nss_solution_type']
    if typ in SB:
        if typ == 'SB1C' and math.isnan(e):
            e = 0
        k = number(r, 'semi_amplitude_primary')
        f = SB_FACTOR * p * k**3 * (1 - e*e)**1.5 if k > 0 and 0 <= e < 1 and p > 0 else math.nan
        return f, f, 'SB1_edge_on'
    a,b,g,f = (number(r, name + '_thiele_innes') for name in ('a','b','g','f'))
    parallax = number(r,'parallax')
    if not (p > 0 and parallax > 0 and all(map(math.isfinite,[a,b,f,g]))):
        return math.nan, math.nan, 'astrometry_dark'
    # Larger singular value of the 2 by 2 Thiele-Innes projection matrix (mas).
    a0 = 0.5*(math.hypot(a+g,b-f) + math.hypot(a-g,b+f))
    f_phot = (a0/parallax)**3 / (p/365.25)**2
    c, h = number(r,'c_thiele_innes'), number(r,'h_thiele_innes')
    f_rv = math.hypot(c,h)**3 / (p/365.25)**2 if math.isfinite(c+h) else math.nan
    return f_phot, f_rv, 'astrometry_dark'


def run(cache, output):
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    metadata = inputs(cache)
    masses = defaultdict(dict)
    for r in records(cache/'binary.csv.gz'):
        if r['m1_ref'] != 'IsocLum':
            continue
        method = r['combination_method']
        masses[r['source_id']][method] = (r['m1'],r['m1_lower'],r['m1_upper'],method)
    counts, valid, with_mass = Counter(), Counter(), Counter()
    sb2_eclipse = set()
    kept, valid_ids, iso_ids = [], set(), set()
    for r in records(cache/'nss.csv.gz'):
        typ, sid = r['nss_solution_type'], r['source_id']
        counts[typ] += 1
        if typ.startswith(('SB2','Eclipsing')):
            sb2_eclipse.add(sid)
        if typ not in SB and typ != 'AstroSpectroSB1' and not typ.startswith('Orbital'):
            continue
        f, f_rv, basis = mass_functions(r)
        if not math.isfinite(f) or f <= 0:
            continue
        valid[typ] += 1
        valid_ids.add(sid)
        if typ in SB:
            methods = ('SB1+M1','Orbital+SB1+M1','Eclipsing+SB1+M1')
        elif typ == 'AstroSpectroSB1':
            methods = ('AstroSpectroSB1+M1',)
        else:
            methods = ('Orbital+M1','Orbital+SB1+M1')
        available = masses.get(sid, {})
        mass = next((available[m] for m in methods if m in available),None)
        if mass and float(mass[0]) > 0:
            with_mass[typ] += 1
            iso_ids.add(sid)
            m1, low, high, method = mass
            m1n = float(m1)
            # M2 exceeds t iff f exceeds t^3/(t + M1)^2.
            if f <= 1.4**3/(1.4+m1n)**2:
                continue
            m2 = mass_ratio_root(f, m1n)
        else:
            if f <= 1.4:
                continue
            m2 = f  # mass-independent floor within the chosen model
            m1, low, high, method = '', '', '', 'M1=0_floor'
        row = {k:r.get(k,'') for k in FMT}
        p = float(r['period'])
        row.update(mass_basis=basis,m1=m1,m1_lower=low,m1_upper=high,mass_method=method,
                   m2_screen_msun=m2,f_msun=f,f_rv_msun=f_rv if math.isfinite(f_rv) else '',
                   period_near_scan_or_year=int(any(abs(p/t - 1) < .03 for t in [62.97,365.25])),
                   esa_withdrawn_astrometry=int(sid in BAD_ASTROMETRY and basis == 'astrometry_dark'))
        kept.append(row)
    for row in kept:
        row['has_sb2_or_eclipse'] = int(row['source_id'] in sb2_eclipse)
    kept.sort(key=lambda r:(int(r['source_id']), r['nss_solution_type']))
    with (output/'mass_screen.csv').open('w',newline='') as file:
        writer=csv.DictWriter(file,FMT,lineterminator='\n')
        writer.writeheader()
        writer.writerows(kept)
    cuts = {}
    for t in [1.4,3]:
        hits = [r for r in kept if r['m2_screen_msun'] > t]
        cuts[str(t)] = dict(rows=len(hits),sources=len({r['source_id'] for r in hits}),
                            by_type=dict(Counter(r['nss_solution_type'] for r in hits)),
                            isoclum_rows=sum(bool(r['m1']) for r in hits),
                            floor_rows=sum(not r['m1'] for r in hits))
    result = dict(inputs=metadata,all_rows_by_type=counts,valid_rows_by_type=valid,
                  isoclum_rows_by_type=with_mass,valid_unique_sources=len(valid_ids),
                  isoclum_unique_sources=len(iso_ids),thresholds_msun=cuts)
    (output/'depth.json').write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', required=True,type=Path)
    parser.add_argument('--out', default=Path(__file__).parent,type=Path)
    args = parser.parse_args()
    run(args.cache,args.out)
