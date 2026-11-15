"""Fetch or replay the withheld Rubin positions against the older JPL target solutions."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
from statistics import median
from urllib.parse import urlencode
from urllib.request import Request, urlopen

TARGETS = [('2017 DO171', 1, 14), ('2022 GV6', 2, 18)]
ROOT = Path(__file__).resolve().parent


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def get_mpc(endpoint, payload):
    request = Request('https://data.minorplanetcenter.net/api/' + endpoint,
                      data=json.dumps(payload).encode(), headers={'Content-Type':'application/json'},
                      method='GET')
    return urlopen(request, timeout=90).read()


def epochs(rubin):
    return [datetime.fromisoformat(r['obstime']).astimezone(timezone.utc).timestamp()/86400 + 2440587.5
            for r in rubin]


def get_horizons(name, rubin):
    query = {'format':'text', 'COMMAND':"'DES="+name+";'", 'CENTER':"'X05'",
             'EPHEM_TYPE':"'OBSERVER'", 'QUANTITIES':"'1,36'", 'TLIST':'\n'.join(map(str, epochs(rubin))),
             'CSV_FORMAT':"'YES'", 'ANG_FORMAT':"'DEG'", 'TIME_DIGITS':"'FRACSEC'",
             'EXTRA_PREC':"'YES'", 'REF_SYSTEM':"'ICRF'", 'APPARENT':"'AIRLESS'"}
    return urlopen('https://ssd.jpl.nasa.gov/api/horizons.api?'+urlencode(query), timeout=90).read()


def residuals(table):
    results = {}
    for name, _, _ in TARGETS:
        subset = [row for row in table if row['target']==name]
        separations = []
        by_night = {}
        east_values, north_values = [], []
        for row in subset:
            east = (float(row['obs_ra_deg'])-float(row['jpl_ra_deg'])) * 3600 * math.cos(math.radians(float(row['jpl_dec_deg'])))
            north = (float(row['obs_dec_deg'])-float(row['jpl_dec_deg'])) * 3600
            sep = math.hypot(east, north)
            east_values.append(east)
            north_values.append(north)
            separations.append(sep)
            by_night.setdefault(row['utc'][:10], []).append(sep)
        results[name] = {'count':len(subset), 'nights':len(by_night),
                         'median_east_arcsec':round(median(east_values), 4),
                         'median_north_arcsec':round(median(north_values), 4),
                         'min_median_max_radial_arcsec':[round(f(separations), 4) for f in (min, median, max)],
                         'night_median_radial_arcsec':{k:round(median(v), 4) for k,v in sorted(by_night.items())}}
    return results


def freeze(directory, offline=None):
    table, provenance = [], {}
    directory.mkdir(parents=True,exist_ok=True)
    for name, solution, old_count in TARGETS:
        base = name.replace(' ', '-')
        if offline:
            observations = (offline/(base+'.dat')).read_bytes()
            orbit = (offline/(base+'-orbit.json')).read_bytes()
        else:
            observations = get_mpc('get-obs', {'desigs':[name], 'output_format':['OBS80','ADES_DF']})
            orbit = get_mpc('get-orb', {'desig':name})
        all_obs = json.loads(observations)[0]['ADES_DF']
        rubin = sorted([r for r in all_obs if r['stn']=='X05' and r['obstime'].startswith('2026') and not r['deprecated']],
                       key=lambda row: row['obstime'])
        fit = json.loads(orbit)[0]['mpc_orb'][0]
        assert fit['orbit_fit_statistics']['nobs_total_sel'] == len(all_obs) == old_count+len(rubin)
        raw_h = (offline/(base+'-horizons.txt')).read_bytes() if offline else get_horizons(name, rubin)
        horizon = raw_h.decode()
        if f'{{source: JPL#{solution}}}' not in horizon:
            raise ValueError(f'{name}: JPL solution superseded, use frozen CSV for replay')
        assert 'Simonyi Survey Telescope, Rubin Obs.' in horizon
        lines = list(csv.reader(horizon.split('$$SOE')[1].split('$$EOE')[0].strip().splitlines()))
        assert len(lines)==len(rubin)
        for obs, fields in zip(rubin, lines):
            api_time = datetime.strptime(fields[0].strip(), '%Y-%b-%d %H:%M:%S.%f').replace(tzinfo=timezone.utc)
            assert abs((datetime.fromisoformat(obs['obstime'])-api_time).total_seconds()) < .002
            table.append(dict(target=name, utc=obs['obstime'], station=obs['stn'], astcat=obs['astcat'],
                              obs_ra_deg=obs['ra'], obs_dec_deg=obs['dec'], obs_rmsra_arcsec=obs['rmsra'],
                              obs_rmsdec_arcsec=obs['rmsdec'], jpl_ra_deg=fields[3].strip(),
                              jpl_dec_deg=fields[4].strip(), jpl_3sigma_ra_arcsec=fields[5].strip(),
                              jpl_3sigma_dec_arcsec=fields[6].strip()))
        provenance[name] = {'sha256_mpc_observations_response':digest(observations),
                            'sha256_mpc_orbit_response':digest(orbit),
                            'sha256_horizons_response':digest(raw_h),
                            'jpl_orbit_id':solution,'old_count':old_count,'latest_observation':rubin[-1]['obstime'],
                            'software_data':fit['software_data'],'epoch_data':fit['epoch_data'],
                            'reference_system':fit['system_data'],'used_count':fit['orbit_fit_statistics']['nobs_total_sel'],
                            'new_2026_x05_records':len(rubin)}
    with (directory/'x05-astrometry.csv').open('w') as f:
        writer = csv.DictWriter(f, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)
    for filename, obj in [('mpc-receipts.json',provenance),('jpl-holdout.json',residuals(table))]:
        (directory/filename).write_text(json.dumps(obj,indent=2)+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fetch',type=Path,help='write fresh evidence to this directory')
    parser.add_argument('--raw-cache',type=Path,help='replay the source receipts in a local cache')
    args=parser.parse_args()
    if args.fetch:
        freeze(args.fetch,args.raw_cache)
    else:
        with (ROOT/'x05-astrometry.csv').open() as f:
            print(json.dumps(residuals(list(csv.DictReader(f))),indent=2))
