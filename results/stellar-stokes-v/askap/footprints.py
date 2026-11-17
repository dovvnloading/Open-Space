"""Inventory released, restored ASKAP Stokes V fields and match their WCS boundaries to hosts."""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from match_racs import frame, ARCSEC

HERE = Path(__file__).resolve().parent
SERVICE = 'https://casda.csiro.au/casda_vo_tools/tap/sync'
CUTOFF = '2026-11-17T11:00:00'
COLUMNS = 'obs_id,obs_publisher_did,filename,s_ra,s_dec,s_region,t_min,t_max,em_min,em_max,s_resolution,quality_level,obs_release_date'
SELECTIONS = {'racs': "obs_collection='The Rapid ASKAP Continuum Survey'", 'vast': "obs_collection IN ('ASKAP Pilot Survey for VAST','VAST')"}


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def archive_query(cache, key, query):
    url = SERVICE + '?' + urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', MAXREC=70000, QUERY=query))
    path = cache / key
    if not path.exists():
        with urlopen(url, timeout=240) as response:
            path.write_bytes(response.read())
    raw = path.read_bytes()
    records = list(csv.DictReader(io.StringIO(raw.decode())))
    if not records or raw.lstrip().startswith(b'<') or len(records) >= 70000:
        raise ValueError(f'Invalid or truncated {key} archive response')
    return records, dict(url=url, bytes=len(raw), sha256=digest(raw), rows=len(records))


def write_gz(name, rows):
    out = io.StringIO(newline='')
    writer = csv.DictWriter(out, fieldnames=rows[0].keys(), lineterminator='\n')
    writer.writeheader()
    writer.writerows(rows)
    raw = out.getvalue().encode()
    (HERE / name).write_bytes(gzip.compress(raw, mtime=0))
    return dict(rows=len(rows), csv_sha256=digest(raw))


def main(cache, hosts):
    cache.mkdir(parents=True, exist_ok=True)
    schedule, times_receipt = archive_query(cache, 'obs_times', 'SELECT sbid, obs_start_mjd, obs_end_mjd FROM casda.observation ORDER BY sbid')
    times = {r['sbid']: r for r in schedule}
    targets = pd.read_csv(hosts, dtype={'gaia_dr3_id': str})
    positions, east, north = frame(targets.ra_deg, targets.dec_deg)
    motion_known = targets[['ref_epoch','pmra_masyr','pmdec_masyr']].notna().all(axis=1).to_numpy()
    motion = np.where(motion_known[:,None], (targets.pmra_masyr.to_numpy()[:,None] * east + targets.pmdec_masyr.to_numpy()[:,None] * north) * ARCSEC / 1000, 0)
    epoch = targets.ref_epoch.fillna(2023).to_numpy()
    central = positions + (2023-epoch)[:,None]*motion
    central /= np.linalg.norm(central,axis=1)[:,None]
    tree = cKDTree(central)
    speed = np.linalg.norm(motion,axis=1).max()
    fields, visits, receipts = [], [], {'schedule': times_receipt}
    for survey, selection in SELECTIONS.items():
        query = (f'SELECT {COLUMNS} FROM ivoa.obscore WHERE {selection} AND pol_states=\'/V/\' '
                 f"AND filename LIKE 'image.v.%taylor.0.restored%.fits' AND obs_release_date < '{CUTOFF}' ORDER BY obs_publisher_did")
        records, receipts[survey] = archive_query(cache, 'v_'+survey+'2', query)
        for row in records:
            sbid = row['obs_id'].removeprefix('ASKAP-')
            timesrc = 'ObsCore'
            tmin,tmax = row['t_min'],row['t_max']
            if not(tmin and tmax):
                tmin,tmax = times[sbid]['obs_start_mjd'],times[sbid]['obs_end_mjd']
                timesrc = 'scheduling block'
            field = dict(survey=survey, **row, used_t_min=tmin, used_t_max=tmax, time_source=timesrc)
            fields.append(field)
            if row['quality_level'] not in ('GOOD','UNCERTAIN') or not(tmin and tmax):
                continue
            year = 2000 + (0.5*(float(tmin)+float(tmax))-51544.5)/365.25
            region = row['s_region'].split()
            if region[:2] != ['POLYGON','ICRS'] or len(region) != 10:
                raise ValueError('Unexpected ObsCore geometry: '+row['s_region'])
            coords = np.array(region[2:],float).reshape(-1,2)
            verts = frame(*coords.T)[0]
            center = verts.mean(axis=0); center /= np.linalg.norm(center)
            radius = np.linalg.norm(verts-center,axis=1).max() + abs(year-2023)*speed + 1e-7
            ix = np.array(tree.query_ball_point(center,radius),dtype=int)
            at = positions[ix] + (year - epoch[ix])[:,None]*motion[ix]
            at /= np.linalg.norm(at,axis=1)[:,None]
            normals = np.cross(verts, np.roll(verts,-1,axis=0))
            normals *= np.sign(normals@center)[:,None]
            inside = np.all(at @ normals.T >= 0,axis=1)
            for j,xyz in zip(ix[inside],at[inside]):
                visits.append(dict(hostname=targets.hostname.iat[j],survey=survey,sbid=sbid,product=row['obs_publisher_did'],
                    midpoint_mjd=round(0.5*(float(tmin)+float(tmax)),6),pm_applied=bool(motion_known[j]),
                    ra_epoch_deg=round(np.degrees(np.arctan2(xyz[1],xyz[0]))%360,8),dec_epoch_deg=round(np.degrees(np.arcsin(xyz[2])),8),quality=row['quality_level']))
    measurements = pd.read_csv(HERE/'racs_matches.csv').query("sample == 'hosts'")
    joined = pd.DataFrame(visits)
    totals = joined.groupby(["hostname", "survey"]).sbid.nunique()
    good = joined[joined.quality == "GOOD"].groupby(["hostname", "survey"]).sbid.nunique()
    census = []
    for target in targets.hostname:
        row = dict(hostname=target)
        for survey in SELECTIONS:
            row[survey+'_distinct_sbids'] = int(totals.get((target,survey),0))
            row[survey+'_good_sbids'] = int(good.get((target,survey),0))
        det=measurements[measurements.target==target]
        row['racs_low2_catalogue_local_5rms_mjy_beam'] = '' if det.empty else round(5*det.v_rms_mjy_beam.min(),4)
        row['pixel_limit_mjy_beam'] = ''
        row['status'] = ('published_detection' if not det.empty else 'footprint_only' if row['racs_distinct_sbids'] or row['vast_distinct_sbids'] else 'outside_selected_polygons')
        census.append(row)
    outputs = {'fields':write_gz('v_image_fields.csv.gz',fields), 'visits':write_gz('host_geometric_visits.csv.gz', sorted(visits,key=lambda r:(r['hostname'], r['midpoint_mjd'],r['product']))), 'hosts':write_gz('host_coverage.csv.gz',census)}
    (HERE/'footprint_receipt.json').write_text(json.dumps(dict(cutoff_utc=CUTOFF+'Z',host_csv_sha256=digest(targets_raw(hosts)),inputs=receipts,outputs=outputs),indent=2)+'\n')


def targets_raw(path):
    return gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()


if __name__=='__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache',type=Path,required=True)
    ap.add_argument('--hosts',type=Path,default=HERE.parent/'samples/hosts.csv.gz')
    args=ap.parse_args();main(args.cache,args.hosts)
