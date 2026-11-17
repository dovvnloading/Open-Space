"""Screen published secondary trajectories against public DR1 long VIS footprints."""
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
from astropy.time import Time
from matplotlib.path import Path as Polygon
from scipy.spatial import cKDTree
from candidates import HERE, act_rows, act_sky, circular_family, other_sky

QUERY = ("SELECT obs_id,t_min,t_max,s_region,s_ra,s_dec FROM ivoa.obscore "
         "WHERE obs_collection='dr1' AND dataproduct_type='image' "
         "AND t_exptime>560 AND t_exptime<561")
URL = "https://eas.esac.esa.int/tap-server/tap/sync?" + urlencode(
    dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=QUERY, MAXREC=150000))


def unit(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(d)*np.cos(a), np.cos(d)*np.sin(a), np.sin(d)], axis=-1)


def tangent(ra, dec, center):
    a, d = np.radians(center)
    xyz = unit(ra, dec)
    east = [-np.sin(a), np.cos(a), 0.]
    north = [-np.sin(d)*np.cos(a), -np.sin(d)*np.sin(a), np.cos(d)]
    denominator = xyz @ unit(*center)
    return np.column_stack([xyz @ east, xyz @ north]) / denominator[:, None]


def main(cache):
    if not cache.exists():
        with urlopen(URL, timeout=180) as r:
            cache.write_bytes(r.read())
    raw = cache.read_bytes()
    products = list(csv.DictReader(io.StringIO(raw.decode())))
    assert 0 < len(products) < 150000 and 's_region' in products[0]
    bounds = np.array([[float(r['s_ra']), float(r['s_dec'])] for r in products])
    centers = unit(bounds[:,0], bounds[:,1])
    mid = np.array([(float(r['t_min'])+float(r['t_max']))/2 for r in products])
    t = Time(mid, format="mjd", scale="utc").tt.jd-2451545.
    days = np.floor(t).astype(int)
    groups = [(day, np.flatnonzero(days == day)) for day in np.unique(days)]
    cases = {**{k: (act_sky,row) for k,row in act_rows().items()},
             **{k: (other_sky,circular_family(k)) for k in ('RR','Chen1','Chen2')}}
    hits, summary = [], []
    for name, (positions, args) in cases.items():
        poss = []
        for day, indices in groups:
            ra, dec = positions(args, day+.5)
            tree = cKDTree(unit(ra, dec))
            sep = tree.query(centers[indices])[0]
            poss.extend(indices[sep < 2*np.sin(np.radians(1.6)/2)])
        covered = set()
        count = 0
        for i in poss:
            r = products[i]
            poly = np.array(r['s_region'].split()[1:], dtype=float).reshape(-1,2)
            ra, dec = positions(args,t[i])
            polygon = Polygon(tangent(*poly.T, bounds[i]))
            inside = np.flatnonzero(polygon.contains_points(tangent(ra,dec,bounds[i]))).tolist()
            if inside:
                count += 1
                covered.update(inside)
                hits.append(dict(target=name, product_row=int(i), obs_id=r['obs_id'],
                                 t_min=r['t_min'], t_max=r['t_max'], covered_track_ids=inside,
                                 s_region=r['s_region']))
        summary.append(dict(target=name, tracks=len(ra), preselected_products=len(poss),
                            matched_products=count, distinct_track_ids=len(covered)))
    with (HERE/'euclid_footprints.csv').open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
    payload=''.join(json.dumps(h,separators=(',',':'))+'\n' for h in hits)
    (HERE/'euclid_intersections.jsonl.gz').write_bytes(gzip.compress(payload.encode(),mtime=0))
    receipt=dict(url=URL,query=QUERY,archive_rows=len(products),sha256=hashlib.sha256(raw).hexdigest(),
                 bytes=len(raw),min_mjd_utc=min(mid),max_mjd_utc=max(mid),
                 tangent_point_approximation='geocentre; exact Euclid parallax and detector gaps not applied')
    (HERE/'euclid_footprint_query.json').write_text(json.dumps(receipt,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('cache',type=Path,help='CSV response cache')
    main(p.parse_args().cache)
