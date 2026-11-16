"""Retrieve the 6 closest type-1 SWORD v16 Zambezi named reaches per station.

Usage: python acquire.py OUTPUT_DIR. Standard library, Python 3.10+.
Archive byte-exact Hydrocron CSV payloads, FTS GeoJSON and CMR granule metadata.
"""
import csv
import hashlib
import io
import json
import math
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

FTS = 'https://fts.podaac.earthdata.nasa.gov/v1/rivers/Zambezi'
CMR = 'https://cmr.earthdata.nasa.gov/search/granules.json'
HYDROCRON = 'https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries'
GAUGES = {'chavuma': (22.68, -13.08), 'nanas_farm': (25.65139, -17.82028)}
FIELDS = ('reach_id,time_str,p_lat,p_lon,wse,wse_u,slope,slope_u,slope2,slope2_u,width,width_u,'
          'reach_q,reach_q_b,partial_f,n_good_nod,obs_frac_n,loc_offset,xtrk_dist,node_dist,'
          'dark_frac,xovr_cal_q,p_n_nodes,p_dist_out,p_length,cycle_id,pass_id,range_start_time,'
          'range_end_time,crid,granuleUR,sword_version,collection_shortname,sos_consensus_q,'
          'sos_hivdi_q,sos_sad_q,sos_metroman_q,sos_momma_q,sos_sic4dvar_q,sos_lakeflow_q')


def fetch(url):
    for attempt in range(3):
        try:
            with urlopen(url, timeout=60) as r:
                return r.read()
        except OSError:
            if attempt == 2:
                raise
            time.sleep(2**attempt)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def distance(row, gauge):
    """Local equirectangular distance (km) to nearest SWORD polyline vertex."""
    lon, lat = gauge
    xscale = 111.32 * math.cos(math.radians(lat))
    return min(math.hypot(xscale * (x - lon), 111.32 * (y - lat))
               for x, y in row['geojson']['coordinates'])


def main(out):
    out.mkdir(parents=True, exist_ok=True)
    manifest = {'retrieval_utc': datetime.now(timezone.utc).isoformat(), 'stations_lon_lat': GAUGES,
                'cmr': {}, 'fts': [], 'hydrocron': []}
    url = CMR + '?' + urlencode(dict(collection_concept_id='C3905028734-POCLOUD', page_size=200))
    data = fetch(url)
    granule = next(x for x in json.loads(data)['feed']['entry'] if x['title'].startswith('af_'))
    archive = (json.dumps(granule, sort_keys=True, indent=2) + '\n').encode()
    (out / 'af_cmr_granule.json').write_bytes(archive)
    manifest['cmr'] = {'url': url, 'response_sha256': sha(data), 'granule_file_sha256': sha(archive),
                       'concept_id': granule['id']}

    def fts_page(i):
        url = FTS + '?' + urlencode(dict(nodes='false', page_size=20, page_number=i))
        data = fetch(url)
        return json.loads(data), {'url': url, 'response_sha256': sha(data)}

    first, record = fts_page(1)
    rows = first['results']
    manifest['fts'].append(record)
    with ThreadPoolExecutor(max_workers=4) as workers:
        for page, record in workers.map(fts_page, range(2, math.ceil(first['hits'] / 20) + 1)):
            rows.extend(page['results'])
            manifest['fts'].append(record)
    features = [r for r in rows if r['shp_origin'].endswith('_v16.shp') and r['type'] == 1]
    selected = []
    for name, point in GAUGES.items():
        for rank, row in enumerate(sorted(features, key=lambda r: distance(r, point))[:6], 1):
            selected.append(dict(type='Feature', geometry=row['geojson'], properties={
                **{k: v for k, v in row.items() if k not in ('geojson', 'geometry')},
                'site': name, 'proximity_rank': rank, 'distance_km': distance(row, point)}))
    archive = (json.dumps(dict(type='FeatureCollection', features=selected), separators=(',', ':'))+'\n').encode()
    (out / 'sword_v16_selected.geojson').write_bytes(archive)
    manifest['selected_geometry_sha256'] = sha(archive)

    def hydro(feature):
        props = feature['properties']
        reach = props['reach_id']
        query = dict(feature='Reach', feature_id=reach, start_time='2023-03-29T00:00:00Z',
                     end_time='2025-05-03T23:59:59Z', output='csv',
                     collection_name='SWOT_L2_HR_RiverSP_2.0', fields=FIELDS)
        url = HYDROCRON + '?' + urlencode(query)
        data = fetch(url)
        response = json.loads(data)
        archive = response['results']['csv'].encode()
        rows = list(csv.DictReader(io.StringIO(archive.decode())))
        assert len(rows) == response['hits'] and all(r['reach_id'] == reach for r in rows)
        path = 'reach_' + reach + '.csv'
        (out / path).write_bytes(archive)
        return dict(url=url, site=props['site'], rank=props['proximity_rank'], reach_id=reach,
                    rows=len(rows), file=path, csv_sha256=sha(archive), response_sha256=sha(data))

    with ThreadPoolExecutor(max_workers=4) as workers:
        manifest['hydrocron'] = list(workers.map(hydro, selected))
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')


if __name__ == '__main__':
    main(Path(sys.argv[1]))
