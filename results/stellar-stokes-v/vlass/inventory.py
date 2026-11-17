"""Freeze the VLASS file inventory and independent public archive indexes."""
import csv
import hashlib
import io
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
CADC = 'https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync'
NRAO = 'https://data-query.nrao.edu/tap/sync'
STORAGE = 'https://ws-cadc.canfar.net/luskan/sync'
PRODUCT_QUERY = """SELECT uri,contentLength,contentChecksum FROM inventory.Artifact
WHERE uri LIKE 'nrao:VLASS/%V_mean_minmaxrej%.fits'
   OR uri LIKE 'nrao:VLASS/%.V.%.fits'
ORDER BY uri"""
QUERIES = {
    'v_products': (STORAGE, PRODUCT_QUERY),
    'cadc_polarization': (CADC, "SELECT pol_states,COUNT(*) AS n FROM ivoa.ObsCore WHERE obs_collection='VLASS' GROUP BY pol_states"),
    'nrao_polarization': (NRAO, "SELECT pol_states,COUNT(*) AS n FROM tap_schema.obscore WHERE dataproduct_type='image' AND project_code LIKE 'VLASS%' GROUP BY pol_states"),
    'cadc_caom_vmean': (CADC, "SELECT COUNT(*) AS n FROM caom2.Artifact WHERE uri LIKE 'nrao:VLASS/%V_mean_minmaxrej%.fits'"),
}
PAGES = {
    'cube_guide': 'https://science.nrao.edu/vlass/vlass-data/single-epoch-cube-users-guide',
    'cache_notice': 'https://vlass-dl.nrao.edu/vlass/se_cube_demo_release/',
    'cube_guide_cache_link': 'https://archive-new.nrao.edu/vlass/se_cube_demo_release/',
}


def run(key):
    if key in QUERIES:
        endpoint, query = QUERIES[key]
        url = endpoint + '?' + urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', MAXREC=1000000, QUERY=query))
    else:
        url = PAGES[key]
    try:
        response = urlopen(url, timeout=180)
    except HTTPError as exc:
        if key in QUERIES:
            raise
        response = exc
    with response:
        payload = response.read()
        record = dict(url=url, utc=datetime.now(timezone.utc).isoformat(), status=response.status,
                      n_bytes=len(payload), sha256=hashlib.sha256(payload).hexdigest())
    if key in QUERIES:
        if payload.lstrip().startswith(b'<'):
            raise RuntimeError(payload[:900])
        rows = list(csv.DictReader(io.StringIO(payload.decode())))
        record.update(query=QUERIES[key][1], rows=len(rows))
        if key == 'v_products':
            (ROOT/'products.csv').write_bytes(payload)
        else:
            record['result'] = rows
    else:
        record['contains_offline_notice'] = b'VLASS image cache offline' in payload
    return key, record


if __name__ == '__main__':
    with ThreadPoolExecutor(max_workers=7) as pool:
        audit = dict(pool.map(run, list(QUERIES) + list(PAGES)))
    (ROOT/'inventory_receipt.json').write_text(json.dumps(audit, indent=2) + '\n')
