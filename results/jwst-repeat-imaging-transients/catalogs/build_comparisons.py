"""Build NEXUS positional comparison lists from source snapshots; Python + lxml + astropy.

usage: python build_comparisons.py CACHE_DIR OUT_DIR
The dated catalog hashes in comparison_inputs.csv pin the exact source bytes used.
"""
import csv
import hashlib
import json
import sys
from pathlib import Path
from urllib.request import Request, urlopen
from math import cos, radians, sin, asin, sqrt
from shutil import copyfileobj

import numpy as np
from astropy.io import fits
from lxml import html

BASE = 'https://ariel.astro.illinois.edu/nexus/'
STONE = 'nexus_wide01_deep01_stacked_sources_allband.fits.gz'
TNS_CONE = 'https://transients.fyi/api/sky/markers?ra=268.53&dec=65.2&r=12&window=1Y&source=tns&measure=new&mag=all'
CACHE, OUT = map(Path, sys.argv[1:3])
CACHE.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
provenance = []
rows = []


def snapshot(name, url, post=None):
    local = CACHE / name
    if not local.exists():
        request = Request(url, data=json.dumps(post).encode() if post else None,
                          headers={'Content-Type': 'application/json'} if post else {})
        with urlopen(request, timeout=750) as response, local.with_suffix('.part').open('wb') as out:
            copyfileobj(response, out)
        local.with_suffix('.part').replace(local)
    with local.open('rb') as source:
        digest = hashlib.file_digest(source, 'sha256').hexdigest()
    provenance.append((name, url, local.stat().st_size, digest))
    return local


def table(filename, columns, body):
    with (OUT / filename).open('w', newline='') as handle:
        writer = csv.writer(handle, lineterminator='\n')
        writer.writerow(columns)
        writer.writerows(body)


article = html.parse(str(snapshot('nexus-year1.html', 'https://arxiv.org/html/2609.06985v1')))
for tr in article.xpath('//figure[@id="S4.T2"]//tbody/tr'):
    cells = tr.xpath('./td')
    for i in (0, 5):
        name = cells[i].text.strip()
        tentative = bool(cells[i].xpath('.//sup'))
        rows.append(('Zhuang2026_Table2', name, cells[i+1].text, cells[i+2].text,
                     'tentative' if tentative else 'robust'))
assert len(rows) == 74 and sum(r[-1] == 'robust' for r in rows) == 68

cone = json.loads(snapshot('transients-fyi-nexus.json', TNS_CONE).read_text())
if cone['truncated'] or cone['source'] != 'tns':
    raise ValueError('Incomplete or incorrect TNS broker cone response')
for obj in cone['objects']:
    rows.append(('TNS_TransientsFYI', 'AT '+obj['object_id'] if not obj.get('tns_type') else obj['object_id'],
                 obj['ra'], obj['dec'], 'one_year_rolling_cone'))

# Validate the full public broker response and include all returned objects within
# the same 12-arcmin cone. The two broker products have different source coverage.
request_args = {'resolver': 'tns', 'name_or_id': '', 'nmax': 1000000}
fink = json.loads(snapshot('fink-all.json', 'https://api.lsst.fink-portal.org/api/v1/resolver', request_args).read_text())
if not fink or len(fink) >= request_args['nmax']:
    raise ValueError('Truncated or empty Fink TNS response')
already = set()
for obj in fink:
    name = obj['f:fullname']
    ra, dec = float(obj['f:ra']), float(obj['f:declination'])
    if abs(dec - 65.2) > 0.2 or name in already:
        continue
    a = sin(radians(dec-65.2)/2)**2 + cos(radians(dec))*cos(radians(65.2))*sin(radians(ra-268.53)/2)**2
    if 2*asin(sqrt(a)) <= radians(.2):
        rows.append(('TNS_Fink', name, ra, dec, 'public_resolver_cone'))
        already.add(name)

edrname = 'wide_ep1_01_src_catalog_v1.fits'
snapshot(edrname, BASE+'edr/nircam/catalog/'+edrname)

stone = fits.getdata(snapshot(STONE, BASE+'paper_data/nuclear_variability/'+STONE))
selected = np.zeros(len(stone), dtype=bool)
for band in ('F200W', 'F444W'):
    selected |= stone['VARIABLE_'+band] & stone['MASK_VISUAL_INSPECTION_'+band]
for obj in stone[selected]:
    rows.append(('Stone2026', 'ID '+str(obj['ID']), obj['RA'], obj['DEC'], 'published_variable'))
assert selected.sum() == 465

for epoch in range(1,9):
    filename = f'nexus_msa_deep_ep{epoch:02}_targets.fits'
    data = fits.getdata(snapshot(filename, BASE+'qdr/nirspec/'+filename))
    for obj in data[(data['Priority Class'] >= 0) & (data['Priority Class'] < 1)]:
        rows.append(('NEXUS_MSA_priority0', 'MSAID '+str(obj['MSAID']), obj['RA'], obj['DEC'],
                     f"Deep{epoch:02}_class{obj['Priority Class']:.1f}"))

# Keep duplicates across sources and epochs for auditing; MSA targets are proximity
# flags, not secure transient identifications.
table('field_comparisons.csv', ['catalog','designation','ra_deg','dec_deg','status'], rows)
table('comparison_inputs.csv', ['local_name','url','bytes','sha256'], provenance)
