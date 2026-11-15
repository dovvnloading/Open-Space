"""Repeat the identifier checks against three published Kepler catalogs."""
import re
from urllib.parse import urlencode
from urllib.request import urlopen

for table in ['J/AJ/151/101', 'J/AJ/151/68', 'J/ApJS/235/38/table8']:
    url='https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?'+urlencode({
        '-source':table,'-out.max':'unlimited','-out.all':1})
    raw=urlopen(url, timeout=120).read()
    lines=[line for line in raw.decode().splitlines() if line and not line.startswith('#')]
    assert len(lines)>100 and any('KIC' in line or 'TCE' in line for line in lines[:10]), table
    print(table, len(lines), [line for line in lines if re.search(r'811075[478]|786[56]\.01',line)])
