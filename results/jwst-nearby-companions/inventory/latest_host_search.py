"""Query SIMBAD identities and the NASA Exoplanet Archive for four release-priority hosts."""
import csv
import hashlib
import json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

OUT = Path(__file__).parent / 'archive'
HOSTS = ['G-99-49','L-788-34','Ross-619','V-V374-Peg']

def run(item):
    name, url, query = item
    params = {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':query}
    with urlopen(url + '?' + urlencode(params),timeout=60) as res:
        content=res.read()
    if b'<html' in content.lower() or b'VOTABLE' in content or not content:
        raise ValueError(f'{name}: not a CSV response: {content[:180]!r}')
    (OUT/(name+'.csv')).write_bytes(content)
    return {'service':url,'params':params,'sha256':hashlib.sha256(content).hexdigest(),
            'file':name+'.csv','retrieved_utc':datetime.now(timezone.utc).isoformat()}

if __name__=='__main__':
    with (OUT/'nearby-fields.csv').open() as f:
        targets={s['targprop']:s for s in csv.DictReader(f) if s['targprop'] in HOSTS and s['program']=='8826'}
    if set(targets)!=set(HOSTS): raise ValueError('incomplete host coordinates')
    cone=' OR '.join(f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{s['targ_ra']},{s['targ_dec']},0.016666667))" for s in targets.values())
    names=','.join("'Gaia DR3 " + s['gcns_edr3_source_id'] + "'" for s in targets.values())
    exo='SELECT DISTINCT hostname,pl_name,gaia_dr3_id,ra,dec FROM ps WHERE gaia_dr3_id IN ('+names+') OR '+cone
    simbad='SELECT main_id,otype,ra,dec,pmra,pmdec,plx_value,ids FROM basic JOIN ids ON oidref=oid WHERE '+cone
    ids='SELECT main_id,otype,ra,dec,pmra,pmdec,plx_value,ids FROM basic JOIN ids ON oidref=oid WHERE oid IN (SELECT oidref FROM ident WHERE id IN (' + names + '))'
    biblio='SELECT DISTINCT main_id,bibcode,title,doi FROM basic JOIN has_ref ON oid=has_ref.oidref JOIN ref ON oidbib=oidbibref WHERE oid IN (SELECT oidref FROM ident WHERE id IN ('+names+'))'
    tasks=[('go8826-host-bibliography','https://simbad.cfa.harvard.edu/simbad/sim-tap/sync',biblio),('go8826-simbad-hosts','https://simbad.cfa.harvard.edu/simbad/sim-tap/sync',ids),('go8826-confirmed-planets','https://exoplanetarchive.ipac.caltech.edu/TAP/sync',exo),('go8826-simbad','https://simbad.cfa.harvard.edu/simbad/sim-tap/sync',simbad)]
    with ThreadPoolExecutor(max_workers=2) as pool:
        result=list(pool.map(run,tasks))
    (OUT/'go8826-host-check-receipt.json').write_text(json.dumps(result,indent=2)+'\n')
