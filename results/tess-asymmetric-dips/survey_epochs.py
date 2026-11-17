"""Query public ASAS-SN v2 and ZTF LC coverage at the Sector 105 trigger positions.

Usage: python survey_epochs.py triggers.csv independent.json
TESS BTJD ~4205-4233 corresponds to MJD ~61204.5-61232.5.
The wider MJD window below is used only for coverage, not precise timing.
"""
import csv
import hashlib
import json
import sys
from datetime import datetime, timezone
from io import StringIO
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ASAS = 'http://asassn-lb01.ifa.hawaii.edu:9006'

def fetch(url, body=None):
    rec = {'url': url, 'utc': datetime.now(timezone.utc).isoformat(), 'payload': body}
    request = Request(url, data=json.dumps(body).encode() if body else None,
                      headers={'Content-Type': 'application/json'} if body else {})
    try:
        with urlopen(request, timeout=150) as response:
            raw = response.read()
            rec.update(status=response.status, bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest())
            return rec, raw
    except (HTTPError, OSError) as error:
        rec['error'] = str(error)
        return rec, None

def main():
    with open(sys.argv[1]) as f:
        stars = {r['tic']: r for r in csv.DictReader(f)}
    query = dict(catalog='stellar_main', cols=['asas_sn_id', 'tic_id', 'ra_deg', 'dec_deg'],
                 format='json', id_col='tic_id', n_rows=len(stars)+1, page_num=0, tar_ids=list(stars))
    rec, raw = fetch(ASAS+'/lookup_targets/catalog_list', query)
    if raw is None:
        raise RuntimeError(rec)
    matches = {str(r['tic_id']):r for r in json.loads(raw)['data']}
    out = {'lookup': rec, 'targets': {}}
    for tic,row in stars.items():
        info = {}
        out['targets'][tic] = info
        match = matches.get(tic)
        if match is None:
            url = ASAS+'/lookup_cone/radius0.002222222222_ra'+row['ra_obj']+'_dec'+row['dec_obj']
            query = dict(catalog='master_list', cols=['asas_sn_id','ra_deg','dec_deg'],
                         format='json', page_num=0, n_rows=50)
            info['cone_request'], raw = fetch(url,query)
            found = json.loads(raw)['data'] if raw else []
            info['cone_matches'] = found
            if len(found)==1:
                match=found[0]
        if match is not None:
            info['match'] = match
            info['asas_request'], raw = fetch(ASAS+'/get_lightcurve/'+str(match['asas_sn_id']))
            if raw:
                doc=json.loads(raw)
                if 'light_curve' in doc:
                    lc=doc['light_curve']
                    jd=lc['columns'].index('jd')
                    times=[r[jd] for r in lc['data']]
                    info['asas_summary']=dict(rows=len(times), last_jd=max(times,default=None),
                        sector_window=sum(61204<=t-2400000.5<=61234 for t in times))
                else:
                    info['asas_service_error']=doc.get('error','no light_curve in response')
        query = dict(POS='CIRCLE '+row['ra_obj']+' '+row['dec_obj']+' 0.001667',
                     TIME='61204 61234', FORMAT='csv')
        info['ztf_request'],raw=fetch('https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?'+urlencode(query))
        if raw:
            table=csv.DictReader(StringIO(raw.decode()))
            if 'mjd' in (table.fieldnames or []):
                info['ztf_in_window']=sum(1 for _ in table)
            else:
                info['ztf_service_error']='missing mjd header'
        Path(sys.argv[2]).write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':
    main()
