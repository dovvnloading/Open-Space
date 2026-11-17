"""Freeze MAST proposal/filter query requests, response hashes, and the matching observation identifiers."""
import json,hashlib
from pathlib import Path
from urllib import request,parse
import argparse

p=argparse.ArgumentParser(description=__doc__);p.add_argument('outputs',type=Path)
p.add_argument('--responses',type=Path,help='directory of previously downloaded unmodified PROGRAM.json responses')
a=p.parse_args(); out=[]
for program in ['16912','17301','3149']:
    mission='JWST' if program=='3149' else 'HST'
    q=dict(service='Mast.Caom.Filtered',params={'columns':'obs_id,obsid,target_name,s_ra,s_dec,filters,t_exptime,obs_collection,proposal_id',
        'filters':[dict(paramName='obs_collection',values=[mission]),dict(paramName='proposal_id',values=[program])]},format='json',pagesize=200)
    url='https://mast.stsci.edu/api/v0/invoke?'+parse.urlencode({'request':json.dumps(q)})
    raw=(a.responses/f'{program}.json').read_bytes() if a.responses else request.urlopen(url,timeout=60).read()
    response=json.loads(raw)
    rows=[r for r in response['data'] if r['target_name'] in ('RCP28','WAKE','WAKEPOS1','WAKEPOS2')]
    out.append(dict(url=url,query=q,response_status=response['status'],rows_returned=len(response['data']),
                    paging=response.get('paging'),response_sha256=hashlib.sha256(raw).hexdigest(),relevant_observations=rows))
a.outputs.write_text(json.dumps(out,indent=2)+'\n')
