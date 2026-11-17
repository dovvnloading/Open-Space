"""Retrieve HST/ACS HRC historical frames covering HD 38206; keep the original archive response and pixel hash."""
import argparse,csv,json,hashlib
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
from datetime import datetime,timezone

def query(out):
    params=dict(service='Mast.Caom.Filtered',format='json',params=dict(columns='obsid,obs_id,target_name,proposal_id,filters,t_min,t_max,calib_level,instrument_name,dataRights,dataURL',filters=[{'paramName':'proposal_id','values':['10539']},{'paramName':'target_name','values':['HD-38206','HD38206','HD 38206','HD-41695','HD41695','HD 41695']},{'paramName':'obs_collection','values':['HST','HLA']},{'paramName':'instrument_name','values':['ACS/HRC']}]),pagesize=3000,page=1)
    url='https://mast.stsci.edu/api/v0/invoke?'+urlencode({'request':json.dumps(params)})
    with urlopen(url,timeout=90) as f:raw=f.read()
    reply=json.loads(raw);assert reply['status']=='COMPLETE'
    (out/'HD38206-hst.json').write_text(json.dumps(dict(utc=datetime.now(timezone.utc).isoformat(),query=url,sha256=hashlib.sha256(raw).hexdigest(),records=reply['data']),indent=2)+'\n')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('out',type=Path);query(p.parse_args().out)
