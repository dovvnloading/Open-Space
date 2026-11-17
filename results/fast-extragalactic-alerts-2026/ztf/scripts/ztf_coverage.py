"""Context counts from Fink, independent of the TAP object sample."""
import datetime as dt
import hashlib
import json
import urllib.request
from ztf_tap import META
def main():
    url='https://api.ztf.fink-portal.org/api/v1/statistics'
    params={'date':'202611','output-format':'json'}
    req=urllib.request.Request(url,data=json.dumps(params).encode(),headers={'Content-Type':'application/json'})
    with urllib.request.urlopen(req,timeout=45) as r:
        raw=r.read()
        receipt=dict(url=url,method='POST',params=params,retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(),http_status=r.status,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))
    (META/'fink_statistics_body.json').write_bytes(raw)
    (META/'fink_statistics_request.json').write_text(json.dumps(receipt,indent=2)+'\n')
    for key in ('basic:raw','basic:sci','basic:exposures'):
        print(key,sum(int(x[key]) for x in json.loads(raw) if 'ztf_20261101' <= x['key:key'] < 'ztf_20261117'))
if __name__=='__main__':
    main()
