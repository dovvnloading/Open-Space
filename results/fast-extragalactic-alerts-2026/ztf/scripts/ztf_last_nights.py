"""Check the empty UW public archives after the final nonempty Fink night."""
import base64
import datetime as dt
import hashlib
import io
import json
import tarfile
from urllib.request import urlopen
from ztf_tap import ROOT

def main():
    receipts = []
    for day in range(13,17):
        url = f'https://ztf.uw.edu/alerts/public/ztf_public_202611{day}.tar.gz'
        with urlopen(url,timeout=30) as response:
            body = response.read(1025)
            if len(body)>1024:
                raise ValueError(f'{url} has an updated nontrivial archive; inspect before interpreting.')
            with tarfile.open(fileobj=io.BytesIO(body),mode='r:gz') as tar:
                members = tar.getnames()
            receipts.append(dict(url=url,method='GET',http_status=response.status,
                                retrieved_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                                bytes=len(body),sha256=hashlib.sha256(body).hexdigest(),
                                base64=base64.b64encode(body).decode(),members=members))
    path = ROOT / 'derived' / 'ztf_latest_archives.json'
    path.write_text(json.dumps(receipts,indent=2)+'\n')
if __name__=='__main__':
    main()
