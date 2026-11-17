"""JPL second-pass asteroid/comet check for the entire Euclid VIS search area."""
import json
from urllib.parse import urlencode
from urllib.request import urlopen
from datetime import datetime,timezone
import hashlib

import erfa
import numpy as np
from vis_pixels import HERE,TREF,SC

def main():
    t=TREF.tdb.jd-2451545.
    _,earth=erfa.epv00(2451545.,t)
    state=np.r_[(SC(t)-earth['p'])*149597870.7,(SC(t,1)-earth['v'])*149597870.7/86400]
    url='https://ssd-api.jpl.nasa.gov/sb_ident.api?'+urlencode({
        'xobs':','.join(f'{x:.9e}' for x in state),
        'obs-time':'2024-12-27T07:30:00','fov-ra-center':'02-19-46.8','fov-dec-center':'M49-41-24.0',
        'fov-ra-hwidth':0.15,'fov-dec-hwidth':0.12,'mag-required':'false','two-pass':'true',
        'suppress-first-pass':'true','req-elem':'false'})
    with urlopen(url,timeout=180) as f:
        data=f.read()
        receipt=dict(url=url,utc=datetime.now(timezone.utc).isoformat(),status=f.status,
                     bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    (HERE/'known_objects.response').write_bytes(data)
    (HERE/'known_objects_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(data[:1000])

if __name__=='__main__':main()
