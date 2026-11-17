"""Query MER/PHZ and DESI for all VIS triggers and two high M/L positions; retain raw responses."""
import csv
import hashlib
import io
import json
import sys
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from infer import sep, write

HERE = Path(__file__).resolve().parent
DESI = 'https://datalab.noirlab.edu/tap/sync'
MER = 'https://eas.esac.esa.int/tap-server/tap/sync'
RADIUS = 3.5/3600

def run(profiles, out):
    centers = []
    for file in sorted(profiles.glob('*.json')):
        centers.extend((r['id'], r['ra'], r['dec']) for r in json.loads(file.read_text())['profiles'] if r['pass_vis'])
    centers.extend([('102158887_2655306065654633189',265.5306065631984,65.46331893363914),
                    ('102158888_2659753895653871052',265.9753895280555,65.38710523882253)])
    results=[]
    for ident, ra, dec in centers:
        w=RADIUS/np.cos(np.deg2rad(dec))
        condition='1=CONTAINS(POINT(\'ICRS\',m.right_ascension,m.declination),CIRCLE(\'ICRS\',%.9f,%.9f,%.10f))' % (ra,dec,RADIUS)
        query='SELECT m.object_id,m.right_ascension,m.declination,p.phz_median,p.phz_95_int,p.phz_flags,p.phz_classification FROM q1.mer_catalogue AS m LEFT OUTER JOIN q1.phz_photo_z AS p ON m.object_id=p.object_id WHERE '+condition+' ORDER BY m.object_id'
        results.append((ident,ra,dec,'euclid',MER,query))
        query=('SELECT targetid,mean_fiber_ra,mean_fiber_dec,z,zerr,zwarn,spectype,deltachi2,zcat_primary FROM desi_dr1.zpix '
               'WHERE mean_fiber_ra BETWEEN %.9f AND %.9f AND mean_fiber_dec BETWEEN %.9f AND %.9f ORDER BY targetid')%(ra-w,ra+w,dec-RADIUS,dec+RADIUS)
        results.append((ident,ra,dec,'desi',DESI,query))
    def request(task):
        ident,ra,dec,survey,endpoint,query=task
        url=endpoint+'?'+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query))
        data=urllib.request.urlopen(url,timeout=150).read()
        if data.startswith(b'<'):
            raise ValueError(data[:300])
        rows=list(csv.DictReader(io.StringIO(data.decode())))
        match=[]
        for r in rows:
            mer=survey=='euclid'
            ra2,dec2= float(r['right_ascension' if mer else 'mean_fiber_ra']),float(r['declination' if mer else 'mean_fiber_dec'])
            dist=sep(ra,dec,ra2,dec2)
            if dist <= 3.50001:
                match.append(dict(id=ident,survey=survey,ra_center=ra,dec_center=dec,offset_arcsec=dist,
                                  object_id=r['object_id' if mer else 'targetid'],ra=ra2,dec=dec2,
                                  redshift=r['phz_median' if mer else 'z'],PHZ_95=r.get('phz_95_int',''),
                                  flags=r.get('phz_flags',r.get('zwarn','')),classification=r.get('phz_classification',r.get('spectype',''))))
        return dict(id=ident,survey=survey,url=url,query=query,sha256=hashlib.sha256(data).hexdigest(),
                    raw_csv=data.decode(),matches=len(match)),match
    with ThreadPoolExecutor(4) as pool:
        responses=list(pool.map(request,results))
    (out/'cone_responses.json').write_text(json.dumps([m for m,_ in responses],indent=1)+'\n')
    write(out/'cone_matches.csv',[row for _,rows in responses for row in rows])

if __name__ == '__main__':
    run(Path(sys.argv[1]), HERE)
