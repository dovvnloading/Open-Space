"""Retrieve independent DR9 Legacy Surveys blended photometric-redshift checks."""
import csv
import hashlib
import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

from infer import sep, write

HERE=Path(__file__).resolve().parent
TARGETS={'102158887_2655306065654633189':(265.5306065631984,65.46331893363914),
         '102158888_2659753895653871052':(265.9753895280555,65.38710523882253)}


def run():
    records, matches=[],[]
    for key,(ra,dec) in TARGETS.items():
        q=('SELECT t.ls_id,t.ra,t.dec,t.type,t.flux_g,t.flux_r,t.flux_z,t.flux_w1,t.flux_w2,'
           'p.z_phot_median,p.z_phot_l68,p.z_phot_u68,p.z_phot_l95,p.z_phot_u95,p.z_spec '
           'FROM ls_dr9.tractor AS t LEFT OUTER JOIN ls_dr9.photo_z AS p ON t.ls_id=p.ls_id '
           'WHERE t.ra BETWEEN %.8f AND %.8f AND t.dec BETWEEN %.8f AND %.8f ORDER BY t.ls_id')%(ra-.003,ra+.003,dec-.0012,dec+.0012)
        url='https://datalab.noirlab.edu/tap/sync?'+urllib.parse.urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))
        b=urllib.request.urlopen(url,timeout=150).read()
        if not b.startswith(b'ls_id,'):
            raise ValueError(b[:400])
        records.append(dict(id=key,query=q,url=url,raw_csv=b.decode(),sha256=hashlib.sha256(b).hexdigest()))
        for row in csv.DictReader(io.StringIO(b.decode())):
            d=sep(ra,dec,float(row['ra']),float(row['dec']))
            if d<4:
                matches.append(dict(id_str=key,delta_arcsec=d,**row))
    (HERE/'legacy_responses.json').write_text(json.dumps(records,indent=1)+'\n')
    write(HERE/'legacy_matches.csv',matches)


if __name__=='__main__':
    run()
