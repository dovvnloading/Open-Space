"""Retrieve the published SFRS galaxy entry and independent SDSS line diagnostics."""
import hashlib
import json
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

OUT=Path(__file__).resolve().parent.parent/'wide'/'catalogs'
BASE='https://skyserver.sdss.org/dr16/SkyServerWS/SearchTools/SqlSearch'
Q="""SELECT TOP 10 s.specObjID,s.ra,s.dec,s.z,s.zWarning,s.class,s.subClass,
g.h_alpha_flux,g.h_alpha_flux_err,g.h_beta_flux,g.h_beta_flux_err,
g.nii_6584_flux,g.nii_6584_flux_err,g.oiii_5007_flux,g.oiii_5007_flux_err
FROM dbo.fGetNearbySpecObjEq(209.0330736,30.0813453,0.25) n
JOIN SpecObj s ON n.specObjID=s.specObjID
LEFT JOIN galSpecLine g ON s.specObjID=g.specObjID ORDER BY n.distance"""
def fetch(name,url):
    with urlopen(url,timeout=90) as r:data=r.read();status=r.status
    if name=='sfrs':
        content=data.decode('utf8','replace').splitlines()
        lines=[k for k in content if k.split()[:1]==['293']]
        if len(lines)!=1:raise ValueError('Cannot identify unique SFRS row')
        with urlopen('https://lweb.cfa.harvard.edu/irac/SFRS/table1.txt',timeout=90) as response:
            catalog=response.read()
        identity=[s for s in catalog.decode().splitlines() if s.split()[:1]==['293']]
        if len(identity)!=1 or '8856' not in identity[0]:
            raise ValueError('Could not cross-identify SFRS 293 with UGC 8856')
        receipt=dict(entry=lines[0],source='Mahajan et al. SFRS III, machine-readable Table 2',
                     columns=['SFRS','log_SFR_FUV','log_SFR_NUV','log_SFR_1.4GHz','log_SFR2_1.4GHz',
                              'log_SFR_FIR','log_SFR_PAH','log_SFR_tot'],
                     table1_entry=identity[0],table1_url='https://lweb.cfa.harvard.edu/irac/SFRS/table1.txt',
                     table1_sha256=hashlib.sha256(catalog).hexdigest())
    else:
        import csv,io
        content=data.decode('utf8','replace')
        if 'h_alpha_flux' not in content or 'error' in content[:220].lower():
            raise ValueError(content[:300])
        lines=content.splitlines()
        offset=next(i for i,x in enumerate(lines) if 'h_alpha_flux' in x)
        rows=list(csv.DictReader(io.StringIO('\n'.join(lines[offset:]))))
        if not rows:raise ValueError('No usable SDSS rows')
        receipt=dict(query=Q,rows=rows)
    receipt.update(url=url,status=status,utc=datetime.now(timezone.utc).isoformat(),bytes=len(data),
                   sha256=hashlib.sha256(data).hexdigest())
    (OUT/('ugc8856_'+name+'.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    print(name,receipt.get('entry') or receipt.get('rows'))
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('source',choices=['sfrs','sdss']);a=p.parse_args()
    fetch(a.source,'https://lweb.cfa.harvard.edu/irac/SFRS/table2.txt' if a.source=='sfrs'
          else BASE+'?'+urlencode({'cmd':Q,'format':'csv'}))
