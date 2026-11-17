"""Query DESI DR1, SDSS-V DR19/20 and LAMOST DR11/12 for the supplied targets.

Positional cones at Gaia J2016 coordinates include 15 arcsec of motion allowance.
Returned rows are possible matches. Consult the data release IDs and motion
before identifying any result with a target.
"""
import csv
from concurrent.futures import ThreadPoolExecutor
import hashlib
import io
import json
import math
from pathlib import Path
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ROOT=Path(__file__).parent
PARENT=ROOT.parent
CACHE=ROOT/'.cache'/'coverage'
CACHE.mkdir(exist_ok=True,parents=True)
GAIA={r['source_id']:r for name in ('control','search') for r in csv.DictReader((PARENT/name/'gaia.csv').open())}
RADIUS=15 / 3600.


def http(url):
    path=CACHE/hashlib.sha256(url.encode()).hexdigest()
    if not path.exists():
        with urllib.request.urlopen(url,timeout=65) as resp:
            raw=resp.read()
        path.write_bytes(raw)
    return path.read_bytes()


def cone(item):
    row, survey, base = item
    ra,dec=row['ra'],row['dec']
    if survey=='Legacy DR17':
        # Query SpecObjAll directly to include repeat plate observations.
        half_ra=RADIUS/math.cos(math.radians(abs(float(dec))+RADIUS))
        sql=("SELECT specObjID,ra,dec,plate,mjd,fiberID,class,subclass,z,zErr FROM SpecObjAll "
             f"WHERE ra BETWEEN {float(ra)-half_ra} AND {float(ra)+half_ra} AND "
             f"dec BETWEEN {float(dec)-RADIUS} AND {float(dec)+RADIUS} AND "
             f"dbo.fDistanceArcMinEq(ra,dec,{ra},{dec}) <= {60*RADIUS}")
        args=dict(cmd=sql,format='csv')
    elif survey.startswith('LAMOST'):
        args=dict(RA=ra,DEC=dec,SR=RADIUS,VERB=1)
    elif survey.startswith('SDSS'):
        args=dict(ra=ra,dec=dec,radius=RADIUS,units='degree',release=survey.split()[-1])
    else:
        # The NOIRLab TAP SQL translator does not map ICRS POINT on this
        # deployment. Use a bounding rectangle and check the circle locally.
        half_ra=RADIUS/math.cos(math.radians(abs(float(dec))+RADIUS))
        query=("SELECT targetid,survey,program,healpix,coadd_numexp,spectype,z,mean_fiber_ra,"
               f"mean_fiber_dec,mean_mjd FROM desi_dr1.zpix WHERE mean_fiber_ra BETWEEN {float(ra)-half_ra} "
               f"AND {float(ra)+half_ra} AND mean_fiber_dec BETWEEN {float(dec)-RADIUS} AND {float(dec)+RADIUS}")
        args=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=query)
    url=base+'?'+urllib.parse.urlencode(args)
    entry=dict(source_id=row['source_id'],survey=survey,url=url)
    try:
        data=http(url)
        entry.update(sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
        if survey=='Legacy DR17':
            table=data.decode().removeprefix('#Table1\n').removeprefix('#Table1\r\n')
            if not table.startswith('specObjID,'):
                raise ValueError(table[:300])
            records=list(csv.DictReader(io.StringIO(table)))
        elif survey.startswith('LAMOST'):
            doc=ET.fromstring(data)
            states=[e.attrib.get('value') for e in doc.findall('.//{*}INFO') if e.attrib.get('name')=='QUERY_STATUS']
            if states!=['OK'] or doc.find('.//{*}TABLEDATA') is None:
                raise ValueError(f'Invalid VOTable status or encoding: {states}')
            keys=[e.get('name') for e in doc.findall('.//{*}TABLE/{*}FIELD')]
            allrows=[dict(zip(keys,[e.text for e in tr])) for tr in doc.findall('.//{*}TR')]
            keep=['obsid','gaia_source_id','obsdate','designation','ra','dec','ra_obs','dec_obs','lmjd','planid','spid','fiberid','class','subclass','snrg']
            records=[{k:r.get(k) for k in keep if k in r} for r in allrows]
        elif survey.startswith('SDSS'):
            records=json.loads(data)
            if not isinstance(records,list):
                raise ValueError('Expected Valis target list')
        else:
            if not data.startswith(b'targetid,'):
                raise ValueError(data.decode()[:400])
            records=[]
            for match in csv.DictReader(io.StringIO(data.decode())):
                dra=math.radians(float(match['mean_fiber_ra'])-float(ra))
                ddec=math.radians(float(match['mean_fiber_dec'])-float(dec))
                a=(math.sin(ddec/2)**2 + math.cos(math.radians(float(dec)))*
                   math.cos(math.radians(float(match['mean_fiber_dec'])))*math.sin(dra/2)**2)
                sep=math.degrees(2*math.asin(math.sqrt(a)))*3600
                if sep <= RADIUS * 3600:
                    records.append({**match,'separation_arcsec':sep})
        entry.update(status='ok',rows=records)
    except Exception as err:
        entry.update(status='unavailable',reason=str(err)[:450])
    return entry


def main():
    surveys=[
        ('Legacy DR17','https://skyserver.sdss.org/dr17/SkyServerWS/SearchTools/SqlSearch'),
        ('SDSS DR20','https://api.sdss.org/valis/query/cone'),
        ('SDSS DR19','https://api.sdss.org/valis/query/cone'),
        ('LAMOST DR12 LRS','https://www.lamost.org/dr12/v2.0/voservice/conesearch'),
        ('LAMOST DR12 MRS','https://www.lamost.org/dr12/v2.0/medvoservice/conesearch'),
        ('LAMOST DR11 LRS','https://www.lamost.org/dr11/v2.0/voservice/conesearch'),
        ('LAMOST DR11 MRS','https://www.lamost.org/dr11/v2.0/medvoservice/conesearch'),
        ('DESI DR1','https://datalab.noirlab.edu/tap/sync')]
    work=[(row, name, url) for _,row in sorted(GAIA.items()) for name,url in surveys]
    with ThreadPoolExecutor(10) as pool:
        entries=list(pool.map(cone,work))
    with (ROOT/'coverage.jsonl').open('w') as f:
        for entry in entries:
            f.write(json.dumps(entry,separators=(',',':'))+'\n')
    for survey,_ in surveys:
        s=[r for r in entries if r['survey']==survey]
        print(survey,sum(r['status']=='ok' for r in s),[(r['source_id'],len(r['rows'])) for r in s if r.get('rows')],flush=True)


if __name__=='__main__':
    main()
