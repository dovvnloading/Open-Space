"""Inspect every public DR1/iron survey/program for redrock coadds in the query HEALPix cell."""
import argparse,csv,hashlib,io,json,re
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib import request
from urllib.error import HTTPError
from astropy_healpix import HEALPix
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits

BASE='https://data.desi.lbl.gov/public/dr1/spectro/redux/iron/healpix/'
CENTER=SkyCoord(53.18577725,-27.78513879,unit='deg')
RADIUS=10.8*u.arcsec

def links(url):
    try:
        with request.urlopen(url,timeout=40) as r: d=r.read()
    except HTTPError as error:
        if error.code==404:return [],dict(url=url,http_status=404)
        raise
    names=[s for s in re.findall(r'href="([^"]+)"',d.decode()) if s!='../']
    return names,dict(url=url,http_status=200,sha256=hashlib.sha256(d).hexdigest(),entries=len(names))

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('receipt',type=Path);a=p.parse_args()
    cells=HEALPix(nside=64,order='nested').cone_search_lonlat(CENTER.ra,CENTER.dec,RADIUS)
    surveys,root=links(BASE);surveys=[s for s in surveys if s.endswith('/')];audit=[root];rows=[]
    for survey in surveys:
        programs,rec=links(BASE+survey);audit.append(rec)
        def run(program):
            records=[];matches=[]
            for cell in cells:
                url=BASE+survey+program+f'{cell//100}/{cell}/'
                try: files,desc=links(url)
                except Exception as e: records.append(dict(url=url,error=str(e)));continue
                records.append(desc)
                for filename in files:
                    if filename.startswith('redrock-') and filename.endswith('.fits'):
                        full=url+filename
                        with request.urlopen(full,timeout=60) as r:d=r.read()
                        with fits.open(io.BytesIO(d)) as h:
                            z=h['REDSHIFTS'].data;fiber=h['FIBERMAP'].data
                            sky=SkyCoord(fiber['TARGET_RA'],fiber['TARGET_DEC'],unit='deg')
                            index=CENTER.separation(sky)<=RADIUS
                            ids=set(fiber['TARGETID'][index].tolist())
                            matches.extend([dict(targetid=int(row['TARGETID']),z=float(row['Z']),zwarn=int(row['ZWARN']),
                                                spectype=str(row['SPECTYPE']),file=filename) for row in z if row['TARGETID'] in ids])
                            records.append(dict(url=full,sha256=hashlib.sha256(d).hexdigest(),redshift_rows=len(z),sources_in_cone=len(ids)))
            return records,matches
        with ThreadPoolExecutor(4) as ex:
            for details,match in ex.map(run,[v for v in programs if v.endswith('/')]):audit.extend(details);rows.extend(match)
    a.receipt.write_text(json.dumps(dict(survey='DESI DR1 iron redrock',radius_arcsec=RADIUS.value,
                    nside=64,nested_cells=cells.tolist(),requests=audit,matches=rows),indent=2)+'\n')
if __name__=='__main__':main()
