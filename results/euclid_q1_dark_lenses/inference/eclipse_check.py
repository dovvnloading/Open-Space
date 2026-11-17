"""Positionally test the public ECLIPSE-D v1 15,072-row DR1 products against both Q1 leads."""
import csv
import hashlib
import io
import urllib.request
import zipfile
import numpy as np
from infer import HERE, sep, write
from legacy_check import TARGETS

URL='https://zenodo.org/api/records/19654904/files/ECLIPSE_D_dr1_redshifts_and_physical_properties.zip/content'
def run():
    raw=urllib.request.urlopen(URL,timeout=120).read()
    if hashlib.md5(raw).hexdigest()!='cb31dbd65a3586ab13f52cff2df67f6b':
        raise ValueError('unexpected ECLIPSE-D v1 file')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        name='eclipse_dr1_catalogue_key_v20261103.csv'
        data=z.read(name)
    rows=list(csv.DictReader(io.StringIO(data.decode())))
    ra2=np.array([float(r['ra']) for r in rows])
    dec2=np.array([float(r['dec']) for r in rows])
    out=[]
    for key,(ra,dec) in TARGETS.items():
        dist=sep(ra,dec,ra2,dec2);idx=int(np.argmin(dist))
        out.append(dict(id_str=key,n_examined=len(rows),n_within_3p5arcsec=int(np.sum(dist<=3.5)),
                        nearest_name=rows[idx]['lens_name'],nearest_ra=ra2[idx],nearest_dec=dec2[idx],
                        separation_deg=dist[idx]/3600,archive_sha256=hashlib.sha256(raw).hexdigest(),
                        csv_sha256=hashlib.sha256(data).hexdigest()))
    write(HERE/'eclipse_nearest.csv',out)

if __name__=='__main__':
    run()
