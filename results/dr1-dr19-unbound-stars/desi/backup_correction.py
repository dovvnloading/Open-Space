"""Apply the published DESI DR1 backup-program correction by DESI target ID."""
import csv, hashlib, io, sys, json
from pathlib import Path
from urllib.request import urlopen
from astropy.io import fits
import numpy as np
OUT=Path(__file__).resolve().parent
URL='https://zenodo.org/api/records/15469272/files/backup_correction.fits/content'
SHA='eb4da91267db39f285a277989489f991ea48371336efedd28bd07d2b58e4a400'
if __name__=='__main__':
    data=Path(sys.argv[1]).read_bytes() if len(sys.argv)>1 else urlopen(URL,timeout=120).read()
    assert hashlib.sha256(data).hexdigest()==SHA
    with fits.open(io.BytesIO(data)) as f:
        r=f[1].data
        source=next(x for x in csv.DictReader((OUT/'los/leads.csv').open()) if x['source_id']=='1717410430709663616')
        off=float(r['VRAD_BIAS'][r['TARGETID']==int(source['targetid'])].item())
    (OUT/'los/backup_correction.json').write_text(json.dumps({'source_id':source['source_id'],'targetid':source['targetid'],'vrad_raw':float(source['vrad']),'vrad_bias':off,'vrad_corrected':float(source['vrad'])-off,'url':URL,'sha256':SHA},indent=2)+'\n')
