"""Retrieve independent LAMOST stellar-parameter estimates for the long-period SB1."""
import hashlib
import json
import sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent.parent/'lamost'))
from fetch import query
COLUMNS = ('mobsid obsid lmjm band snr gaia_source_id teff_lasp teff_lasp_err logg_lasp logg_lasp_err '
           'feh_lasp feh_lasp_err teff_cnn logg_cnn feh_cnn coadd fibermask').split()

if __name__=='__main__':
    url,raw,rows=query(['275582835635329024'],'med_stellar',COLUMNS)
    (HERE/'lasp_visits.csv').write_bytes(raw)
    (HERE/'lasp_input.json').write_text(json.dumps(dict(url=url,rows=len(rows),sha256=hashlib.sha256(raw).hexdigest()),indent=2)+'\n')
