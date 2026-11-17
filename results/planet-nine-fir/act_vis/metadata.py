"""Obtain DR1 exposure metadata and the Euclid observer's barycentric ephemeris."""
import csv
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent


def retrieve(job):
    name,url=job
    path=HERE/(name+'.response')
    if not path.exists():
        with urlopen(url,timeout=120) as resp:
            data=resp.read()
            path.write_bytes(data)
            receipt=dict(url=url,status=resp.status,bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),utc=datetime.now(timezone.utc).isoformat())
            (HERE/(name+'_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    return path


def main():
    sql=("SELECT obs_id,t_min,t_max,t_exptime,s_ra,s_dec,s_region,access_url FROM ivoa.obscore "
         "WHERE obs_collection='dr1' AND obs_id='6578' AND t_exptime>560 AND t_exptime<561")
    tap='https://eas.esac.esa.int/tap-server/tap/sync?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=sql,MAXREC=100))
    horizons='https://ssd.jpl.nasa.gov/api/horizons.api?'+urlencode(dict(
        format='json',COMMAND="'-680'",CENTER="'500@0'",MAKE_EPHEM="'YES'",EPHEM_TYPE="'VECTORS'",START_TIME="'2025-02-24'",STOP_TIME="'2025-02-26'",STEP_SIZE="'4 h'",OUT_UNITS="'AU-D'",REF_PLANE="'FRAME'",VEC_TABLE="'1'",CSV_FORMAT="'YES'",VEC_CORR="'NONE'"))
    listing='https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/?'+urlencode({'list-type':2,'prefix':'dr1/VIS_QUAD/6578/'})
    with ThreadPoolExecutor(max_workers=3) as pool:
        list(pool.map(retrieve,[('euclid_dr1',tap),('horizons',horizons),('s3_listing',listing)]))

if __name__=='__main__':main()
