"""Resolve the pre-2024 Gaia variability designation at KMT-2024-BLG-1516."""
import csv
import gzip
import io
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode

from fetch import SHA, OUT, get, write_csv

SOURCE = '4063555424531757184'
TAP = 'https://gea.esac.esa.int/tap-server/tap/sync?'
CACHE = Path('.survey-cache/catalogs')


def retrieve(job):
    name,url=job
    path=CACHE/(SOURCE+'_'+name)
    result=dict(event='KMT-2024-BLG-1516',gaia_id=SOURCE,query=name,url=url)
    try:
        if not path.exists():
            path.write_bytes(get(url))
        payload=path.read_bytes()
        result.update(bytes=len(payload),sha256=SHA(payload))
        if payload[:2]==b'\x1f\x8b':payload=gzip.decompress(payload)
        if b'<VOTABLE' in payload or b'<html' in payload.lower() or b'QUERY_STATUS' in payload:raise ValueError('non-CSV response')
        rows=list(csv.DictReader(io.StringIO(payload.decode())))
        result.update(status='ok',n_rows=len(rows))
        return result,rows
    except Exception as err:
        result.update(status='error',error=str(err))
        return result,[]


if __name__ == '__main__':
    tasks=[]
    for tab in ['vari_classifier_result','vari_summary','vari_long_period_variable','vari_eclipsing_binary','vari_microlensing']:
        q=f'SELECT * FROM gaiadr3.{tab} WHERE source_id={SOURCE}'
        tasks.append((tab,TAP+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))))
    epoch_url='https://gea.esac.esa.int/data-server/data?'+urlencode(dict(ID='Gaia DR3 '+SOURCE,RETRIEVAL_TYPE='EPOCH_PHOTOMETRY',FORMAT='CSV',DATA_STRUCTURE='INDIVIDUAL',RELEASE='Gaia DR3'))
    tasks.append(('epoch',epoch_url))
    out={};receipts=[]
    with ThreadPoolExecutor(6) as workers:
        for receipt,rows in workers.map(retrieve,tasks):
            receipts.append(receipt)
            name=receipt['query']
            if name=='vari_classifier_result' and rows:
                out[name]=rows
            elif name=='vari_summary' and rows:
                out[name]=[{k:v for k,v in r.items() if (k in ('num_selected_g_fov','time_duration_g_fov','mean_mag_g_fov','std_dev_mag_g_fov','range_mag_g_fov','trimmed_range_mag_g_fov') or (k.startswith('in_') and v in ('true','True','t','1')))} for r in rows]
            elif rows and name!='epoch':
                out[name]=rows
            elif name=='epoch' and rows:
                import numpy as np
                out[name]={}
                for band in ('G','BP','RP'):
                    series=[r for r in rows if r.get('band')==band and r.get('rejected_by_photometry','false').lower()!='true' and r.get('rejected_by_variability','false').lower()!='true' and r.get('mag') not in (None,'','NaN')]
                    if not series:continue
                    a=np.array([[float(r['time']),float(r['mag'])] for r in series])
                    out[name][band]=dict(epochs=len(a),gaia_time_min=float(min(a[:,0])),gaia_time_max=float(max(a[:,0])),p05_mag=float(np.quantile(a[:,1],.05)),p95_mag=float(np.quantile(a[:,1],.95)))
    write_csv(OUT/'long_fs_gaia_queries.csv',receipts)
    (OUT/'long_fs_gaia.json').write_text(json.dumps(out,indent=2)+'\n')
