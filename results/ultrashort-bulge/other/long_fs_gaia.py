"""Resolve the pre-2024 Gaia variability and recover quality-filtered archive epochs."""
import csv
import io
import json
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode

from lxml import etree
import numpy as np

from fetch import SHA, OUT, get, write_csv

SOURCE = '4063555424531757184'
TAP = 'https://gea.esac.esa.int/tap-server/tap/sync?'
AIP = 'https://gaia.aip.de/tap/sync?'
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
        if name=='aip':
            doc=etree.fromstring(payload)
            assert doc.xpath('//*[local-name()="INFO" and @name="QUERY_STATUS" and @value="OK"]')
            fields=[f.get('name') for f in doc.xpath('//*[local-name()="FIELD"]')]
            rows=[dict(zip(fields,[cell.text for cell in row])) for row in doc.xpath('//*[local-name()="TR"]')]
        else:
            if not payload.startswith(b'solution_id,'):raise ValueError('unexpected TAP result')
            rows=list(csv.DictReader(io.StringIO(payload.decode())))
        result.update(status='ok',n_rows=len(rows))
        return result,rows
    except Exception as err:
        result.update(status='error',error=str(err))
        return result,[]


def summarize_epochs(row):
    out={}
    for band in ('g','bp','rp'):
        prefix='g_transit' if band=='g' else band
        t=np.array(row['g_transit_time' if band=='g' else band+'_obs_time'].split(),float)
        m=np.array(row[prefix+'_mag'].split(),float)
        accept=np.array([f.lower()=='false' for f in row['variability_flag_'+band+'_reject'].split()])
        use=np.isfinite(t)&np.isfinite(m)&accept
        if band!='g':
            use &= np.array([f.lower()=='false' for f in row['photometry_flag_'+band+'_reject'].split()])
        data=m[use];dates=t[use]
        out[band]=dict(epochs=len(data),min_time_gaia_bjd_minus_2455197p5=float(min(dates)),
                       max_time_gaia_bjd_minus_2455197p5=float(max(dates)),
                       range_mag=round(float(data.max()-data.min()),6),
                       p05_mag=round(float(np.quantile(data,.05)),6),p95_mag=round(float(np.quantile(data,.95)),6))
    return out


if __name__ == '__main__':
    tasks=[]
    for tab in ['vari_classifier_result','vari_summary','vari_long_period_variable','vari_eclipsing_binary','vari_microlensing']:
        q=f'SELECT * FROM gaiadr3.{tab} WHERE source_id={SOURCE}'
        tasks.append((tab,TAP+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))))
    q='select * from gaiadr3.epoch_photometry where source_id='+SOURCE
    tasks.append(('aip',AIP+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=q))))
    out={};receipts=[]
    with ThreadPoolExecutor(6) as workers:
        for receipt,rows in workers.map(retrieve,tasks):
            receipts.append(receipt)
            name=receipt['query']
            if not rows:
                continue
            if name=='aip':
                out['gaia_aip_epoch_photometry']=summarize_epochs(rows[0])
            elif name=='vari_summary':
                keys=('num_selected_g_fov','time_duration_g_fov','mean_mag_g_fov','std_dev_mag_g_fov','range_mag_g_fov','trimmed_range_mag_g_fov')
                out[name]=[{k:v for k,v in r.items() if k in keys or (k.startswith('in_') and v.lower()=='true')} for r in rows]
            else:
                out[name]=rows
    write_csv(OUT/'long_fs_gaia_queries.csv',receipts)
    (OUT/'long_fs_gaia.json').write_text(json.dumps(out,indent=2)+'\n')
