"""Re-acquire the Gaia parent, official crossmatches, WISE detections and bandpasses."""
import csv
import hashlib
import io
import json
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
GAIA='https://gea.esac.esa.int/tap-server/tap/sync'
IRSA='https://irsa.ipac.caltech.edu/TAP/sync'
RSR='https://irsa.ipac.caltech.edu/data/WISE/docs/release/All-Sky/expsup/figures/'
WISE_COLUMNS=('designation,ra,dec,w1mpro,w1sigmpro,w1snr,w1rchi2,w1sat,'
              'w2mpro,w2sigmpro,w3mpro,w3sigmpro,w3flux,w3sigflux,w3rchi2,w3sat,'
              'w4mpro,w4sigmpro,w4flux,w4sigflux,w4rchi2,w4sat,cc_flags,ph_qual,ext_flg,nb,na')


def rows(path):
    with (HERE/path).open() as stream:
        return list(csv.DictReader(stream))


def acquire(spec):
    name,url,q,col=spec
    for attempt in range(3):
        try:
            with urlopen(url,timeout=180) as response:
                data=response.read()
                code=response.status
            if col:
                reader=csv.DictReader(io.StringIO(data.decode()))
                if col not in (reader.fieldnames or []) or len(list(reader))>=10000:
                    raise ValueError('Response schema/truncation error')
            elif not data.startswith(b'# Relative Spectral Response per photon'):
                raise ValueError('Invalid bandpass response')
            (HERE/name).write_bytes(data)
            rec=dict(url=url,utc=datetime.now(timezone.utc).isoformat(),http=code,
                     bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
            if q: rec['adql']=q
            return name,rec
        except Exception:
            if attempt==2: raise
            time.sleep(2**attempt)


def tap(name,host,q,col):
    url=host+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC='10000',QUERY=q))
    return name,url,q,col


def main():
    receipts={}
    def run(tasks):
        with ThreadPoolExecutor(4) as pool:
            receipts.update(dict(pool.map(acquire,tasks)))
        (HERE/'catalog_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
    q=(HERE/'stars.adql').read_text()
    count='SELECT COUNT(*) AS n FROM gaiadr3.gaia_source AS g WHERE '+q.split('WHERE ',1)[1].split('ORDER BY')[0]
    run([tap('stars.csv',GAIA,q,'source_id'),tap('field_count.csv',GAIA,count,'n')])
    stars=rows('stars.csv')
    assert len({s['source_id'] for s in stars})==len(stars)==int(rows('field_count.csv')[0]['n'])
    tasks=[]
    for i in range(0,len(stars),100):
        ids=','.join(s['source_id'] for s in stars[i:i+100])
        q=('SELECT source_id,original_ext_source_id,angular_distance,number_of_neighbours,number_of_mates '
           f'FROM gaiadr3.allwise_best_neighbour WHERE source_id IN ({ids}) ORDER BY source_id')
        tasks.append(tap(f'xm_{i//100:03d}.csv',GAIA,q,'original_ext_source_id'))
    run(tasks)
    names=sorted({r['original_ext_source_id'] for task in tasks for r in rows(task[0])})
    tasks=[]
    for i in range(0,len(names),70):
        ids=','.join("'"+s+"'" for s in names[i:i+70])
        q=(f"SELECT {WISE_COLUMNS} FROM allwise_p3as_psd WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
           f"CIRCLE('ICRS',180,30,2.05)) AND designation IN ({ids}) ORDER BY designation")
        tasks.append(tap(f'wise_{i//70:03d}.csv',IRSA,q,'designation'))
    run(tasks)
    actual=[r['designation'] for task in tasks for r in rows(task[0])]
    if sorted(actual)!=names: raise ValueError('Incomplete IRSA designation set')
    run([(f'RSR-W{n}.txt',RSR+f'RSR-W{n}.txt',None,None) for n in [1,3,4]])


if __name__=='__main__':
    main()
