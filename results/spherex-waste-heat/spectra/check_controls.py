"""AllWISE W1–W3 quality, colors, and PM-aware associations for comparison stars."""
import csv
import hashlib
import io
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'wise'))
from fetch import TABLES,TAP,position,separation

def retrieve(target):
    table,cols,epoch=TABLES['allwise']
    ra,dec=position(target,epoch)
    q=(f"SELECT {cols} FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),"
       f"CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.003333333))=1")
    url=TAP+'?'+urlencode(dict(QUERY=q,FORMAT='csv',LANG='ADQL',REQUEST='doQuery'))
    with urlopen(url,timeout=120) as resp: raw=resp.read()
    rows=list(csv.DictReader(io.StringIO(raw.decode())))
    matches=sorted([(separation(float(r['ra']),float(r['dec']),ra,dec),r) for r in rows],key=lambda a:a[0])
    result=dict(label=target['label'],source_id=target['source_id'],n_within12=len(rows),
                nearest_sep_arcsec='',designation='',w1_w2='',w2_w3='',w3_snr='',w4_snr='',status='no_match')
    if matches:
        sep,r=matches[0]
        result.update(nearest_sep_arcsec=round(sep,3),designation=r['designation'])
        if sep<=2.0 and len([x for x in matches if x[0]<=6])==1:
            w12=float(r['w1mpro'])-float(r['w2mpro'])
            result.update(w1_w2=round(w12,3),w3_snr=r['w3snr'],w4_snr=r['w4snr'])
            if r['w3mpro'] and r['w3snr'] and float(r['w3snr'])>=3:
                result['w2_w3']=round(float(r['w2mpro'])-float(r['w3mpro']),3)
            clean=(r['cc_flags'][:2]=='00' and r['ph_qual'][:2]=='AA' and int(r['ext_flg'])==0 and int(r['nb'])==1)
            result['status']='ambiguous'
            if clean and -.2<=w12<=.45:
                if result['w2_w3']=='': result['status']='W12_photospheric_W3_unconstrained'
                elif r['cc_flags'][2]=='0' and -.3<=result['w2_w3']<=.8: result['status']='W123_photospheric'
            elif clean: result['status']='color_outlier'
    receipt=dict(label=target['label'],url=url,sha256=hashlib.sha256(raw).hexdigest(),rows=len(rows),bytes=len(raw))
    return result,receipt

def main():
    fields=[]
    for name in ('positions.csv','positions_fgk.csv'):
        fields.extend(r for r in csv.DictReader((ROOT/name).open()) if '_cmp' in r['label'])
    existing={}
    if (ROOT/'comparison_wise.csv').exists() and (ROOT/'comparison_queries.json').exists():
        meta={r['label']:r for r in json.loads((ROOT/'comparison_queries.json').read_text())}
        existing={r['source_id']:(r,meta[r['label']]) for r in csv.DictReader((ROOT/'comparison_wise.csv').open()) if r['label'] in meta}
    new=[row for row in fields if row['source_id'] not in existing]
    with ThreadPoolExecutor(6) as pool: results=list(existing.values())+list(pool.map(retrieve,new))
    with (ROOT/'comparison_wise.csv').open('w') as out:
        writer=csv.DictWriter(out,fieldnames=results[0][0]); writer.writeheader()
        writer.writerows(r for r,_ in results)
    (ROOT/'comparison_queries.json').write_text(json.dumps([v for _,v in results],indent=2)+'\n')
if __name__=='__main__': main()
