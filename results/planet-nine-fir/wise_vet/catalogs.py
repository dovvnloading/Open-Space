"""Query independent AllWISE and Gaia DR3 catalogues near both W2 marks."""
import hashlib,json
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
HERE=Path(__file__).resolve().parent

def query(job):
    label,u=job
    receipt=dict(url=u,utc=datetime.now(timezone.utc).isoformat())
    try:
        with urlopen(u,timeout=150 if 'jpl' in label else 60) as r:
            z=r.read();receipt.update(status=r.status,bytes=len(z),sha256=hashlib.sha256(z).hexdigest())
            (HERE/(label+'.response')).write_bytes(z)
    except Exception as exc:receipt['error']=str(exc)
    (HERE/(label+'_receipt.json')).write_text(json.dumps(receipt,indent=2)+'\n')
    print(label,receipt.get('status'),receipt.get('error'),flush=True)

def main():
    trials=json.loads((HERE/'positions.json').read_text()); jobs=[]
    for v in range(3):
        a=trials[0]['visits'][v];ra=a['ra'];de=a['dec']
        for label,tbl,select,b in [
            ('allwise','allwise_p3as_psd','designation,ra,dec,w1mpro,w2mpro,w2sigmpro,w2snr,cc_flags,var_flg','https://irsa.ipac.caltech.edu/TAP/sync'),
            ('gaia','gaiadr3.gaia_source','source_id,ra,dec,phot_g_mean_mag,pmra,pmdec,phot_variable_flag','https://gea.esac.esa.int/tap-server/tap/sync')]:
            sql=f"SELECT {select} FROM {tbl} WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra},{de},0.065))"
            jobs.append((label+str(v),b+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=20000,QUERY=sql))))
    sql = "SELECT * FROM gaiadr3.vari_classifier_result WHERE source_id=2193806661294321792"
    jobs.append(('gaia_variable', 'https://gea.esac.esa.int/tap-server/tap/sync?' + urlencode(
        dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY=sql))))
    with ThreadPoolExecutor(max_workers=6) as e: list(e.map(query,jobs))
if __name__=='__main__':main()
