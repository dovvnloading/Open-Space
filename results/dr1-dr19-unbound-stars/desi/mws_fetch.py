"""Query DESI DR1 MWS stellar-fit coadds and retain the original, lossless CSV responses."""
import argparse, csv, hashlib, io, json, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
from urllib.error import HTTPError
COLS = ('targetid source_id ref_id ref_cat target_ra target_dec glon glat healpix survey program '
        'zcat_primary vrad vrad_err vrad_skew vrad_kurt rvs_warn success teff teff_err logg logg_err '
        'feh feh_err rr_z rr_spectype sn_b sn_r sn_z chisq_tot chisq_c_tot chisq_b chisq_c_b '
        'chisq_r chisq_c_r chisq_z chisq_c_z rv_adop rv_err').split()
def extract(ids, out, name):
    adql = 'SELECT '+','.join(COLS)+' FROM desi_dr1.mws WHERE targetid IN ('+','.join(ids)+')'
    params = {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','MAXREC':50000,'QUERY':adql}
    base = 'https://datalab.noirlab.edu/tap/sync'
    path = out/(name+'.csv')
    for attempt in range(3):
        try:
            with urlopen(base+'?'+urlencode(params), timeout=150) as reply: raw=reply.read()
            assert raw.startswith(','.join(COLS).encode()), raw[:250]
            path.write_bytes(raw)
            break
        except Exception:
            if attempt==2: raise
            time.sleep(2**attempt)
    (out/(name+'.adql')).write_text(adql+'\n')
    return list(csv.DictReader(io.StringIO(raw.decode())))
def main():
    p=argparse.ArgumentParser();p.add_argument('screen',type=Path);p.add_argument('out',type=Path);args=p.parse_args()
    args.out.mkdir(exist_ok=True,parents=True)
    with args.screen.open() as f: stars=list(csv.DictReader(f))
    # Fetch nearby candidates first so that their exposure checks can begin independently.
    ids=list(dict.fromkeys(r['TARGETID'] for r in sorted(stars,key=lambda r:-int(r['nearby']))))
    chunks=[ids[i:i+64] for i in range(0,len(ids),64)]
    with ThreadPoolExecutor(6) as pool: batches=list(pool.map(lambda it: extract(it[1],args.out,'coadds_%02d'%it[0]),enumerate(chunks)))
    with (args.out/'coadds.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=COLS);w.writeheader()
        for batch in batches:w.writerows(batch)
    print('requested IDs',len(ids),'returned rows',sum(map(len,batches)),'unique',len({r['targetid'] for b in batches for r in b}),flush=True)
if __name__=='__main__':main()
