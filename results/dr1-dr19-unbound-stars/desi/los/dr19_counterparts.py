"""Retrieve all Gaia-ID matches to priority DESI targets without DR19 RV quality prefilters."""
from pathlib import Path
from urllib.request import urlopen
from urllib.parse import urlencode
from concurrent.futures import ThreadPoolExecutor
import csv
out=Path(__file__).resolve().parent
ids=set(r['source_id'] for r in csv.DictReader((out/'leads.csv').open()))
ids.update(['2752545725694798336','6340918677332110720','2294414169958935296','4435527217525195264','4550108286326235520'])
def fetch(tab):
    sql=f"SELECT gaia_dr3_source_id,sdss_id,spectrum_pk,ra,dec,snr,v_rad,e_v_rad FROM {tab} WHERE gaia_dr3_source_id IN ({','.join(sorted(ids))})"
    url='https://skyserver.sdss.org/dr19/SkyServerWS/SearchTools/SqlSearch?'+urlencode({'cmd':sql,'format':'csv'})
    with urlopen(url,timeout=150) as f:b=f.read()
    assert b.startswith(b'#Table1\ngaia_dr3_source_id'),b[:250]
    (out/(tab+'.sql')).write_text(sql+'\n');(out/(tab+'.csv')).write_bytes(b)
    print(tab,len(b.splitlines())-2,flush=True)
if __name__=='__main__':
    with ThreadPoolExecutor(2) as pool:list(pool.map(fetch,['mwm_boss_allstar','mwm_apogee_allstar']))
