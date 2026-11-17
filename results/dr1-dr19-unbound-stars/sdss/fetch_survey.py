"""Fetch public DR19 MWM tables. Run with Python 3, from this file's directory."""
import csv
import hashlib
import io
import sys
import time
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor

HERE = Path(__file__).resolve().parent
WORK = HERE / 'raw'
WORK.mkdir(exist_ok=True)
ENDPOINT = 'https://skyserver.sdss.org/dr19/SkyServerWS/SearchTools/SqlSearch'
TABLES = ['mwm_boss_allstar', 'mwm_apogee_allstar']
BASE = ("release='sdss5' AND gaia_dr3_source_id>0 AND snr>=15 "
        "AND e_v_rad>0 AND e_v_rad<=25 AND ABS(v_rad)<5000")
# Loose heliocentric speed component bounds to reduce the downloaded table size.
# The no-parallax radial branch retains distant, high-radial-velocity objects.
PREFILTER = """(ABS(v_rad)>=275 OR
 (plx>0 AND e_plx>0 AND plx>=4*e_plx AND
 ABS(pmra)<10000 AND ABS(pmde)<10000 AND
 22.47206*(SQUARE(pmra)+SQUARE(pmde))>=SQUARE(275.0*plx)) OR
 (r_med_photogeo>0 AND ABS(pmra)<10000 AND ABS(pmde)<10000 AND
 0.00002247206*(SQUARE(pmra)+SQUARE(pmde))*SQUARE(r_med_photogeo)>=75625))"""
COMMON = ('sdss_id,gaia_dr3_source_id,spectrum_pk,ra,dec,l,b,'
          'plx,e_plx,pmra,e_pmra,pmde,e_pmde,r_med_photogeo,r_lo_photogeo,r_hi_photogeo,'
          'r_med_geo,r_lo_geo,r_hi_geo,gaia_v_rad,gaia_e_v_rad,g_mag,bp_mag,rp_mag,'
          'telescope,snr,v_rad,e_v_rad,std_v_rad,min_mjd,max_mjd,n_good_visits,n_good_rvs')
EXTRA = {'mwm_boss_allstar': 'zwarning_flags,xcsao_meanrxc',
         'mwm_apogee_allstar':'spectrum_flags,doppler_flags,n_components'}

def query(sql):
    url = ENDPOINT + '?' + urlencode({'cmd': sql, 'format':'csv'})
    for i in range(3):
        try:
            with urlopen(url, timeout=240) as response:
                body = response.read()
            if not body.startswith(b'#Table1'):
                raise RuntimeError(repr(body[:400]))
            return body
        except Exception:
            if i==2: raise
            time.sleep(2**i)

def fetch_table(table):
    columns = COMMON + ',' + EXTRA[table]
    # DR19 SDSS-V APOGEE allStar e_v_rad is commonly zero. Inspect visit e_v_rel later.
    base = BASE if table == 'mwm_boss_allstar' else BASE.replace('e_v_rad>0 AND e_v_rad<=25', 'n_good_rvs>=1')
    tasks = {'counts':f'SELECT COUNT_BIG(*) AS n_rows, COUNT(DISTINCT sdss_id) AS n_sources FROM {table} WHERE {base}',
             'prefilter': f'SELECT {columns} FROM {table} WHERE {base} AND {PREFILTER} ORDER BY spectrum_pk',
             'control': f'SELECT {columns} FROM {table} WHERE {base} AND (ABS(CHECKSUM(CONVERT(VARCHAR(30),sdss_id))) % 251)=37 ORDER BY spectrum_pk'}
    for label, sql in tasks.items():
        p=WORK / f'{table}_{label}.csv'
        body=query(sql)
        p.write_bytes(body)
        print(p.name, len(body), body.count(b'\n')-2, hashlib.sha256(body).hexdigest(),flush=True)

if __name__ == '__main__':
    with ThreadPoolExecutor(2) as pool: list(pool.map(fetch_table,TABLES))
