"""Rebuild the reduced tables from public data. Run from any directory."""
import argparse
import concurrent.futures
import gzip
import io
from pathlib import Path
import tempfile
import urllib.parse
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parent
HAPI = 'https://cdaweb.gsfc.nasa.gov/hapi/data'
FIRST = pd.Timestamp('2022-11-01', tz='UTC')
STOP = pd.Timestamp('2024-07-01', tz='UTC')
MONTHS = pd.date_range(FIRST, STOP, freq='MS')
ARCHIVE = {
    'stereo': ('STA_L2_MAGPLASMA_1M', 'BFIELDRTN', ['br','bt','bn']),
    'wind': ('WI_K1-RTN_MFI', 'BRTN', ['br','bt','bn']),
    'ace': ('AC_H0_MFI', 'BGSEc', ['bx_gse','by_gse','bz_gse']),
    'stereoplasma': ('STA_L2_MAGPLASMA_1M', 'BTOTAL,Np,Vp,Tp,Beta', ['b','density','speed','temperature','beta']),
    'windplasma': ('WI_H1_SWE','Proton_V_nonlin,Proton_W_nonlin,Proton_Np_nonlin',['speed','thermal_speed','density']),
}

def read(cache, filename, url):
    dest = cache / filename
    if not dest.exists():
        with urllib.request.urlopen(url, timeout=120) as resp:
            content = resp.read()
        if not content.strip() or content.lstrip().startswith(b'{'):
            raise ValueError(f'No usable data from {url}')
        dest.write_bytes(content)
    return dest.read_bytes()

def hapi(cache, key, start, end, filename):
    dataset, params, names = ARCHIVE[key]
    args = {'id': dataset, 'parameters': params, 'time.min': start.strftime('%Y-%m-%dT00:00:00Z'),
            'time.max': end.strftime('%Y-%m-%dT00:00:00Z'), 'format':'csv'}
    raw = read(cache, filename, HAPI + '?' + urllib.parse.urlencode(args))
    data = pd.read_csv(io.BytesIO(raw), names=['time']+names)
    data.time = pd.to_datetime(data.time, format='ISO8601', utc=True)
    values = data.set_index('time')
    return values.mask(abs(values)>1e8)

def save(frame, name):
    # Deterministic compression and explicit decimal precision.
    with (ROOT/'data'/name).open('wb') as out:
        with gzip.GzipFile(filename='', fileobj=out, mode='wb', mtime=0) as archive:
            archive.write(frame.to_csv(float_format='%.7g', lineterminator='\n').encode())

def build(cache):
    def month(pair):
        start,end=pair
        return hapi(cache, 'stereo', start, end, f'STEREO-A{start:%Y-%m}.csv').resample('10min').median()
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        stereo = pd.concat(list(pool.map(month, zip(MONTHS[:-1], MONTHS[1:]))))
    save(stereo, 'stereo_search_10min.csv.gz')
    start, end = pd.Timestamp('2023-09-28'), pd.Timestamp('2023-09-30')
    wind = hapi(cache, 'wind', pd.Timestamp('2023-09-01'), pd.Timestamp('2023-10-01'), 'Wind2023-09.csv')
    sta = hapi(cache, 'stereo', pd.Timestamp('2023-09-01'), pd.Timestamp('2023-10-01'), 'STEREO-A2023-09.csv')
    def interval(frame):
        return frame.loc['2023-09-28':'2023-09-29']
    w = interval(wind)
    save(w.resample('10min').median(), 'wind_search_10min.csv.gz')
    tables = [w.resample('5min').median().add_prefix('wind_'),
              interval(sta).resample('5min').median().add_prefix('stereo_')]
    for key, filename in [('ace','acesep.csv'),('windplasma','wsep.csv'),('stereoplasma','STEREO-Aseptplasma.csv')]:
        df = hapi(cache,key,start,end,filename)
        if key.endswith('plasma'):
            df=df.mask(df<0)
        tables.append(df.resample('5min').median().add_prefix(key+'_'))
    save(pd.concat(tables,axis=1),'encounter_5min.csv.gz')
    rows=[]
    for name, ident in [('wind','-8'),('stereo','-234')]:
        args={'format':'text','COMMAND':repr(ident),'EPHEM_TYPE':"'VECTORS'",'CENTER':"'500@10'",
              'START_TIME':"'2023-09-29 09:00'",'STOP_TIME':"'2023-09-29 17:00'",'STEP_SIZE':"'1 h'",
              'REF_PLANE':"'ECLIPTIC'",'OUT_UNITS':"'AU-D'",'VEC_TABLE':"'2'",'CSV_FORMAT':"'YES'",'OBJ_DATA':"'NO'"}
        raw=read(cache,'horizons_'+name+'.txt','https://ssd.jpl.nasa.gov/api/horizons.api?'+urllib.parse.urlencode(args))
        body=raw.decode().split('$$SOE')[1].split('$$EOE')[0].strip()
        for line in body.splitlines():
            cols=[s.strip() for s in line.split(',')]
            rows.append([name,cols[1],*map(float,cols[2:5])])
    pd.DataFrame(rows,columns=['spacecraft','calendar_tdb','x_au','y_au','z_au']).to_csv(ROOT/'data'/'positions.csv',index=False)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache',type=Path,help='Keep downloaded science files in this directory')
    opt=parser.parse_args()
    with tempfile.TemporaryDirectory() as tmp:
        cache=opt.cache or Path(tmp)
        cache.mkdir(parents=True,exist_ok=True)
        build(cache)
