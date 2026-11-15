"""Reproduce the exploratory Wind/STEREO-A vector-depression selection."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.ndimage import binary_closing, label

PRODUCTS = {'Wind': ('WI_H3-RTN_MFI@0', 'BRTN'), 'STEREO-A': ('STA_L2_MAGPLASMA_1M', 'BFIELDRTN')}
DATES = pd.date_range('2023-05-01', '2023-12-31', freq='8D', tz='UTC')


def load_chunk(args):
    sc, t, cache = args
    dest = cache / (sc + '_' + t.strftime('%Y%m%d') + '.csv')
    if not dest.exists():
        product, parameter = PRODUCTS[sc]
        query = {'id':product, 'parameters':parameter, 'time.min':t.strftime('%Y-%m-%dT00:00:00Z'),
                 'time.max':(t+pd.Timedelta('8D')).strftime('%Y-%m-%dT00:00:00Z'), 'format':'csv'}
        with urlopen('https://cdaweb.gsfc.nasa.gov/hapi/data?'+urlencode(query), timeout=100) as response:
            raw = response.read()
        if not raw.startswith(b'2023'):
            raise ValueError(f'Unexpected response for {sc} on {t}')
        dest.write_bytes(raw)
    frame = pd.read_csv(dest, names=['UTC','R','T','N'])
    frame.index = pd.to_datetime(frame.pop('UTC'), utc=True)
    return frame.mask(frame.abs() > 500)


def candidates(sc, cache):
    with ThreadPoolExecutor(max_workers=4) as pool:
        frames = list(pool.map(load_chunk, [(sc,t,cache) for t in DATES]))
    grid = pd.date_range('2023-05-01','2024-01-04',freq='5min',tz='UTC')
    vectors = pd.concat(frames).resample('5min').mean().reindex(grid)
    magnitude = np.sqrt(vectors.pow(2).sum(axis=1,min_count=3))
    background = magnitude.rolling(288, center=True, min_periods=205).median()
    ratio = magnitude / background
    mask = (ratio < .45) & (magnitude < 4)
    features, count = label(binary_closing(mask, np.ones(3)))
    out = []
    for idx in range(1, count+1):
        sample = ratio[features == idx]
        if len(sample) >= 5 and magnitude.loc[sample.index].min() < 2:
            out.append((sample.idxmin(), len(sample)*5, round(sample.min(),2)))
    return out


def run(cache):
    wind, stereo = [candidates(sc, cache) for sc in PRODUCTS]
    rows = []
    for w, wd, wr in wind:
        for s, sd, sr in stereo:
            dt = (w-s).total_seconds()/3600
            if -4 < dt < 10:
                rows.append([w,s,dt,wd,sd,wr,sr,min(wd,sd)/(max(wr,sr)+.05)])
    cols = ['Wind_min_utc','STEREO_min_utc','lag_hours','Wind_duration_min','STEREO_duration_min',
            'Wind_vector_ratio','STEREO_vector_ratio','ranking_value']
    return pd.DataFrame(rows,columns=cols).sort_values('ranking_value',ascending=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, required=True, help='Directory for the 62 archive chunks')
    args = parser.parse_args()
    args.cache.mkdir(exist_ok=True, parents=True)
    print(run(args.cache).to_csv(index=False), end='')
