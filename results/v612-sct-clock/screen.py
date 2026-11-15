"""Re-run the exploratory clock search on the frozen target list (network required)."""
import argparse
import csv
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

from replay import DATA, clean


def screen_target(table, year):
    good=table[(table.catflags==0)&(table.magerr<.18)&(table.mag>12)].sort_values('hjd')
    good=good.drop_duplicates(['expid','filtercode'])
    parts=[clean(block) for _,block in good.groupby('oid') if len(block)>=12]
    if not parts:
        return None
    rows=pd.concat(parts).sort_values('hjd')
    if 2017<=year<2026:
        eruption=2451544.5+(year-2000)*365.25
        rows=rows[(rows.hjd>eruption+400)|(rows.hjd<eruption-200)]
    if len(rows)<120:
        return None
    division=rows.hjd.quantile(.65)
    early=rows[rows.hjd<division];late=rows[rows.hjd>=division]
    if early.hjd.max()-early.hjd.min()<100 or late.hjd.max()-late.hjd.min()<90:
        return None
    def periodogram(data):
        return LombScargle(data.hjd-2458000,data.res,np.hypot(data.magerr,.06))
    train=periodogram(early);test=periodogram(late)
    f,power=train.autopower(method='fast',minimum_frequency=.30,maximum_frequency=22,samples_per_peak=5)
    use=np.abs(f/1.0027379-np.round(f/1.0027379))>.025
    candidates=[]
    for _ in range(4):
        peak=np.argmax(np.where(use,power,-1));freq=f[peak]
        pred=train.model(late.hjd-2458000,freq)
        corr=np.corrcoef(pred,late.res)[0,1]
        score=np.sign(corr)*np.sqrt(max(0,power[peak]*test.power(freq)))
        candidates.append((score,1/freq,corr));use&=abs(f-freq)>.023
    return (len(rows),*max(candidates))


def run(directory):
    directory.mkdir(parents=True,exist_ok=True)
    with (DATA/'screen.csv').open() as stream:
        for row in csv.DictReader(stream):
            name=row['name'];path=directory/(name.replace(' ','_')+'.csv')
            if not path.exists():
                args={'POS':f"CIRCLE {row['ra_deg']} {row['dec_deg']} 0.00050000",
                      'FORMAT':'CSV','BAD_CATFLAGS_MASK':32768}
                url='https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?'+urlencode(args)
                with urlopen(url,timeout=90) as response:
                    packet=response.read()
                if not packet.startswith(b'oid,'):
                    raise ValueError(f'Unexpected response for {name}')
                path.write_bytes(packet)
            tab=pd.read_csv(path)
            result=screen_target(tab,float(row['catalog_eruption_year'] or 0)) if len(tab) else None
            if result:
                print(name,*result,sep=',',flush=True)


if __name__=='__main__':
    args=argparse.ArgumentParser(description=__doc__)
    args.add_argument('cache',type=Path)
    run(args.parse_args().cache)
