"""Evaluate the v2.0 stage–drop relationship on subsequent repeat SWOT visits in D.

Only RiverSP geometry enters this test. SoS v3 does not extend into the holdout.
"""
import argparse
import csv
import hashlib
import json
from io import StringIO
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, theilslopes

BASE = Path(__file__).resolve().parent
API = 'https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries'
# African translation table v17b/AF_ReachIDs_v17b_vs_v16.csv, referenced
# in ../../zambezi/reprocessing/vD_manifest.json. All three boundary flags are zero.
MAP = {'lower':('13217000061','13217000061'),
       'upper':('13217000091','13217000071'),
       'chenal':('13219000041','13219000021')}
FIELDS = ('reach_id,time_str,cycle_id,pass_id,wse,wse_u,reach_q,partial_f,obs_frac_n,'
          'loc_offset,node_dist,xtrk_dist,dark_frac,width,width_u,slope,slope_u,'
          'crid,sword_version,granuleUR,p_lon,p_lat,p_length,p_dist_out')


def acquire():
    manifest = []
    for name,(v16,v17) in MAP.items():
        url = API + '?' + urlencode(dict(collection_name='SWOT_L2_HR_RiverSP_D',feature='Reach',feature_id=v17,
                                      start_time='2023-07-01T00:00:00Z',end_time='2026-11-01T00:00:00Z',
                                      fields=FIELDS,output='csv'))
        with urlopen(url,timeout=60) as response:
            body = response.read()
        data = json.loads(body)
        if data['status'] != '200 OK':
            raise ValueError(data['status'])
        raw = data['results']['csv'].encode()
        assert len(list(csv.DictReader(StringIO(raw.decode())))) == data['hits']
        path = BASE / f'{name}.csv'
        path.write_bytes(raw)
        manifest.append(dict(v16=v16,v17b=v17,file=path.name,url=url,rows=data['hits'],
                             response_sha256=hashlib.sha256(body).hexdigest(),csv_sha256=hashlib.sha256(raw).hexdigest()))
    (BASE/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')


def summarize():
    items = json.loads((BASE/'manifest.json').read_text())
    pair = None
    for rec,(name,ids) in zip(items,MAP.items(),strict=True):
        assert (rec['v16'],rec['v17b']) == ids
        raw = (BASE/rec['file']).read_bytes()
        assert hashlib.sha256(raw).hexdigest()==rec['csv_sha256']
        d = pd.read_csv(StringIO(raw.decode())).replace(-999999999999.,np.nan)
        d['date'] = pd.to_datetime(d.time_str,utc=True,errors='coerce').dt.normalize()
        d['revision'] = d.granuleUR.str.extract(r'_(\d+)\.zip$').astype(int)
        d = d.sort_values('revision').drop_duplicates(['cycle_id','pass_id'],keep='last')
        d = d[d.date.notna()&d.reach_q.isin([0,1])&d.partial_f.eq(0)&(d.obs_frac_n>=.8)
              &d.wse_u.between(0,.30,inclusive='right')&d.wse.notna()]
        d = d[['cycle_id','pass_id','date','wse','crid']].rename(columns={'wse':name,'crid':f'crid_{name}'})
        if pair is None:
            pair = d
        else:
            pair = pair.merge(d,on=['cycle_id','pass_id','date'],validate='one_to_one',how='left' if name=='chenal' else 'inner')
    pair['drop_m'] = pair.upper-pair.lower
    train = pd.read_csv(BASE.parent/'paired-reach.csv').dropna(subset=['chenal_wse'])
    slope,intercept,*_ = theilslopes(train.drop_m,train.chenal_wse)
    pair['predicted_drop_v2_m'] = intercept+slope*pair.chenal
    pair['residual_m'] = pair.drop_m-pair.predicted_drop_v2_m
    pair['period'] = np.where(pair.date>pd.Timestamp('2025-05-02',tz='UTC'),'holdout','overlap')
    pair.to_csv(BASE/'visits.csv',index=False,float_format='%.8g')
    summary = []
    for period,g in pair.dropna(subset=['chenal']).groupby('period'):
        b=theilslopes(g.drop_m,g.chenal).slope if len(g)>=3 else np.nan
        summary.append(dict(period=period,n=len(g),first=g.date.min().date(),last=g.date.max().date(),
                            tracks=','.join(map(str,sorted(g.pass_id.unique()))),r=spearmanr(g.drop_m,g.chenal).statistic,
                            response_slope=b,median_abs_residual_m=g.residual_m.abs().median(),
                            median_signed_residual_m=g.residual_m.median(),train_slope=slope,train_n=len(train)))
    pd.DataFrame(summary).to_csv(BASE/'summary.csv',index=False,float_format='%.8g')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--retrieve',action='store_true')
    if parser.parse_args().retrieve:
        acquire()
    summarize()
