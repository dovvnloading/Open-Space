"""Test achromatic excursions using identical host, test and reference epochs in both bands."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

KEYS=['id','start','n','reference_mask','sample']
BANDS=['f200w','f444w']


def joint(screen):
    parts=[pd.read_csv(screen/b/'contrasts.csv.gz')[KEYS+['score']].rename(columns={'score':b}) for b in BANDS]
    pairs=parts[0].merge(parts[1],on=KEYS,validate='one_to_one')
    pairs['joint_scale']=0.
    for n,grp in pairs.groupby('n'):
        control=grp.loc[grp['sample']=='control',BANDS].sum(axis=1)
        if len(control)<100:raise ValueError('Insufficient comparisons for joint calibration')
        scale=max(np.sqrt(2.),1.4826*np.median(np.abs(control-np.median(control))))
        pairs.loc[grp.index,'joint_scale']=scale
    pairs['joint_score']=pairs[BANDS].sum(axis=1)/pairs.joint_scale
    pairs['passes']=(pairs.f200w*pairs.f444w>0)&(pairs[BANDS].abs().min(axis=1)>=3.)&(pairs.joint_score.abs()>=7.)
    pairs.to_csv(screen/'joint_contrasts.csv.gz',index=False,float_format='%.7g',compression={'method':'gzip','mtime':0})
    base=pairs[(pairs.n==1)&(pairs['sample']=='high_z')].copy()
    base['epoch']=base.start+1
    a,b=[pd.read_csv(screen/ch/'injections.csv.gz')[['id','epoch','injected_ab','injected_score']].rename(columns={'injected_score':ch+'_injected'}) for ch in BANDS]
    inj=a.merge(b,on=['id','epoch','injected_ab'],validate='one_to_one').merge(base[['id','epoch','joint_scale','joint_score']],on=['id','epoch'],validate='many_to_one')
    inj['joint_injected']=(inj.f200w_injected+inj.f444w_injected)/inj.joint_scale
    inj['recovered']=(inj.f200w_injected>=3.)&(inj.f444w_injected>=3.)&(inj.joint_injected>=7.)
    inj.to_csv(screen/'joint_injections.csv.gz',index=False,float_format='%.7g',compression={'method':'gzip','mtime':0})
    inj[inj.joint_score.abs()<7].groupby('injected_ab').recovered.agg(['count','sum','mean']).to_csv(screen/'joint_recovery.csv')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('results',type=Path);args=p.parse_args()
    joint(args.results)
