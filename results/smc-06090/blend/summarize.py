"""Tabulate pointwise held-block predictions and a paired observing-season bootstrap."""
import json
from pathlib import Path
import numpy as np
from compare import get_data,fold_id,solve_model,HERE

def main():
    records=sum((json.loads((HERE/name).read_text())['fitted'] for name in ['basic.json','blend.json']),[])
    models=['static','two_clocks','blend','apsidal']
    D=get_data();folds=fold_id(D[:,0]);out={}
    for model in models:
        loss=np.empty(len(D));sq=np.empty(len(D))
        for fit in (x for x in records if x['model']==model and x['fold']!='all'):
            keep=folds==fit['fold']
            dx=solve_model(fit['params'],D[keep],model)-D[keep,4]
            r=dx/D[keep,5]
            loss[keep]=2*1.8**2*(np.sqrt(1+(r/1.8)**2)-1);sq[keep]=dx**2
        out[model]={'loss':loss,'sq':sq,'by_block':[float(loss[folds==i].sum()) for i in range(6)]}
    # Seasons group observations across gaps longer than 60 days, before resampling.
    order=np.argsort(D[:,0]);seasons=np.zeros(len(D),int)
    seasons[order]=np.cumsum(np.r_[0,np.diff(D[order,0])>60])
    rng=np.random.default_rng(77714)
    seasons_unique=np.unique(seasons)
    size=np.bincount(seasons);draws=rng.integers(0,len(size),(10000,len(size)))
    summary={}
    for m in models:
        diff=out[m]['loss']-out['apsidal']['loss']
        sums=np.bincount(seasons,weights=diff)
        bootstrap=sums[draws].sum(axis=1)/size[draws].sum(axis=1)
        summary[m]={'total_test_loss':float(out[m]['loss'].sum()),
                    'pooled_test_rms_flux':float(np.sqrt(out[m]['sq'].mean())),
                    'test_loss_by_time_block':out[m]['by_block'],
                    'n_time_blocks_with_antiphased_lower':sum(int(out[m]['by_block'][i]>out['apsidal']['by_block'][i]) for i in range(6)),
                    'mean_excess_loss_per_observation_vs_antiphased':float(diff.mean()),
                    'observing_season_gap_days':60,'season_count':len(size),
                    'seasons_with_antiphased_lower':int((sums>0).sum()),
                    'season_resample_interval_mean_loss_difference':np.quantile(bootstrap,[.025,.975]).tolist()}
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
