"""Replay same-orbit same-stage and anniversary comparisons from frozen source rows."""
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis'
DIRE = 14277700101
MOPTI = 14279000471


def visits(reach, version, cover):
    if version == '2.0':
        data = pd.read_csv(ROOT / f'swot/v3/{reach}.csv')
    else:
        data = pd.read_csv(ROOT / f'swot/geometry/reach_{reach}_D.csv')
    mask = (data.reach_q.isin([0,1]) & (data.partial_f==0) &
            (data.obs_frac_n>=cover) & (data.width>0) & (data.width_u>=0) &
            (data.slope>-1) & (data.slope_u>=0) & (data.wse>-1000) & (data.wse_u>0))
    if 'xovr_cal_q' in data:
        mask &= data.xovr_cal_q.isin([0,1])
    data = data[mask].copy()
    data['date'] = pd.to_datetime(data.time_str,utc=True)
    data['season'] = data.date.dt.year-(data.date.dt.month<6)
    return data.set_index(['cycle_id','pass_id'],verify_integrity=True).sort_index()


def hydraulic_loop():
    samples, comparisons = [], []
    for version, floor in product(('2.0','D'), (.8,.75)):
        site = visits(MOPTI,version,floor).join(visits(DIRE,version,floor),lsuffix='_m',rsuffix='_d',how='inner')
        for (cycle,track), row in site.iterrows():
            if abs((row.date_m-row.date_d).total_seconds())>120:
                raise ValueError('Acquisitions differ by more than two minutes')
            samples.append(dict(version=version,floor=floor,cycle=cycle,track=track,
                                time_mopti=row.time_str_m,time_dire=row.time_str_d,
                                mopti_m=row.wse_m,dire_m=row.wse_d,season=row.season_m))
        for (track,season), group in site.groupby(['pass_id','season_m']):
            rise = group[group.date_m.dt.month.between(7,9)]
            fall = group[group.date_m.dt.month.between(1,3)]
            for (_, a), (_, b) in product(rise.iterrows(),fall.iterrows()):
                if abs(b.wse_m-a.wse_m)>.2:
                    continue
                row=dict(version=version,floor=floor,track=track,season=season,
                         rising=a.time_str_m,receding=b.time_str_m,
                         mopti_change_m=b.wse_m-a.wse_m,dire_change_m=b.wse_d-a.wse_d,
                         dire_width_change_m=b.width_d-a.width_d,
                         dire_slope_change_cm_km=(b.slope_d-a.slope_d)*1e5,
                         sigma_dire_change_m=np.hypot(a.wse_u_d,b.wse_u_d),
                         sigma_width_change_m=np.hypot(a.width_u_d,b.width_u_d),
                         sigma_slope_change_cm_km=np.hypot(a.slope_u_d,b.slope_u_d)*1e5,
                         crid_rise=a.crid_d,crid_fall=b.crid_d,
                         min_frac=min(a.obs_frac_n_m,a.obs_frac_n_d,b.obs_frac_n_m,b.obs_frac_n_d))
                comparisons.append(row)
    pd.DataFrame(samples).to_csv(OUT/'concurrent-stages.csv',index=False)
    pd.DataFrame(comparisons).to_csv(OUT/'matched-limbs.csv',index=False)


def anniversary_pairs():
    ids = [p.stem for p in sorted((ROOT/'swot/v3').glob('[0-9]*.csv'))]
    rows, summaries = [], []
    for version,floor,reach,stage_limit in product(('2.0','D'),(.8,.75),ids,(None,.20)):
        data=visits(reach,version,floor).reset_index()
        for season in (2023,2024):
            matches=[]
            for track, group in data.groupby('pass_id'):
                left=group[group.season==season]
                right=group[group.season==season+1]
                if left.empty or right.empty:
                    continue
                # Match by anniversary at this pass, with no reused visit for this year pair.
                # Ineligible cells cost more than the total of all valid cells.
                anniv=left.date+pd.DateOffset(years=1)
                days=np.abs((anniv.to_numpy()[:,None]-right.date.to_numpy()[None,:]))
                days=np.array([[v.total_seconds()/86400 for v in line] for line in days])
                valid=days<=14
                if stage_limit is not None:
                    valid &= abs(left.wse.to_numpy()[:,None]-right.wse.to_numpy()) <= stage_limit
                ii,jj=linear_sum_assignment(np.where(valid,days,1e6))
                for i,j in zip(ii,jj):
                    if not valid[i,j]:
                        continue
                    a,b=left.iloc[i],right.iloc[j]
                    entry=dict(reach=reach,version=version,floor=floor,stage_limit=stage_limit,
                               season=season,track=track,t0=a.time_str,t1=b.time_str,
                               date_offset_days=days[i,j],wse_delta_m=b.wse-a.wse,
                               width_delta_m=b.width-a.width,slope_delta_cm_km=(b.slope-a.slope)*1e5,
                               crid0=a.crid,crid1=b.crid,frac0=a.obs_frac_n,frac1=b.obs_frac_n)
                    q0=a.get('sos_consensus_q',np.nan); q1=b.get('sos_consensus_q',np.nan)
                    entry['q_delta_m3s']=q1-q0 if min(q0,q1)>0 else np.nan
                    start=pd.Timestamp(f'{season}-06-01',tz='UTC')
                    entry['block42']=(a.date-start).days//42
                    rows.append(entry);matches.append(entry)
            frame=pd.DataFrame(matches)
            for metric in ('width_delta_m','slope_delta_cm_km','q_delta_m3s','wse_delta_m'):
                if (version=='D' and metric=='q_delta_m3s') or (stage_limit and metric=='wse_delta_m'):
                    continue
                counts=pd.Series(dtype=float) if frame.empty else frame.groupby('block42')[metric].median().dropna()
                nonzero=counts[counts!=0]; n=len(nonzero)
                p=binomtest(int((nonzero>0).sum()),n,.5).pvalue if n else np.nan
                summaries.append(dict(reach=reach,version=version,floor=floor,season=season,stage_limit=stage_limit,
                                      metric=metric,n_pairs=0 if frame.empty else frame[metric].notna().sum(),
                                      same_crid_pairs=0 if frame.empty else (frame[metric].notna()&(frame.crid0==frame.crid1)).sum(),
                                      n_nonzero_blocks=n,median_of_block_medians=counts.median(),
                                      nominal_sign_p=p))
    pd.DataFrame(rows).to_csv(OUT/'anniversary-pairs.csv',index=False)
    results=pd.DataFrame(summaries)
    # Primary discovery family: 9 sites x [4 calendar-matched + 3 stage-matched endpoints].
    primary=(results.version=='2.0')&(results.floor==.8)&(results.season==2023)
    results.loc[primary,'primary_bonferroni_p']=(63*results.loc[primary,'nominal_sign_p']).clip(upper=1.)
    results.to_csv(OUT/'anniversary-tests.csv',index=False)


if __name__=='__main__':
    hydraulic_loop()
    anniversary_pairs()
