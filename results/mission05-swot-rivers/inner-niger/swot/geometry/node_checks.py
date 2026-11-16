"""Fixed-node location and WSE consistency at three RiverSP-v2.0 and D target reaches."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from fetch_geom import out
index=json.loads((out/'sword_node_crosswalk.json').read_text())['selected']
summary=[]
event=[]
for item in index:
    rid=item['v16_reach_id']; nid=item['v16_node_id']
    for version in ['2.0','D']:
        x=pd.read_csv(out/f'node_{nid}_{version}.csv')
        assert not x.duplicated(['cycle_id','pass_id']).any()
        y=pd.read_csv(out/f'reach_{rid}_{version}.csv')
        z=x.merge(y[['cycle_id','pass_id','reach_q','obs_frac_n','partial_f','wse']],on=['cycle_id','pass_id'],suffixes=('_node','_reach'),validate='one_to_one')
        z['prior_distance_m']=np.hypot((z.lat-float(item['lat']))*111195,(z.lon-float(item['lon']))*(111195*np.cos(np.deg2rad(float(item['lat'])))))
        node_ok=(z.node_q<=1)&(z.n_good_pix>0)&(z.wse_node>-1000)&(z.wse_u>0)
        reach_ok=(z.reach_q<=1)&(z.obs_frac_n>=.7)&(z.partial_f==0)&(z.wse_reach>-1000)
        z['residual']=z.wse_node-z.wse_reach
        for orbit,grp in z[node_ok & reach_ok].groupby('pass_id'):
            summary.append(dict(reach_v16=rid,node_v16=nid,collection=version,pass_id=orbit,n=len(grp),median_wse_node_minus_reach_m=grp.residual.median(),median_abs_wse_node_minus_reach_m=grp.residual.abs().median(),median_reported_minus_prior_distance_m=grp.prior_distance_m.median(),p95_reported_minus_prior_distance_m=grp.prior_distance_m.quantile(.95),median_cross_track_m=grp.xtrk_dist.median()))
        for row in z[z.time_str.str.startswith(('2024-12-08','2024-12-15','2024-12-29'),na=False)].to_dict('records'):
            event.append(dict(reach_v16=rid,node_v16=nid,collection=version,time_str=row['time_str'],pass_id=row['pass_id'],node_q=row['node_q'],wse_node_m=row['wse_node'],wse_u_m=row['wse_u'],reach_q=row['reach_q'],wse_reach_m=row['wse_reach'],obs_frac=row['obs_frac_n'],partial_f=row['partial_f'],width_node_m=row['width'],prior_distance_m=row['prior_distance_m'],node_dist_m=row['node_dist'],cross_track_m=row['xtrk_dist']))
pd.DataFrame(summary).to_csv(out/'node_quality_by_track.csv',index=False)
d=pd.DataFrame(event); d.to_csv(out/'node_event.csv',index=False)
s=d[(d.collection=='D')&(d.reach_v16.str.startswith('142777')) & (d.node_q<=1)].copy()
print(s.groupby(['reach_v16',s.time_str.str[:10]]).agg(N=('wse_node_m','size'),median_fixed_node=('wse_node_m','median'),min_m=('wse_node_m','min'),max_m=('wse_node_m','max'),median_u=('wse_u_m','median'),max_centroid_offset_m=('prior_distance_m','max')).to_string())
print(s.groupby(s.time_str.str[:10]).node_q.value_counts())
