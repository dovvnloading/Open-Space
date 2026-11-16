"""Compare matched cycle/pass geometry after crosswalking SWORD reaches."""
import pandas as pd
from fetch_geom import reaches, out
import json
def available(d,fr):
 return ((d.reach_q<=1)&(d.obs_frac_n>=fr)&(d.partial_f==0)&(d.wse>-1e6)&(d.width>0)&(d.slope>-1.0)&(d.wse_u>=0)&(d.slope_u>=0)&(d.width_u>=0))
def read(ID,vr):
 d=pd.read_csv(out/f'reach_{ID}_{vr}.csv')
 assert not d.duplicated(['cycle_id','pass_id']).any()
 return d
rows=[]
for rid in reaches:
 a=read(rid,'2.0');b=read(rid,'D')
 for f in [.8,.7]:
  aa=available(a,f)
  bb=available(b,f)
  a['usable']=aa;b['usable']=bb
  c=a.merge(b,on=['cycle_id','pass_id'],suffixes=('_c','_d'),validate='one_to_one')
  good=c[c.usable_c&c.usable_d]
  row=dict(reach_v16=rid,reach_v17=str(b.reach_id.iloc[0]),obs_frac_floor=f,n_vC=len(a),eligible_vC=int(aa.sum()),n_vD=len(b),eligible_vD=int(bb.sum()),paired_orbits=len(c),paired_both_eligible=len(good),paired_c_only=int((c.usable_c&~c.usable_d).sum()),paired_d_only=int((~c.usable_c&c.usable_d).sum()),versionD_eligible_since_may2025=int((bb & b.time_str.str.match('2025-(0[5-9]|1)|2026')).sum()))
  for var,units,factor in [('wse','m',1),('width','m',1),('slope','cm_per_km',1e5)]:
   delta=(good[var+'_d']-good[var+'_c'])*factor
   row.update({f'median_D_minus_C_{var}_{units}': delta.median(),f'median_abs_D_minus_C_{var}_{units}':delta.abs().median()})
  rows.append(row)
res=pd.DataFrame(rows);res.to_csv(out/'release_sensitivity.csv',index=False)
print(res.to_string(index=False))
