"""Measure independent-pipeline RV disagreements in the DESI velocity tail and controls."""
import gzip, hashlib, json
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import beta
ROOT=Path(__file__).resolve().parent
C=299792.458
def read(name):
    return pd.read_csv(ROOT/name,dtype={'targetid':str,'source_id':str})
def enrich(d):
    l,b=np.radians(d.glon),np.radians(d.glat)
    d['v_gsr_los']=d.vrad+(11.1*np.cos(l)+245.6*np.sin(l))*np.cos(b)+7.25*np.sin(b)
    # Relativistic redshift relation is negligible (~0.4 km/s at 500 km/s).
    d['delta_rv']=abs(d.vrad-C*d.rr_z)
    return d
def rate(d):
    k=int((d.delta_rv>60).sum());n=len(d)
    lo=beta.ppf(.025,k,n-k+1) if k else 0
    hi=beta.ppf(.975,k+1,n-k) if k<n else 1
    return dict(rows=n,disagreements=k,rate=k/n,binomial_95pct=[lo,hi])
def main():
    tail=enrich(read('losrv.csv'))
    tail=tail.query("rvs_warn==0 and vrad_err>0 and vrad_err<20 and rr_spectype=='STAR' and abs(v_gsr_los)>580").copy()
    control=enrich(read('control.csv.gz'))
    # Select a distinct moderate-velocity match for each tail row within the
    # same observing program, survey when available, and closest log-S/N.
    pool=control[abs(control.v_gsr_los)<300].copy()
    matches=[]
    for _,r in tail.sort_values(['survey','program','sn_r','targetid']).iterrows():
        sub=pool[pool.program==r.program]
        same=sub[sub.survey==r.survey]
        if not same.empty: sub=same
        i=(np.log1p(np.maximum(sub.sn_r,0))-np.log1p(max(r.sn_r,0))).abs().idxmin()
        match=pool.loc[i].copy();match['tail_targetid']=r.targetid;matches.append(match);pool=pool.drop(i)
    m=pd.DataFrame(matches)
    m[['targetid','tail_targetid','survey','program','sn_r','vrad','vrad_err','rr_z','v_gsr_los','delta_rv']].to_csv(ROOT/'matched.csv',index=False)
    report={'definition':'abs(vrad - c * rr_z) > 60 km/s; RVS and Redrock use same spectrum','quality_coadd_rows':int(read('cohort.csv').n.iloc[0]),
           'tail_abs_los_gt580': rate(tail),'tail_abs_los_gt1000':rate(tail[abs(tail.v_gsr_los)>1000]),
           'random_control':rate(control),'matched_moderate_los_control':rate(m),'coadds_in_tail':len(read('fastlos.csv')),
           'unique_tail_targetid':tail.targetid.nunique(),'unique_random_targetid':control.targetid.nunique()}
    (ROOT/'rates.json').write_text(json.dumps(report,indent=2)+'\n')
    h={}
    for f in sorted(ROOT.glob('*.csv*')):
        data=gzip.decompress(f.read_bytes()) if f.suffix=='.gz' else f.read_bytes()
        h[f.name]={'sha256_uncompressed':hashlib.sha256(data).hexdigest(),'bytes_uncompressed':len(data)}
    (ROOT/'hashes.json').write_text(json.dumps(h,indent=2)+'\n')
if __name__=='__main__':main()
