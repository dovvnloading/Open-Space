"""Audit annual OPIDIN peak relationships and the fixed-node checks on the paired dates."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t

ROOT = Path(__file__).resolve().parents[1]


def annual_stages():
    s = pd.DataFrame(json.loads((ROOT/'akka/opidin_maxima.json').read_text())['data'])
    s = s[s.year>=1955].copy()
    rows = []
    # These fits concern annual crests, which occur on different days at each station.
    for station in ['mopti','akka']:
        for start,stop in [(1955,2017),(1955,2023),(1982,2023),(2000,2023)]:
            predictor = (s['level_'+station]-500)/100
            for degree in [1,2]:
                X = np.vander(predictor,degree+1,increasing=True)
                mask = s.year.between(start,stop)
                x=X[mask]; y=s.loc[mask,'level_dire'].to_numpy()
                beta = np.linalg.lstsq(x,y,rcond=None)[0]
                pred=X@beta
                dof = len(y)-len(beta)
                mse = np.sum((y-x@beta)**2)/dof
                se = np.sqrt(mse*(1+np.einsum('ij,jk,ik->i',X,np.linalg.inv(x.T@x),X)))
                for i,(_,r) in enumerate(s.iterrows()):
                    rows.append(dict(predictor=station,model_degree=degree,start=start,end=stop,year=r.year,
                                     predictor_cm=r['level_'+station],dire_cm=r.level_dire,predicted_cm=pred[i],
                                     residual_cm=r.level_dire-pred[i],
                                     nominal_lo95_cm=pred[i]-t.ppf(.975,dof)*se[i],
                                     nominal_hi95_cm=pred[i]+t.ppf(.975,dof)*se[i]))
    pd.DataFrame(rows).to_csv(ROOT/'analysis/annual-stage-fits.csv',index=False)


def fixed_nodes():
    matched=pd.read_csv(ROOT/'analysis/matched-limbs.csv')
    rows=[]
    for path in sorted((ROOT/'swot/geometry').glob('node_1427770010*_D.csv')):
        x=pd.read_csv(path).set_index(['pass_id','cycle_id'])
        for _,r in matched.iterrows():
            starts=pd.to_datetime(x.time_str,errors='coerce',utc=True).dt.date
            left=x[starts.eq(pd.Timestamp(r.rising).date())]
            right=x[starts.eq(pd.Timestamp(r.receding).date())]
            if len(left)!=1 or len(right)!=1:
                raise ValueError('Missing or duplicate fixed-node visit')
            a=left.iloc[0]; b=right.iloc[0]
            ok=all(z.node_q in (0,1) and z.wse_u>0 and z.n_good_pix>0 for z in (a,b))
            rows.append(dict(prior_node_v16=path.name.split('_')[1],season=r.season,rising=a.time_str,
                             receding=b.time_str,q_rising=a.node_q,q_receding=b.node_q,
                             wse_rising=a.wse,wse_receding=b.wse,usable=ok,change_m=b.wse-a.wse,
                             rss_u_m=np.hypot(a.wse_u,b.wse_u)))
    pd.DataFrame(rows).to_csv(ROOT/'analysis/matched-nodes.csv',index=False)


if __name__=='__main__':
    annual_stages()
    fixed_nodes()
