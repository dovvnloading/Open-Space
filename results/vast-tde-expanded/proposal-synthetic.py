"""Prospective proposal only: fixed metadata grouping plus synthetic arithmetic.
No real radio flux is loaded, and no discovery threshold is defined.
"""
from datetime import datetime
import numpy as np

def sessions(rows):
    if len({r['image_id'] for r in rows})!=len(rows):raise ValueError('duplicate image ID')
    ordered=sorted(rows,key=lambda x:x['time']); result=[]
    for row in ordered:
        t=datetime.fromisoformat(row['time'])
        if not result or (t-datetime.fromisoformat(result[-1][-1]['time'])).total_seconds()>86400:result.append([])
        result[-1].append(row)
    if len(result)<4:raise ValueError('insufficient distinct sessions')
    return result[:-2],result[-2:]

def diagnostic(early,late):
    # Inputs are independent session-level values and 1-sigma errors; additive
    # calibration/forced-fit covariance is NOT modeled in this toy.
    e=np.asarray(early,float);h=np.asarray(late,float)
    if e.ndim!=2 or h.shape!=(2,2) or e.shape[1]!=2 or len(e)<2 or not np.isfinite(e).all() or not np.isfinite(h).all() or np.any(e[:,1]<=0) or np.any(h[:,1]<=0):raise ValueError('invalid toy sample')
    weights=1/e[:,1]**2; baseline=np.sum(e[:,0]*weights)/sum(weights);var=1/sum(weights)
    return (h[:,0]-baseline)/np.sqrt(h[:,1]**2+var)

if __name__=='__main__':
    rows=[{'image_id':str(i),'time':t} for i,t in enumerate(['2025-01-01T00:00:00','2025-02-01T00:00:00','2025-03-01T00:00:00','2025-03-01T00:15:00','2025-04-01T00:00:00','2025-04-01T00:15:00'])]
    early,hold=sessions(rows);assert list(map(len,hold))==[2,2] and len(early)==2
    assert np.allclose(diagnostic([(1,.1),(1,.1)],[(1,.1),(1,.1)]),[0,0])
    assert np.all(diagnostic([(1,.1),(1,.1)],[(1.5,.1),(1.5,.1)])>0)
    try:sessions(rows+[rows[0]])
    except ValueError:pass
    else:raise AssertionError('duplicate accepted')
    print('synthetic grouping and counterfactual checks passed; no actual measurements read')
