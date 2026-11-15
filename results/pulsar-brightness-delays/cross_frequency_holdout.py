#!/usr/bin/env python3
"""A locked cross-frequency and time-transfer test of B1822-09, ADDITION to PR #21.

The original published pair slopes used channels 0,3 and 4,7.  Independently
test cross-band pairs 1,6 and 2,5 (no shared arrival-time channels) and 0,7.
Inputs are the immutable 52-visit derived values already in this PR.  This is
not a new observatory or independent epoch sample.
"""
import csv
import hashlib
import json
import math
import random
from pathlib import Path

INPUT=Path(__file__).with_name("visits.csv")
EXPECTED_BLOB="ef04158ff595e5cc377658a6f94175d78daf58be"
K=("0_3","4_7","0_7","1_6","2_5")
N=52

def read():
    raw=INPUT.read_bytes()
    blob=hashlib.sha1(b"blob "+str(len(raw)).encode()+b"\0"+raw).hexdigest()
    assert blob==EXPECTED_BLOB, ("Source data unexpectedly changed",blob)
    records=[]
    for x in csv.DictReader(raw.decode().splitlines()):
        if x["pulsar"]=="J1825-0935":
            records.append({k:(v if k in ("pulsar","observation") else float(v)) for k,v in x.items()})
    records.sort(key=lambda x:x["mjd"])
    assert len(records)==N
    return records

def linear(a,b):
    n=len(b)
    z=[r[:]+[v] for r,v in zip(a,b)]
    for j in range(n):
        k=max(range(j,n),key=lambda t:abs(z[t][j]))
        z[j],z[k]=z[k],z[j]
        assert abs(z[j][j])>1e-14
        scale=z[j][j]
        for c in range(j,n+1): z[j][c]/=scale
        for i in range(n):
            if i!=j:
                factor=z[i][j]
                for c in range(j,n+1):z[i][c]-=factor*z[j][c]
    return [v[-1] for v in z]

def slope(rows,y,adjust=False):
    t0=sum(r["mjd"] for r in rows)/len(rows)
    f0=sum(r["flux_mjy"] for r in rows)/len(rows)
    x=[]
    for r in rows:
        t=(r["mjd"]-t0)/365.25
        a=[1.,t,r["flux_mjy"]-f0]
        if adjust:
            a += [r["duration_s"]/500,math.sin(2*math.pi*t),math.cos(2*math.pi*t),
                  (r["f0_mhz"]+r["f7_mhz"])/2000,
                  (r["f7_mhz"]-r["f0_mhz"])/700]
        x.append(a)
    n=len(x[0]); a=[[0.]*n for _ in range(n)];b=[0.]*n
    for row,y0 in zip(x,y):
        for j in range(n):
            b[j]+=row[j]*y0
            for k in range(n):a[j][k]+=row[j]*row[k]
    return linear(a,b)[2]

def values(rows,pair):
    return [r["relative_dm_"+pair] for r in rows]

def xsum(r,i,j):
    return 1/r["f%d_mhz"%i]**2+1/r["f%d_mhz"%j]**2

def predict(rows,pair):
    i,j=map(int,pair.split("_"))
    out=[]
    for r in rows:
        lo=xsum(r,0,3);hi=xsum(r,4,7)
        weight=(xsum(r,i,j)-lo)/(hi-lo)
        out.append(r["relative_dm_0_3"]+weight*(r["relative_dm_4_7"]-r["relative_dm_0_3"]))
    return out

def contrast(rows,name):
    if name=="upper_minus_lower":
        return [r["relative_dm_4_7"]-r["relative_dm_0_3"] for r in rows]
    if name=="inner_minus_outer":
        return [r["relative_dm_1_6"]-r["relative_dm_0_7"] for r in rows]
    if name=="center_minus_outer":
        return [r["relative_dm_2_5"]-r["relative_dm_0_7"] for r in rows]
    pair=name.split("error_")[1]
    return [x-y for x,y in zip(values(rows,pair),predict(rows,pair))]

def analyze(rows):
    fit={k:slope(rows,values(rows,k)) for k in K}
    holdout={}
    for k in K[2:]:
        estimated=slope(rows,predict(rows,k))
        observed=fit[k]
        holdout[k]={"observed":round(observed,9),"predicted":round(estimated,9),
                    "error":round(observed-estimated,9),
                    "error_percent":round(100*abs(observed-estimated)/abs(observed),3)}
    rng=random.Random(1715)
    check={}
    for name in ("upper_minus_lower","inner_minus_outer","center_minus_outer",
                 "error_0_7","error_1_6","error_2_5"):
        y=contrast(rows,name)
        observed=slope(rows,y)
        boot=[]
        for _ in range(4000):
            ids=[]
            for j in range(13):
                off=rng.randrange(N)
                ids.extend((off+k)%N for k in range(4))
            ids=ids[:N]
            boot.append(slope([rows[i] for i in ids],[y[i] for i in ids]))
        boot.sort()
        shifts=[]
        for k in range(1,N):
            shifted=[dict(row,flux_mjy=rows[(i+k)%N]["flux_mjy"]) for i,row in enumerate(rows)]
            shifts.append(slope(shifted,y))
        check[name]={"slope":round(observed,9),
                     "block95":[round(boot[99],9),round(boot[3899],9)],
                     "nonzero_circular_shifts_at_least_as_extreme":sum(abs(v)>=abs(observed) for v in shifts),
                     "first26":round(slope(rows[:26],y[:26]),9),
                     "last26":round(slope(rows[26:],y[26:]),9)}
    # Freeze early training coefficients, challenge them with later unseen channels.
    low=slope(rows[:26],values(rows[:26],"0_3"))
    high=slope(rows[:26],values(rows[:26],"4_7"))
    temporal={}
    for k in ("1_6","2_5"):
        i,j=map(int,k.split("_"))
        w=[(xsum(r,i,j)-xsum(r,0,3))/(xsum(r,4,7)-xsum(r,0,3)) for r in rows[26:]]
        pred=low+(high-low)*sum(w)/len(w)
        actual=slope(rows[26:],values(rows[26:],k))
        temporal[k]={"early_frozen_prediction":round(pred,9),
                     "late_measured":round(actual,9),
                     "late_error_pct":round(100*abs(actual-pred)/abs(actual),3)}
    trimmed=sorted(rows,key=lambda r:r["flux_mjy"])[5:-5]
    extra={"trimmed_n":len(trimmed),
           "trimmed_upper_minus_lower":round(slope(trimmed,contrast(trimmed,"upper_minus_lower")),9),
           "trimmed_inner_minus_outer":round(slope(trimmed,contrast(trimmed,"inner_minus_outer")),9),
           "covariate_upper_minus_lower":round(slope(rows,contrast(rows,"upper_minus_lower"),True),9),
           "covariate_inner_minus_outer":round(slope(rows,contrast(rows,"inner_minus_outer"),True),9)}
    loo=[slope(rows[:i]+rows[i+1:],contrast(rows[:i]+rows[i+1:],"inner_minus_outer")) for i in range(N)]
    extra["leave_one_epoch_out_range"]=[round(min(loo),9),round(max(loo),9)]
    result={"sample":N,"target":"B1822-09","regression":"intercept + years + simultaneous flux",
            "pairs":{k:round(v,9) for k,v in fit.items()},
            "heldout_channel_pairs":holdout,"controls":check,
            "frozen_early_to_late":temporal,"robustness":extra}
    return result

def frozen_checks(o):
    assert o["sample"]==52
    assert abs(o["pairs"]["0_3"]+.011850671)<.000005
    assert abs(o["pairs"]["4_7"]+.024561510)<.000005
    for pair,tolerance in (("0_7",3),("1_6",5),("2_5",8)):
        assert o["heldout_channel_pairs"][pair]["error_percent"]<tolerance
        lo,hi=o["controls"]["error_"+pair]["block95"]
        assert lo<0<hi,(pair,lo,hi)
    for name in ("upper_minus_lower","inner_minus_outer","center_minus_outer"):
        a=o["controls"][name]
        lo,hi=a["block95"]
        assert lo<hi<0
        assert a["nonzero_circular_shifts_at_least_as_extreme"]==0
    assert o["robustness"]["covariate_inner_minus_outer"]<0
    assert o["robustness"]["trimmed_inner_minus_outer"]<0
    assert o["robustness"]["leave_one_epoch_out_range"][1]<0
    for k in ("1_6","2_5"):
        assert o["frozen_early_to_late"][k]["late_error_pct"]<12

if __name__=="__main__":
    result=analyze(read())
    frozen_checks(result)
    print(json.dumps(result,indent=2,sort_keys=True))
