"""Lagged cross-image template with flexible slow extrinsic spline; seasonal holdouts.

A(t) = external_spline_A(t) + alpha * smoothed_B(t - delta).
Each B prediction uses >=3 actual B nights within +/- 3 bandwidths.
Compare lags on identical A test nights to avoid sampling-opportunity bias.
The B template uses ALL years, including B dates corresponding to held-out A
seasons: this is cross-image reconstruction, NOT prospective prediction.
The fitted amplitude may be negative for null comparators (unphysical for
macro-lensing). This is a *competing statistical model* rather than a lens.
"""
import argparse
import csv
import json
import numpy as np
from scipy.interpolate import BSpline
from gaussian_process import ROOT, TARGET

START, END = 58190., 61000.


def curves():
    data = {"A": [], "B": []}
    with (ROOT / "ztf-nights.csv").open(newline="") as f:
        for r in csv.DictReader(f):
            if r["system"] == TARGET and r["band"] == "zg":
                data[r["image"]].append([float(r["mjd"]),float(r["magnitude"]),
                                           float(r["median_catalog_error"])])
    return {im:np.array(sorted(v)).T for im,v in data.items()}


def smooth_B(B, target_times, width):
    t,m,e = B
    dt = target_times[:,None] - t[None,:]
    w = np.exp(-.5*(dt/width)**2)/(e[None,:]**2+.035**2)
    w[abs(dt) > 3*width] = 0
    value = (w@m) / np.maximum(w.sum(axis=1), 1e-99)
    return value, ((w > 0).sum(axis=1) >= 3)


def design(A,B,delay,spacing,width):
    t = A[0]
    internal = np.arange(START+spacing,END,spacing)
    knots = np.r_[np.repeat(START,4),internal,np.repeat(END,4)]
    slow = BSpline.design_matrix(np.clip(t,START,END),knots,3,extrapolate=False).toarray()
    fast,valid = smooth_B(B,t-delay,width)
    X = np.c_[slow,fast-np.median(B[1])]
    return X,valid,slow.shape[1]


def fit_model(X,y,e,nslow,ridge=4.):
    # Shrink slow extrinsic coefficients to training median to stabilize holdouts.
    center=np.median(y)
    scale=1/np.sqrt(e**2+.06**2)
    Z=X*scale[:,None]
    reg=np.eye(X.shape[1])*1e-6
    reg[:nslow,:nslow] += np.eye(nslow)*3.
    reg[-1,-1] += .01
    for i in range(1,nslow-1):
        d=np.zeros(X.shape[1]);d[i-1:i+2]=[1,-2,1]
        reg+=ridge*np.outer(d,d)
    co=np.linalg.solve(Z.T@Z+reg,Z.T@((y-center)*scale))
    co[:nslow]+=center
    return co


def validate(spacing=550,width=20,candidates=(0,210,365)):
    series=curves()
    A,B=series["A"],series["B"]
    t,y,e=A
    mm={d:design(A,B,d,spacing,width) for d in candidates}
    all_support=np.logical_and.reduce([v[1] for v in mm.values()])
    seasons=np.floor((t-58200)/365.25).astype(int)
    output=[]
    for fold in sorted(set(seasons)):
        test=all_support & (seasons==fold)
        if test.sum()<6:continue
        row={"season":int(fold),"n":int(test.sum())}
        for d,(X,valid,nslow) in mm.items():
            train=(seasons!=fold)&valid
            if train.sum()<50: continue
            co=fit_model(X[train],y[train],e[train],nslow)
            err=(X[test]@co-y[test])
            row[str(d)]={"squared_error_sum":round(float(err@err),7),
                         "rms":round(float(np.sqrt(np.mean(err**2))),5),
                         "alpha":round(float(co[-1]),4)}
        output.append(row)
    total=sum(r["n"] for r in output)
    pooled={str(d):round(float(np.sqrt(sum(r[str(d)]["squared_error_sum"]
                     for r in output)/total)),5) for d in candidates}
    return {"spacing":spacing,"bandwidth":width,"holdout_A_seasons":output,
            "n_common_A_nights":int(total),"pooled_holdout_rms":pooled}


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--quick",action="store_true")
    args=parser.parse_args()
    for spacing in ((550,) if args.quick else (350,550,850)):
        print(json.dumps(validate(spacing),sort_keys=True))
