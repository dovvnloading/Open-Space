"""Compare the archived Collins et al. run 198 with measured HiRISE fields.

Run from the result directory; see MODEL.md. Distances in metres.
"""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import rankdata
from analyze import cluster_rows, equal_count

MODEL_SHA256 = 'cf35f4b28fc9df2f7653c92bc78d05e961aeaa76d1d52040b4363767788858a4'


def pair_gap(x, y, same):
    """Sum |x-y| with prefix sums to avoid quadratic memory on rich fields."""
    y = np.sort(y)
    prefix = np.r_[0., np.cumsum(y)]
    k = np.searchsorted(y, x)
    total = x*(2*k-len(y)) + prefix[-1] - 2*prefix[k]
    return total.sum() / (len(x) * (len(y) - int(same)))


def measure(xy, d):
    n = len(d)
    centered = xy - xy.mean(axis=0)
    axes = centered @ np.linalg.svd(centered, full_matrices=False)[2].T
    along = axes[:, 0]
    if along[np.argmax(abs(along))] < 0:
        along = -along
    labels = equal_count(along, min(4, n//4))
    dist, near = cKDTree(xy).query(xy, k=3)
    if np.any(dist[:,1] == 0):
        raise ValueError('Duplicate crater centers')
    nn = near[:,1]
    # Match the observational implementation's stable input-index tie break.
    for i in np.flatnonzero(np.isclose(dist[:,1], dist[:,2], atol=1e-10, rtol=0)):
        ds = ((xy-xy[i])**2).sum(axis=1)
        ds[i] = np.inf
        nn[i] = np.argmin(ds)
    marks = (rankdata(d)-1)/(n-1)
    groups = [marks[labels==i] for i in np.unique(labels)]
    ref = np.array([[pair_gap(a,b,i==j) for j,b in enumerate(groups)]
                    for i,a in enumerate(groups)])
    gap = np.abs(marks-marks[nn]).mean()
    expected = ref[labels,labels[nn]].mean()
    return dict(n=n, effective_diameter=float(np.cbrt((d**3).sum())),
                extent=float(np.sqrt((centered**2).sum()/n)), observed=gap,
                expected=expected, difference=expected-gap)


def model_sites(path):
    if hashlib.sha256(path.read_bytes()).hexdigest() != MODEL_SHA256:
        raise ValueError('Unexpected model archive checksum')
    with zipfile.ZipFile(path) as z:
        member = next(n for n in z.namelist() if n.endswith('/crater-data/craters-198.csv'))
        with z.open(member) as stream:
            runs = {}
            for row in csv.DictReader(io.TextIOWrapper(stream)):
                runs.setdefault(int(row['ID']), []).append(tuple(float(row[k]) for k in ('x','y','r')))
        for key, rows in runs.items():
            a = np.array(rows)
            yield key, a[:,:2], 2*a[:,2]


def balanced(o,m,keys,rng):
    """With-replacement local comparison, each observed site gets equal weight."""
    def cov(d):
        return np.log([[r[k] for k in keys] for r in d])
    oc, mc = cov(o), cov(m)
    scale = np.subtract(*np.percentile(mc,[75,25],axis=0))
    od = np.array([r['difference'] for r in o])
    md = np.array([r['difference'] for r in m])
    near = cKDTree(mc/scale).query(oc/scale, k=10)[1]
    weights = np.bincount(near.ravel(), minlength=len(m))/(len(o)*near.shape[1])
    paired = (od - md[near].mean(axis=1))
    # Site bootstrap for the empirical observations, with the simulated library fixed.
    ci = np.quantile(rng.choice(paired, (10000,len(o)), replace=True).mean(axis=1), [.025,.975])
    # Also vary the finite simulation library. Resampling a model site twice gives
    # two seats in a ten-neighbor match. Reuse the sorted distance ordering.
    candidates = cKDTree(mc/scale).query(oc/scale, k=min(len(m),100))[1]
    boot = []
    for _ in range(1999):
        draw = rng.integers(len(m), size=len(m))
        counts = np.bincount(draw, minlength=len(m))[candidates]
        cumul = counts.cumsum(axis=1)
        if np.min(cumul[:,-1]) < 10:
            raise ValueError('Too few candidate neighbors in bootstrap')
        seats = np.clip(10 - cumul + counts, 0, counts)
        ref = (seats * md[candidates]).sum(axis=1)/10
        boot.append(rng.choice(od-ref, size=len(o)).mean())
    supported = np.max(np.abs(oc[:,None,:]-mc[near]), axis=(1,2)) <= np.log(1.5)
    return dict(observed_mean=od.mean(), matched_model_mean=np.dot(weights,md),
                ci95_both_site_libraries=list(np.quantile(boot,[.025,.975])),
                close_match_sites=int(supported.sum()),
                close_match_mean_difference=paired[supported].mean(),
                observed_minus_model=paired.mean(), ci95_fixed_library=list(ci),
                model_effective_sample_size=1/np.sum(weights**2),
                observed_median_covariates={k:np.median([r[k] for r in o]) for k in keys},
                matched_median_covariates={k:np.median(np.array([r[k] for r in m])[near]) for k in keys},
                max_local_median_abs_log_ratio={k:float(np.max(np.median(np.abs(oc[:,i,None]-mc[:,i][near]),axis=1))) for i,k in enumerate(keys)})


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('observations', type=Path)
    parser.add_argument('model', type=Path)
    parser.add_argument('output', type=Path)
    args=parser.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    rows=[]
    for kind, data in [('HiRISE',((f['HiRiseID'],xy,d) for f,xy,d in cluster_rows(args.observations))),('run198',model_sites(args.model))]:
        for key, xy, d in data:
            for cutoff in [1,1.5,2,3]:
                use = d >= cutoff
                if use.sum() >= 12:
                    rows.append(dict(source=kind,id=key,cutoff=cutoff,**measure(xy[use],d[use])))
    with (args.output/'model-comparison-fields.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    rng=np.random.default_rng(6566345)
    out={}
    for c in [1,1.5,2,3]:
        o=[r for r in rows if r['cutoff']==c and r['source']=='HiRISE']
        m=[r for r in rows if r['cutoff']==c and r['source']=='run198']
        out[c]=dict(observed_sites=len(o),model_sites=len(m),observed_mean=np.mean([r['difference'] for r in o]),
                    model_mean=np.mean([r['difference'] for r in m]))
        for keys in [('n',),('n','effective_diameter'),('n','effective_diameter','extent')]:
            out[c]['+'.join(keys)]=balanced(o,m,keys,rng)
    (args.output/'model-comparison.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
