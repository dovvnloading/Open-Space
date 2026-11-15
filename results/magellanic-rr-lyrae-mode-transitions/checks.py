"""Independent early EROS bands, mixed-frequency tests and period-shift controls."""
import argparse
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
from measure import NAME, PERIODS, RA, DEC, fit, refine
from screen import fetch, model, photo


def eros(cache):
    url = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?" + urlencode({
        "-source":"II/390/xtime", "-c":f"{RA:.7f} {DEC:.7f}", "-c.rs":"2",
        "-out.max":"3000", "-out":"EROS,HJD,Rmag,e_Rmag,Bmag,e_Bmag"})
    rows = [l.split("\t") for l in fetch(url,cache).decode().splitlines() if l.startswith("lm0561m6445")]
    canonical = "\n".join("\t".join(v.strip() for v in r) for r in rows)
    result = {}
    for band, col in [("Red",2),("Blue",4)]:
        good = [(float(r[1]),float(r[col]),float(r[col+1])) for r in rows
                if r[col].strip() and r[col+1].strip() and 0 < float(r[col+1]) <= .16]
        d = np.array(good)
        if np.median(d[:,0]) > 2400000:
            d[:,0] -= 2450000
        p = refine(d, PERIODS)
        # A single radial mode and fixed catalog clocks avoid inflation of weak secondary peaks.
        sol = fit(d,p)
        fixed = fit(d,[p[0],.485254])
        result[band] = dict(n=len(d), first=d[:,0].min(), last=d[:,0].max(), periods=p.tolist(),
                            amps=sol["amp"].tolist(), fixed_early_fund_amps=fixed["amp"].tolist())
    return dict(url=url,rows=len(rows),scientific_rows_sha256=hashlib.sha256(canonical.encode()).hexdigest(),bands=result)


def flux_check(data, p, rng, count):
    # Unresolved independent stars add linearly in flux. Harmonics of each mode are fitted
    # separately, so non-sinusoidal individual waveforms are in the null.
    t, y, e = data.T
    t = t - np.mean(t)
    flux = 10**(-.4*(y-np.median(y)))
    w = 1/(np.maximum(.015,e)*flux*.4*np.log(10))
    X = model(t,p,3)*w[:,None]
    yc = flux*w
    group = np.unique(np.floor(data[:,0]),return_inverse=True)[1]
    projector = np.linalg.pinv(X)
    prediction = X @ (projector @ yc)
    residual = yc-prediction
    freqs = np.array([1/p[0]+1/p[1],1/p[0]-1/p[1]])
    Y = np.column_stack([fn(2*np.pi*t*f)*w for f in freqs for fn in [np.sin,np.cos]])
    Y -= X @ (projector @ Y)
    proj = np.linalg.pinv(Y)
    mixed = proj @ residual
    gain = np.sum((Y @ mixed)**2)
    exceed = 0
    for _ in range(count):
        signs = rng.choice([-1,1],size=group.max()+1)
        noise = residual*signs[group]
        exceed += np.sum((Y @ (proj @ noise))**2) >= gain
    return dict(n=len(t),gain=gain,semiamp_flux=np.hypot(mixed[::2],mixed[1::2]).tolist(),
                wild_night_p=(1+exceed)/(1+count),trials=count)


def period_sensitivity(d, p):
    results = []
    # In each block, compare an early and late fundamental clock at a common
    # overtone period, with one arbitrary phase per mode and no cross-block phase link.
    for trial in [.485254,.485766]:
        a = fit(d, [p[0], trial])
        results.append(dict(period=trial, chi2=a["chi2"], fundamental_amp=a["amp"][1]))
    return results


def period_bootstrap(d, p, rng, trials):
    from scipy.optimize import minimize_scalar
    groups = np.unique(np.floor(d[:,0]),return_inverse=True)[1]
    blocks = [np.where(groups==i)[0] for i in np.unique(groups)]
    recovered = []
    # Hold the other mode at the per-window value. Each bootstrap re-estimates
    # the fundamental inside a wider bracket than the measured shift.
    for _ in range(trials):
        draw = np.concatenate([blocks[i] for i in rng.integers(0,len(blocks),len(blocks))])
        sample = d[draw]
        def objective(v):
            return fit(sample,[p[0],v])["chi2"]
        g = np.linspace(.4846,.4864,17)
        k = np.argmin([objective(v) for v in g])
        res = minimize_scalar(objective,bounds=(g[max(0,k-1)],g[min(16,k+1)]),
                              method="bounded",options={"xatol":1.e-9})
        recovered.append(res.x)
    return np.percentile(recovered,[2.5,50,97.5]).tolist()


def main(cache, trials):
    rng = np.random.default_rng(8972)
    tests = {}
    for band,start,stop in [("I",5200,5750),("I",5750,6120),("I",6120,6485),
                             ("V",6485,6850),("I",6850,7215),("I",7215,7580)]:
        d = photo(4,NAME,band,cache)
        d = d[(d[:,0]>=start)&(d[:,0]<stop)]
        p = refine(d,PERIODS)
        tests[f"{band}_{start}_{stop}"] = dict(periods=p.tolist(),mixed=flux_check(d,p,rng,5000),
                                               fundamental_clock=period_sensitivity(d,p),
                                               fundamental_period_night_interval=period_bootstrap(d,p,rng,trials))
    return dict(eros=eros(cache),windows=tests,period_trials=trials)


if __name__ == "__main__":
    a = argparse.ArgumentParser(description=__doc__)
    a.add_argument("--cache",type=Path,default=Path("raw-cache"))
    a.add_argument("--out",type=Path,default=Path("independent-checks.json"))
    a.add_argument("--trials",type=int,default=400)
    args = a.parse_args()
    args.out.write_text(json.dumps(main(args.cache,args.trials),indent=2,allow_nan=False)+"\n")
