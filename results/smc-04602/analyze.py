"""Reproduce the 04602 target selection, folded eclipse timings and uncertainty tests."""
import argparse
import csv
import gzip
import json
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

BASE = Path(__file__).resolve().parent
PERIOD, EPOCH = 1.7295874, 2460001.4679
POINTINGS = ("sm0021k11605", "sm0057l21161")
SEED = 4618


def selection():
    groups = defaultdict(list)
    with (BASE.parent / "smc-06090/population/epoch_measurements.csv").open() as stream:
        for row in csv.DictReader(stream):
            groups[row["name"]].append(row)
    ranked = []
    complete = 0
    for name, rows in groups.items():
        rows.sort(key=lambda r: (int(r["survey"]), int(r["from_hjd"])))
        if len(rows) != 4 or any(r["quality"] != "usable" or float(r["se"]) > .015 for r in rows):
            continue
        complete += 1
        sep = [float(r["sep"]) for r in rows]
        se = [float(r["se"]) for r in rows]
        old, new = sep[1] - sep[0], sep[3] - sep[2]
        if old * new > 0:
            score = min(abs(old) / np.hypot(*se[:2]), abs(new) / np.hypot(*se[2:]))
            ranked.append((name, float(score), old, new))
    ranked.sort(key=lambda x: -x[1])
    with (BASE / "results/ranking.csv").open("w") as stream:
        w = csv.writer(stream)
        w.writerow(["id", "weaker_standardized_change", "OGLE_III_cycles", "OGLE_IV_cycles"])
        w.writerows(ranked)
    return {"total": len(groups), "four_usable": complete, "same_sign": len(ranked),
            "first_two": ranked[:2]}


def wrap(value):
    return (value + .5) % 1 - .5


def measure(data, start=None, shape="gaussian", period=PERIOD):
    # Phase relative to the modern primary, keeping like eclipses labeled across all epochs.
    t, mag, err = data.T
    phase = wrap((t - EPOCH) / period)
    harmonic = np.cos(4 * np.pi * phase)
    sigma = np.hypot(err, .012)
    baseline = np.percentile(mag, 35)
    lower = [baseline-.6, -.25, .001, -.23, .006, .001, .27, .006]
    upper = [baseline+.6, .25, 1.8, .23, .135, 1.8, .73, .135]

    def evaluate(x, jac=False):
        base, ell, a1, c1, w1, a2, c2, w2 = x
        z1, z2 = wrap(phase-c1)/w1, wrap(phase-c2)/w2
        if shape == "gaussian":
            g1, g2 = np.exp(-.5*z1**2), np.exp(-.5*z2**2)
            k1, k2 = z1*g1/w1, z2*g2/w2
        else:
            # Set triangle base half-width to 2.5 times the scale parameter.
            g1, g2 = np.maximum(0, 1-abs(z1)/2.5), np.maximum(0, 1-abs(z2)/2.5)
            k1 = np.sign(z1)*(g1 > 0)/(2.5*w1)
            k2 = np.sign(z2)*(g2 > 0)/(2.5*w2)
        if jac:
            return np.column_stack((np.ones(len(t)), harmonic, g1, a1*k1, a1*k1*z1,
                                    g2, a2*k2, a2*k2*z2))/sigma[:, None]
        return (base + ell*harmonic + a1*g1 + a2*g2 - mag)/sigma

    if start is None:
        def locate(center):
            grid = center + np.arange(-.18,.181,.012)
            weights = np.exp(-.5*(wrap(phase[:, None]-grid)/.025)**2) / sigma[:, None]**2
            return grid[np.argmax(weights.T @ np.clip(mag-baseline,-.5,1.5)/weights.sum(axis=0))]
        start = [baseline,.005,.42,locate(0),.024,.3,locate(.5),.024]
    res = least_squares(evaluate, np.clip(start, np.array(lower)+1e-6, np.array(upper)-1e-6),
                        jac=lambda x: evaluate(x,True), bounds=(lower,upper),
                        loss="soft_l1", f_scale=1.8, max_nfev=450, ftol=1e-7, xtol=1e-7)
    if not res.success:
        raise RuntimeError(res.message)
    return res.x


def ogle(survey):
    dat = np.loadtxt(BASE / f"data/{survey}-04602-I.dat", usecols=(0,1,2))
    dat = dat[np.isfinite(dat).all(axis=1) & (dat[:,2]>0) & (dat[:,2]<=.10)]
    if survey == "iii":
        dat[:,0] += 2450000
    return dat


def eros():
    path = BASE / "data/04602_xtime.tsv.gz"
    with gzip.open(path, "rt") as stream:
        lines = [line.rstrip("\n") for line in stream if line.strip() and not line.startswith("#")]
    rows = list(csv.DictReader([lines[0]]+lines[3:], delimiter="\t"))
    for src in POINTINGS:
        for band in ("B","R"):
            columns = ["HJD",band+"mag","e_"+band+"mag"]
            array = np.array([[float(r[k].strip() or "nan") for k in columns]
                              for r in rows if r["EROS"] == src])
            array[:,0] += 2450000
            keep = np.isfinite(array).all(axis=1) & (array[:,2]>0) & (array[:,2]<=.4)
            yield f"{src}-{band}", array[keep]


def samples():
    yield "OGLE-II", ogle("ii")
    survey3, survey4 = ogle("iii"), ogle("iv")
    for label, dat, bounds in [
        ("OGLE-III-early", survey3, (0,2453500)), ("OGLE-III-late",survey3,(2453500,2455200)),
        ("OGLE-IV-early",survey4,(0,2457000)), ("OGLE-IV-middle",survey4,(2457000,2459700)),
        ("OGLE-IV-late",survey4,(2459700,2500000))]:
        yield label, dat[(dat[:,0]>=bounds[0]) & (dat[:,0]<bounds[1])]
    # Historical EROS check uses disjoint dates within each pointing and band.
    # The HJD 2452000 round-number split was exploratory. Bands share exposures.
    for label, dat in eros():
        yield label+"-early", dat[dat[:,0]<2452000]
        yield label+"-late", dat[dat[:,0]>=2452000]


def analyze_group(task):
    label, data, draws, seed = task
    rng = np.random.default_rng(seed)
    centers = measure(data)
    phase = wrap((data[:,0]-EPOCH)/PERIOD)
    nights, indices = np.unique(np.floor(data[:,0]+.5), return_inverse=True)
    blocks = [np.flatnonzero(indices==i) for i in range(len(nights))]
    # Resample nights within fixed year-long strata to keep long-term coverage similar.
    seasons = np.floor((nights-2451700)/365.25).astype(int)
    strata = [np.flatnonzero(seasons==s) for s in np.unique(seasons)]
    sep = lambda x: float((x[6]-x[3])*PERIOD*24)
    bootstrap = []
    failed = 0
    while len(bootstrap) < draws:
        if failed > .03 * draws:
            raise ValueError(f"Unstable epoch {label}: {failed} refits failed.")
        chosen = np.concatenate([rng.choice(ix,len(ix),replace=True) for ix in strata])
        rows = np.concatenate([blocks[i] for i in chosen])
        try:
            bootstrap.append(sep(measure(data[rows],start=centers)))
        except RuntimeError:
            failed += 1
    triangle = measure(data,start=centers,shape="triangle")
    # Sensitivity to a one-part-per-million change over the full extrapolation.
    alternative = PERIOD*(1+1e-6)
    shifted = measure(data, period=alternative)
    return label, dict(n=len(data),nights=len(nights),median_hjd=float(np.median(data[:,0])),
        hjd_range=[float(data[:,0].min()),float(data[:,0].max())],
        primary_phase=float(centers[3]),secondary_phase=float(centers[6]),hours=sep(centers),
        interval95=np.percentile(bootstrap,[2.5,97.5]).tolist(),triangle_hours=sep(triangle),
        alternative_period_hours=float((shifted[6]-shifted[3])*alternative*24),
        core_nights=[len(np.unique(np.floor(data[abs(wrap(phase-x))<.035,0]+.5))) for x in centers[[3,6]]],
        failed_bootstrap_fits=failed,parameters=centers.tolist(),bootstrap=bootstrap)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--draws", type=int, default=600)
    parser.add_argument("--workers", type=int, default=4)
    args=parser.parse_args()
    (BASE/"results").mkdir(exist_ok=True)
    sel=selection()
    tasks=[(label,data,args.draws,SEED+i) for i,(label,data) in enumerate(samples())]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        measurements=dict(pool.map(analyze_group,tasks))
    pairs=[("OGLE-II","OGLE-IV-late"),("OGLE-III-early","OGLE-III-late"),
           ("OGLE-IV-early","OGLE-IV-late"),("OGLE-IV-middle","OGLE-IV-late")]
    pairs.extend((key+"-early",key+"-late") for key,_ in eros())
    comparisons={}
    for earlier,later in pairs:
        a,b=measurements[earlier],measurements[later]
        delta=np.subtract(b["bootstrap"],a["bootstrap"])
        comparisons[f"{earlier} to {later}"]={"hours":b["hours"]-a["hours"],
            "interval95":np.percentile(delta,[2.5,97.5]).tolist(),
            "triangle_hours":b["triangle_hours"]-a["triangle_hours"]}
    for measurement in measurements.values():
        del measurement["bootstrap"]
    out={"selection":sel,"period_d":PERIOD,"epoch_hjd":EPOCH,"draws_per_epoch":args.draws,
         "seed":SEED,"measurements":measurements,"comparisons":comparisons}
    (BASE/"results/timing.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps(comparisons,indent=2))


if __name__ == "__main__":
    main()
