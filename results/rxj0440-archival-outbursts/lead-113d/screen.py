"""Replay the approximate 113-day search; see README for sample and limitations."""
import argparse
import csv
import hashlib
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
SKIP = {"GX1p4", "H1538-522", "SMCX-1"}


def emit(path, header, rows):
    with path.open("w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def glitch_screen(out):
    groups = defaultdict(list)
    with (HERE / "glitch-epochs.csv").open() as f:
        for row in csv.DictReader(f):
            groups[row["jname"]].append(float(row["mjd"]))
    freq = np.linspace(1/118, 1/108, 401)
    results = []
    rng = np.random.default_rng(113)
    for name, epochs in groups.items():
        t = np.sort(epochs)
        split = len(t) * 2 // 3
        phases = np.exp(2j*np.pi*(t[:, None]-t[0])*freq)
        train = phases[:split].mean(axis=0)
        i = abs(train).argmax()
        direction = train[i] / abs(train[i])
        held = (phases[split:,i] * direction.conjugate()).real.mean()
        null = np.cos(rng.uniform(0, 2*np.pi, (30000, len(t)-split))).mean(axis=1)
        p = (1+np.count_nonzero(null>=held))/(1+len(null))
        results.append((name, len(t), split, 1/freq[i], abs(train[i]), held, p, min(1.,p*len(groups))))
    emit(out / "glitch-scores.csv", ["jname","n","train_n","period_days","train_resultant",
         "heldout_projection","uniform_phase_p","bonferroni_p"], sorted(results, key=lambda r: -r[5]))


def get_bat(cache, item):
    path = cache / (item["name"]+".lc.txt")
    if not path.exists():
        with urlopen(item["url"], timeout=90) as inp:
            path.write_bytes(inp.read())
    if hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
        raise ValueError(f"Swift archive has updated {item['name']}. Supply the frozen input.")
    return item["name"], path


def bat_screen(cache, out):
    freq = np.linspace(1/118, 1/108, 181)
    results = []
    with (HERE/"bat-inputs.csv").open() as f:
        sources = list(csv.DictReader(f))
    cache.mkdir(exist_ok=True, parents=True)
    with ThreadPoolExecutor(6) as ex:
        paths = list(ex.map(lambda item:get_bat(cache,item), sources))
    rng = np.random.default_rng(20261114)
    for name,path in paths:
        if name in SKIP:
            continue
        a = np.loadtxt(path)
        a = a[(a[:,7]==0) & (a[:,2]>0) & (a[:,2]<.003)]
        t = a[:,0]
        train = t < 57400
        if len(a)<1200 or min(train.sum(), (~train).sum())<500:
            continue
        y = a[:,1] - np.median(a[train,1])
        limit = np.quantile(abs(y[train]), .98)
        y = np.clip(y, -limit, limit)
        w = np.minimum(a[:,2]**-2, np.quantile(a[train,2]**-2, .8))
        phase = np.exp(2j*np.pi*(t[:,None]-53000)*freq)
        def terms(mask):
            return w[mask]*(y[mask]-np.average(y[mask], weights=w[mask]))
        first = terms(train)
        z = (first @ phase[train]) / np.linalg.norm(first)
        i = abs(z).argmax()
        uu = z[i]/abs(z[i])
        last = terms(~train)
        contrib = last * phase[~train,i]*uu.conjugate()/np.linalg.norm(last)
        score = float(contrib.real.sum())
        # Each year retains its entire internal signal and gaps. Only its phase
        # relative to the frozen early ephemeris is randomized.
        years = a[~train,3].astype(int)
        blocks = np.array([contrib[years==year].sum() for year in np.unique(years)])
        angles = rng.uniform(0, 2*np.pi, (30000,len(blocks)))
        null = (np.exp(1j*angles) @ blocks).real
        p = (1+np.count_nonzero(null>=score))/(1+len(null))
        results.append([name,int(train.sum()),int((~train).sum()),1/freq[i],abs(z[i]),score,p])
    for row in results:
        row.append(min(1.,row[-1]*len(results)))
    emit(out/"bat-scores.csv", ["name","train_days","heldout_days","period_days","train_score",
         "heldout_projection","year_rotation_p","bonferroni_p"],
         sorted(results,key=lambda r:-min(r[4],r[5])))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    args.output.mkdir(exist_ok=True,parents=True)
    glitch_screen(args.output)
    bat_screen(args.cache,args.output)
