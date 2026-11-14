"""Re-run the 1,800-target selection on OGLE bulk data laid out as described in README."""
import argparse
import csv
import gzip
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from astropy.timeseries import LombScargle
from scipy.ndimage import median_filter
from scipy.signal import find_peaks
from scipy.spatial import cKDTree

FIELDS = ["id","P","mag","n","scatter","p2","power","high","low","comment",
          "pref","a3","c3","c4","cv","av","neighbor","rank_score"]


def phase_residual(file, period):
    time, mag, err = np.loadtxt(file, unpack=True)
    good = (err > 0) & (err < max(0.03, 2.5*np.median(err))) & np.isfinite(mag)
    time, mag, err = time[good], mag[good], err[good]
    order = np.argsort((time / period) % 1)
    width = max(5, min(25, round(len(time)/150)))//2*2 + 1
    model = median_filter(mag[order], size=width, mode="wrap")[np.argsort(order)]
    residual = mag - model
    labels = np.floor((time-2455100)/365.25)
    for lab in np.unique(labels):
        residual[labels == lab] -= np.median(residual[labels == lab])
    scatter = max(1.483*np.median(abs(residual-np.median(residual))), 0.7*np.median(err), 0.0015)
    return time, err, np.clip(residual, -5*scatter, 8*scatter), scatter


def phase_bins(time, resid, err, f, size=24):
    b = ((time*f % 1)*size).astype(int)
    weights = 1/np.maximum(err, np.median(err)*0.85)**2
    return np.bincount(b, weights=weights*resid, minlength=size)/np.maximum(
        np.bincount(b, weights=weights, minlength=size), 1)


def match(time, resid, err, freq, profile):
    phase = time*freq % 1
    model = np.interp(phase,(np.arange(len(profile))+0.5)/len(profile),profile,period=1)
    weights = 1/np.maximum(err, np.median(err)*0.85)**2
    model -= np.average(model, weights=weights)
    resid = resid - np.average(resid, weights=weights)
    return (np.sum(weights*model*resid)/np.sum(weights*model**2), np.corrcoef(model,resid)[0,1])


def preliminary(name, period, mag, comment, path):
    time, err, resid, scatter = phase_residual(path, period)
    time -= 2450000
    if len(time) < 210:
        return []
    freq = np.arange(0.05, 5, 1/(3*np.ptp(time)))
    ls = LombScargle(time,resid,dy=np.maximum(err,np.median(err)*0.8))
    power = ls.power(freq, method="fast", assume_regular_frequency=True)
    peaks = find_peaks(power)[0]
    ordered = peaks[np.argsort(power[peaks])[-160:][::-1]]
    result, seen = [], []
    for j in ordered:
        f = freq[j]
        if any(abs(f-g)<0.004 or min(abs(f-g-i*1.00274) for i in (-2,-1,1,2))<0.002 for g in seen):
            continue
        if any(abs(abs(f-k/period)-i*1.00274)<0.002 for k in range(-12,13) for i in range(6)):
            continue
        if any(abs(f-i*0.002738)<0.0025 for i in range(1,16)):
            continue
        seen.append(f)
        b = ((time*f % 1)*32).astype(int)
        profile = np.array([np.mean(resid[b==i]) if np.any(b==i) else np.nan for i in range(32)])
        if np.isfinite(profile).sum()<27:
            continue
        smooth = (profile+np.roll(profile,1))/2
        result.append(dict(zip(FIELDS[:10],(name, period, mag, len(time), round(scatter,5),
                                            round(1/f,9), round(power[j],5),
                                            round(np.nanmax(smooth)/scatter,4),
                                            round(-np.nanmin(smooth)/scatter,4), comment))))
        if len(result)==4:
            break
    return result


def direction(ra, dec):
    def degrees(s):
        a,b,c=map(float,s.split(":"))
        return (-1 if s.startswith("-") else 1)*(abs(a)+b/60+c/3600)
    a,d=np.radians([15*degrees(ra),degrees(dec)])
    return [np.cos(a)*np.cos(d), np.sin(a)*np.cos(d), np.sin(d)]


def score(record):
    f = 1/record["pref"]
    ratio = record["P"]*f
    if (record["c3"]<0.06 or record["c4"]<0.10 or record["power"]<0.027 or
        min(abs(f-x) for x in range(7))<0.06 or
        min(abs(ratio-i/j) for i in range(1,43) for j in range(1,4))<0.002):
        return 0.
    return round(max(record["high"],0)**2/max(record["low"],0.05)*record["c3"]*np.sqrt(record["c4"]),4)


def run_target(args):
    name, period, magnitude, comment, older, neighbor, archives = args
    file = archives / "iv/phot/I" / (name+".dat")
    if not file.exists():
        return []
    rows = preliminary(name,period,magnitude,comment,file)
    if not rows:
        return []
    t,e,r,_ = phase_residual(file,period)
    t -= 2450000
    for row in rows:
        if row["power"]<=0.023 or row["high"]<=0.30:
            continue
        f = 1/row["p2"]
        grid = np.linspace(f-0.00010,f+0.00010,101)
        ls = LombScargle(t,r,dy=np.maximum(e,np.median(e)*0.85))
        f = grid[np.argmax(ls.power(grid))]
        profile = phase_bins(t,r,e,f)
        early = t<np.median(t)
        _,c4 = match(t[~early],r[~early],e[~early],f,phase_bins(t[early],r[early],e[early],f))
        a3=c3=av=cv=0.0
        if older and (file3:=archives/"iii/phot/I"/(older+".dat")).exists():
            t3,e3,r3,_ = phase_residual(file3,period)
            a3,c3 = match(t3,r3,e3,f,profile)
        if (v:=archives/"iv/phot/V"/(name+".dat")).exists():
            tv,ev,rv,_ = phase_residual(v,period)
            if len(tv)>12:
                av,cv = match(tv-2450000,rv,ev,f,profile)
        row.update(pref=round(1/f,9),a3=round(a3,3),c3=round(c3,3),c4=round(c4,3),
                   av=round(av,3),cv=round(cv,3),neighbor=neighbor)
        row["rank_score"]=score(row)
    return rows


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("archives",type=Path)
    ap.add_argument("output",type=Path)
    ap.add_argument("--workers",type=int,default=2)
    args=ap.parse_args()
    def lines(n):
        return (args.archives/n).read_text().splitlines()
    params={v[0]:(float(v[3]),float(v[1])) for line in lines("iv_ecl.dat")
            if (v:=line.split())[3]!="-"}
    notes={v[0]:" ".join(v[1:]) for line in lines("iv_remarks.txt") if (v:=line.split())}
    targets=sorted(params,key=lambda n:(n not in notes,params[n][1]))[:1800]
    ids4={v[0]:direction(*v[2:4]) for line in lines("iv_ident.dat") if (v:=line.split())}
    ids3={v[0]:direction(*v[4:6]) for line in lines("iii_ident.dat") if (v:=line.split())}
    keys3,keys4=list(ids3),list(ids4)
    tree3,tree4=cKDTree(list(ids3.values())),cKDTree(list(ids4.values()))
    tasks=[]
    for name in targets:
        d,i=tree3.query(ids4[name])
        close,indices=tree4.query(ids4[name],k=3)
        neighbor=" ".join(f"{keys4[i]}@{d*206265:.2f}@{params.get(keys4[i],(0,))[0]}" for d,i in zip(close[1:],indices[1:]))
        tasks.append((name,*params[name],notes.get(name,""),keys3[i] if d*206265<1.2 else "",neighbor,args.archives))
    with Pool(args.workers) as pool, gzip.open(args.output,"wt",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=FIELDS,lineterminator="\n");writer.writeheader()
        for rows in pool.imap(run_target,tasks,chunksize=4):
            writer.writerows(rows)


if __name__ == "__main__":
    main()
