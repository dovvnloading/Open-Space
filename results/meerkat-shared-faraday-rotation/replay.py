"""Cross-pulsar tests for residual Faraday rotation. Run with NumPy."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from zipfile import ZipFile

import numpy as np

ROOT = Path(__file__).resolve().parent
EXAMPLE = "J1057-5226"


def pack(original, destination):
    """Transcribe the original measured rows and the source positions only."""
    if hashlib.md5(original.read_bytes()).hexdigest() != "f0eb4dffe8793a36069dbe5725bb3028":
        raise ValueError("Unexpected Zenodo tpa_data.zip checksum")
    with ZipFile(original) as archive, destination.open("w") as out:
        print("psr,ra_deg,dec_deg,mjd,rm,rm_err,iono,iono_err", file=out)
        for name in sorted(archive.namelist()):
            if not (name.startswith("tpa_data/") and name.endswith(".rmtable")):
                continue
            psr = Path(name).stem
            par = archive.read(name.replace(".rmtable", "_full.par")).decode()
            params = {s[0]: s[1] for line in par.splitlines() if len(s := line.split()) >= 2}
            if "RAJ" in params:
                def degrees(value):
                    fields = list(map(float, value.lstrip("+-").split(":")))
                    return (-1 if value.startswith("-") else 1) * sum(v / 60**i for i,v in enumerate(fields))
                ra = 15 * degrees(params["RAJ"])
                dec = degrees(params["DECJ"])
            else:
                lon, lat = np.radians([float(params["ELONG"]), float(params["ELAT"])])
                obliq = np.radians(23.4392911)
                x = np.cos(lon)*np.cos(lat)
                y = np.sin(lon)*np.cos(lat)*np.cos(obliq)-np.sin(lat)*np.sin(obliq)
                z = np.sin(lon)*np.cos(lat)*np.sin(obliq)+np.sin(lat)*np.cos(obliq)
                ra, dec = np.degrees(np.arctan2(y, x)) % 360, np.degrees(np.arcsin(z))
            for line in archive.read(name).decode().splitlines():
                if line.strip() and not line.startswith("#"):
                    print(",".join([psr, f"{ra:.5f}", f"{dec:.5f}", *line.split()]), file=out)


def read_rows():
    groups = defaultdict(list)
    with (ROOT / "measurements.csv").open() as f:
        for row in csv.reader(f):
            if row[0] != "psr":
                groups[row[0]].append(list(map(float,row[1:])))
    stars = {}
    for name, rows in groups.items():
        a = np.array(rows)
        a = a[np.all(np.isfinite(a),axis=1)]
        if len(a) < 18:
            continue
        a = a[np.argsort(a[:,2],kind="stable")]
        ra, dec = np.radians(a[0,:2])
        sky = np.array([np.cos(ra)*np.cos(dec), np.sin(ra)*np.cos(dec), np.sin(dec)])
        t, rm, rmerr, ion, ionerr = a[:,2:].T
        err = np.hypot(rmerr,ionerr)
        if (err < 2).sum() < 12:
            continue
        stars[name] = (t,rm-ion,err,sky,ion,rm)
    return stars


def peer_table(stars):
    sessions = defaultdict(list)
    for n,(t,y,e,sky,*_) in stars.items():
        med = np.median(y[e<2])
        for i in np.flatnonzero(e<1.8):
            sessions[int(t[i])].append((n,t[i],y[i]-med,sky))
    return sessions


def line(t,y,e):
    x = (t-np.mean(t))/365.25
    X = np.column_stack([np.ones(len(x)),x])
    w = e**-2
    cov = np.linalg.inv(X.T @ (w[:,None]*X))
    b = cov @ (X.T @ (w*y))
    res = y-X@b
    err = np.sqrt(cov[1,1]*np.sum(w*res**2)/(len(w)-2))
    return b[1],err,res


def experiment(stars, radius, hours, min_peers, split=False):
    days = peer_table(stars)
    output = []
    cases = {}
    tracks = []
    fold = {n:hashlib.sha256(n.encode()).digest()[0]%2 for n in stars}
    min_dot=np.cos(np.radians(radius))
    for n,(t,y,e,sky,*_) in stars.items():
        if split and fold[n] != 1:
            continue
        pred = np.full(len(t),np.nan)
        counts = np.zeros(len(t),int)
        for i,ti in enumerate(t):
            if e[i]>=2:
                continue
            peers = defaultdict(list)
            for neighbour,tp,r,p in days[int(ti)]:
                if neighbour!=n and (not split or fold[neighbour]==0) and abs(tp-ti)<hours/24 and sky@p>min_dot:
                    peers[neighbour].append(r)
            group=[v[0] if len(v)==1 else np.median(v) for v in peers.values()]
            if len(group)>=min_peers:
                pred[i] = np.median(group)
                counts[i] = len(group)
        ok = np.isfinite(pred)
        if ok.sum()<18 or np.ptp(t[ok])<800:
            continue
        tt,yy,ee,pp = t[ok],y[ok],e[ok],pred[ok]
        raw,rawerr,res = line(tt,yy,ee)
        adj,adjerr,adjres = line(tt,yy-pp,ee)
        full,fullerr,_ = line(t,y,e)
        output.append((n,len(t),len(tt),full,fullerr,raw,rawerr,adj,adjerr,
                       np.std(yy),np.std(yy-pp),np.std(res),np.std(adjres)))
        if n==EXAMPLE:
            cases[n]=list(zip(tt,yy,ee,pp,counts[ok]))
        # Detrend both series to test same-visit information beyond a common slope.
        _,_,p_res = line(tt,pp,ee)
        tracks.append((n,tt,res/np.std(res),p_res/np.std(res)))
    return output,cases,tracks


def summary(rows):
    a=np.array([r[1:] for r in rows])
    sig=abs(a[:,4])>3*a[:,5]
    corrected=abs(a[:,6])>3*a[:,7]
    return dict(pulsars=len(a), epochs=int(a[:,1].sum()),
                lower_scatter=int(sum(a[:,9]<a[:,8])),
                lower_detrended_scatter=int(sum(a[:,11]<a[:,10])),
                median_scatter_ratio=float(np.median(a[:,9]/a[:,8])),
                sig_before=int(sig.sum()), sig_after=int(corrected.sum()),
                lost=int((sig&~corrected).sum()), gained=int((~sig&corrected).sum()))


def coherence(tracks):
    # Random circular lags per target test phase alignment conditional on the observed
    # target and predicted waveforms. Shared sessions invalidate a global p-value.
    observed = np.mean([np.mean(y*p) for _,_,y,p in tracks])
    rng=np.random.default_rng(59860)
    draws=[]
    for _,t,y,p in tracks:
        draws.append([np.mean(y*np.roll(p,k)) for k in range(1,len(t))])
    null = np.mean([rng.choice(d,size=999) for d in draws], axis=0)
    # Flip entire 60-day bins together across all targets to preserve shared visits.
    blocks=defaultdict(float)
    for _,t,y,p in tracks:
        for ti,term in zip(t,y*p):
            blocks[int((ti-58000)//60)]+=term/(len(t)*len(tracks))
    block_draws=rng.choice([-1,1],size=(9999,len(blocks)))@np.array(list(blocks.values()))
    return dict(detrended_mean_cross_product=float(observed),
                signflip_60day_blocks=len(blocks),
                signflip_60day_quantile_99=float(np.quantile(block_draws,.99)),
                signflip_60day_exceedances=int(sum(block_draws>=observed)),
                signflip_60day_trials=len(block_draws),
                conditional_lag_range=[float(null.min()),float(null.max())],
                trials=len(null))


def independent_dates(stars):
    sessions = defaultdict(list)
    for name,(t,y,e,_,ion,raw) in stars.items():
        keep = e < 1.8
        if keep.sum() < 15:
            continue
        fold = hashlib.sha256(name.encode()).digest()[0] % 2
        for i in np.flatnonzero(keep):
            baseline = keep & (abs(t-t[i]) > 6) & (abs(t-t[i]) < 180)
            if baseline.sum() < 5:
                continue
            difference = [z[i]-np.median(z[baseline]) for z in (y,ion,raw)]
            sessions[int(t[i])].append((name,fold,difference))
    dates = []
    for day,visits in sorted(sessions.items()):
        halves = []
        for h in [0,1]:
            by_name=defaultdict(list)
            for name,fold,values in visits:
                if fold==h:
                    by_name[name].append(values)
            if len(by_name)<12:
                break
            vals=[np.median(v,axis=0) for v in by_name.values()]
            halves.extend([len(vals),*np.median(vals,axis=0)])
        if len(halves)==8:
            dates.append((day,*halves))
    a=np.array(dates)
    r=np.corrcoef(a[:,2],a[:,6])[0,1]
    return dates, dict(days=len(dates),pearson_r=float(r))


def write_csv(name,cols,rows):
    with (ROOT/name).open("w",newline="") as out:
        writer=csv.writer(out,lineterminator="\n")
        writer.writerow(cols.split(","))
        writer.writerows([[f"{v:.8g}" if isinstance(v,(np.floating,float)) else v for v in r] for r in rows])


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--extract",type=Path,help="Recreate measurements.csv from Zenodo tpa_data.zip")
    args=ap.parse_args()
    if args.extract:
        pack(args.extract,ROOT/"measurements.csv")
    stars=read_rows()
    settings=[("primary",30,2.16,5,False),("narrow",20,1.2,5,False),
              ("wide",40,3.,5,False),("ten_peers",30,2.16,10,False),
              ("disjoint_pulsars",30,2.16,5,True)]
    report={}
    for label,rad,h,k,split in settings:
        rows,cases,tracks=experiment(stars,rad,h,k,split)
        report[label]=summary(rows)
        if label in ["primary","disjoint_pulsars"]:
            report[label].update(coherence(tracks))
        if EXAMPLE in cases:
            match=next(r for r in rows if r[0]==EXAMPLE)
            report[label][EXAMPLE]=[float(v) for v in match[5:9]]
        if label=="primary":
            write_csv("slopes.csv","psr,n_full,n_matched,full_slope,full_se,matched_slope,matched_se,differential_slope,differential_se,rms_base,rms_diff,detrended_rms_base,detrended_rms_diff",rows)
            write_csv("J1057-5226.csv","mjd,corrected_rm,formal_error,other_pulsars_offset,n_other_pulsars",cases[EXAMPLE])
    dates,report["disjoint_dates"]=independent_dates(stars)
    write_csv("dates.csv","mjd_utc,n_A,resid_A,iono_A,measured_A,n_B,resid_B,iono_B,measured_B",dates)
    (ROOT/"checks.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report['primary']))
if __name__=="__main__":
    main()
