#!/usr/bin/env python3
"""Frozen-finder injection, shuffled and local BACK proxy checks (NOT empirical archive FPR).

Also tests whether range requests can fetch the checksum-pinned 4XMM FITS
in parallel. No change to the in-sample control thresholds.
"""
import datetime as dt
import hashlib
import json
from pathlib import Path
import time

import numpy as np
import requests
from control_validation import BASE, OUT, TARGETS, find_eruptions
from replay_controls import from_csv

RNG = np.random.default_rng(14002026)

def found(a):
    return bool(a.get("n_peaks",0)>=2 and a.get("range_match",False))

def copy(curve, rate, error=None):
    d=dict(curve)
    d["rate"]=rate.copy()
    if error is not None: d["error"]=error.copy()
    return d

def cat_probe():
    urls=[
        "https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/4XMM_DR14cat_v1.0.fits.gz",
        "https://xmm-ssc.irap.omp.eu/Catalogue/4XMM-DR14/4XMM_DR14cat_v1.0.fits.gz",
        "https://xmmssc.aip.de/Catalogue/4XMM-DR14/4XMM_DR14cat_v1.0.fits.gz"
    ]
    ss=requests.Session()
    rec=[]
    for u in urls:
        r={"url":u}
        t=time.monotonic()
        try:
            with ss.get(u,headers={"Range":"bytes=0-1023"},
                        stream=True,timeout=(12,15)) as response:
                block=next(response.iter_content(1024),b"")
                r.update(http=response.status_code,content_range=response.headers.get("Content-Range"),
                         content_length=response.headers.get("Content-Length"),
                         accept_ranges=response.headers.get("Accept-Ranges"),
                         content_type=response.headers.get("Content-Type"),
                         first_block_sha256=hashlib.sha256(block).hexdigest(),bytes_read=len(block))
        except Exception as exc:
            r["error"]=repr(exc)[:300]
        r["elapsed_s"]=round(time.monotonic()-t,2)
        rec.append(r)
    (OUT/"archive_http_range_probe.json").write_text(json.dumps(rec,indent=2)+"\n")
    return rec

def main():
    cfg=json.loads((BASE/"control_batch.json").read_text())
    assert cfg.get("stage")=="archive_diagnostics"
    scans=[]; injections=[]; losses=[]
    for target,info in TARGETS.items():
        name=target.replace(" ","_").replace("+","p")
        path=OUT/f"control_curve_{name}_{info['obsid']}.csv"
        assert path.is_file(), f"missing saved control curve: {path.name}"
        c=from_csv(path)
        valid=c["good"] & np.isfinite(c["rate"]) & np.isfinite(c["error"]) & np.isfinite(c["background"])
        good=np.flatnonzero(valid)
        y=c["rate"].copy(); e=c["error"].copy()
        if len(good)<30:
            losses.append({"control":target,"reason":"few good bins"});continue
        standard=find_eruptions(c,50814.0,(0.5,24))
        shuffle_flags=0;shuffle_peakcount=0
        for trial in range(35):
            shuffled=y.copy()
            shuffled[good]=RNG.permutation(y[good])
            result=find_eruptions(copy(c,shuffled),50814.0,(0.5,24))
            shuffle_flags+=int(found(result))
            shuffle_peakcount+=result.get("n_peaks",0)
        proxy=copy(c,c["background"])
        proxy_result=find_eruptions(proxy,50814.0,(0.5,24))
        scans.append({"control":target,"n_real_peaks":standard.get("n_peaks"),
            "time_shuffles_with_recurrence_flag":shuffle_flags,
            "time_shuffle_trials":35,"time_shuffled_total_detected_peaks":shuffle_peakcount,
            "local_scaled_background_proxy_recurrence":found(proxy_result),
            "local_scaled_background_proxy_n_peaks":proxy_result.get("n_peaks",0),
            "caveat":"Full-data permutations include actual QPE bins; the shuffled positives are NOT independent false positives. BACK proxy is scaled source background, not a separate region."})
        # Real-bin injection into an explicitly quiet-resampled baseline.
        quiet=np.quantile(y[good],.40)
        quiescent=y[good][y[good]<=np.quantile(y[good],.60)]
        if len(quiescent)<12:
            losses.append({"control":target,"reason":"too few quiet bins"});continue
        valid_centers=[int(i) for i in good if i>=6 and i+6<len(y)
                       and bool(np.all(valid[i-4:i+5]))]
        errtyp=float(np.nanmedian(e[good]))
        valid_centers=np.array(valid_centers,dtype=int)
        for width in (450.,850.,1300.):
            for height_sigma in (4,8,16):
                for repeat in range(6):
                    if len(valid_centers)<2:continue
                    ip=int(RNG.choice(valid_centers))
                    possibles=valid_centers[
                        (np.abs(c["t_s"][valid_centers]-c["t_s"][ip])>=2*3600) &
                        (np.abs(c["t_s"][valid_centers]-c["t_s"][ip])<=22*3600)]
                    if not len(possibles):
                        losses.append({"control":target,"reason":"no 2-22h exposure pair"})
                        continue
                    jp=int(RNG.choice(possibles))
                    baseline=y.copy()
                    baseline[good]=RNG.choice(quiescent,size=len(good),replace=True)
                    amp=max((height_sigma-1)*max(float(quiet),.005),height_sigma*errtyp)
                    injected=amp*(np.exp(-.5*((c["t_s"]-c["t_s"][ip])/width)**2)+
                                  np.exp(-.5*((c["t_s"]-c["t_s"][jp])/width)**2))
                    shifted=baseline+injected
                    ee=np.sqrt(e**2 + injected/np.maximum(c["exposure"],1))
                    detected=find_eruptions(copy(c,shifted,ee),50814.,(0.5,24.))
                    predicted=[c["t_s"][ip],c["t_s"][jp]]
                    actual=[(float(v["mjd_tt"])-50814.)*86400 for v in detected.get("peaks",[])]
                    one_to_one=[]
                    used=set()
                    for t in predicted:
                        options=sorted((abs(d-t),k) for k,d in enumerate(actual) if k not in used)
                        if options and options[0][0]<=900:
                            one_to_one.append(options[0][0]);used.add(options[0][1])
                    injections.append({"control":target,"width_s":width,"height_sigma_input":height_sigma,
                          "repeat":repeat,"n_injected":2,"n_time_matched":len(one_to_one),
                          "recurrent_train":found(detected),
                          "both_recovered":len(one_to_one)==2 and found(detected)})
    summary={"as_of_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "method":"frozen control_validation.find_eruptions; 0.2-2keV 300s saved NET PN curves; fixed seed 14002026",
        "controls":scans,"injection_grid":injections,
        "n_shuffled_trials":sum(r["time_shuffle_trials"] for r in scans),
        "n_shuffled_recurrent_flags":sum(r["time_shuffles_with_recurrence_flag"] for r in scans),
        "n_injection_trials":len(injections),
        "n_two_injection_recovery":sum(x["both_recovered"] for x in injections),
        "coverage_losses":losses,
        "false_positive_rate":None,"archive_screen_permitted":False,
        "limitations":["In-sample six positive controls only; no independent stellar, blank, AGN or archive-wide negative population",
          "Rate permutations contain source eruptions, destroying correlated flare morphology; not a calibrated FDR",
          "Local BACK1-3 is a scaled proxy, not separately extracted region or particle flare GTI",
          "Injected Gaussian rates are bright positive-control-style surrogates, not Poisson event-list injections",
          "Injection excludes non-detections/host incompleteness/pile-up; cannot yield survey occurrence rate"]}
    OUT.joinpath("archive_statistical_checks.json").write_text(json.dumps(summary,indent=2,allow_nan=False)+"\n")
    lines=["# DR14 frozen-finder diagnostic checks — NOT calibrated survey FPR","",
        f"Time-shuffles including control eruptions: **{summary['n_shuffled_recurrent_flags']}/{summary['n_shuffled_trials']}** recurrent-train flags (not independent nulls).",
        f"Two-pulse injections into resampled quiet real PN bins: **{summary['n_two_injection_recovery']}/{summary['n_injection_trials']}** paired recoveries (synthetic in-sample, not completeness).",
        "Independent archive FPR: **UNESTIMATED**. Background-region controls: **NOT RUN**; only local scaled BACK proxy tested.",
        "Archive gate: **CLOSED**. No observed new-QPE sample was screened.","",
        "## Limitations"]+["- "+x for x in summary["limitations"]]
    OUT.joinpath("archive_statistical_checks.md").write_text("\n".join(lines)+"\n")
    probes=cat_probe()
    print(json.dumps({"n_shuffled_flags":summary["n_shuffled_recurrent_flags"],
        "n_shuffle_trials":summary["n_shuffled_trials"],
        "n_injections_recovered":summary["n_two_injection_recovery"],
        "n_injection_trials":summary["n_injection_trials"],
        "range_probes":[{k:x.get(k) for k in ("url","http","content_range","content_length","error")}
                        for x in probes]},indent=2),flush=True)

if __name__=="__main__":
    main()
