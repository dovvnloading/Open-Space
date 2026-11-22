"""Measure an instrumental Ca II H&K proxy in individual overlapping HARPS spectra."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import io
from pathlib import Path
import re
import time
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.io import fits
from scipy.stats import f as f_dist, spearmanr

from holdout import key
from runtime import require_survey_workflow


MAX_STARS = 3
MAX_NIGHTS = 24
MAX_SECONDS = 360
MAX_FILE = 16*1024**2
MAX_TOTAL = 400*1024**2
C_MS = 299792458.
TIME = re.compile(r"\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}(?:\.\d+)?")
ADP = re.compile(r"ADP\.\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d+")
BANDS = {"h":(3968.470,1.09), "k":(3933.664,1.09),
         "v":(3901.070,10.), "r":(4001.070,10.)}


def read(path):
    return pd.read_csv(path,dtype={"star":str},keep_default_na=False,low_memory=False)


def timestamp(value):
    token = TIME.search(str(value))
    return pd.to_datetime(token[0],utc=True,errors="coerce") if token else pd.NaT


def select_peaks(out):
    s = read(out/"signals.csv")
    s = s[s["flags"].eq("missing_sindex") & s.rv_status.eq("inconclusive") &
          s.fap.le(.01) & s.delta_bic.ge(10)].copy()
    v = read(out/"vetting.csv")
    cols = ["star","signal_index","identity","aliases","external_status","prior_signal",
            "alias_matches","rotation_check","binary_check","minimum_mass_mjup"]
    s = s.merge(v[cols],on=["star","signal_index"],validate="one_to_one")
    m = pd.to_numeric(s.minimum_mass_mjup,errors="coerce")
    elig = (s.external_status.isin(["inconclusive","requires literature review"]) &
            s.prior_signal.eq("") & s.alias_matches.eq("") & (m.isna() | m.lt(80)) &
            s.rotation_check.isin(["clear of catalog periods","unavailable"]) &
            s.binary_check.isin(["no catalog flag","incomplete"]))
    s = s[elig].copy()
    s["priority"] = s.star.map(key).eq("HD17970")
    return s.sort_values(["priority","delta_bic","star","signal_index"],ascending=[False,False,True,True])


def match_exposures(primary,esp):
    p = primary.sort_values("bjd").reset_index(drop=True)
    e = esp.sort_values("bjd").reset_index(drop=True)
    left,right = p.bjd.to_numpy(),e.bjd.to_numpy()
    tol = 2/86400
    lo,hi = np.searchsorted(right,left-tol),np.searchsorted(right,left+tol,"right")
    rl,rh = np.searchsorted(left,right-tol),np.searchsorted(left,right+tol,"right")
    rows = []
    for i in np.flatnonzero(hi-lo==1):
        j = lo[i]
        if rh[j]-rl[j]!=1 or e.era.iloc[j]!=p.era.iloc[i]:
            continue
        raw, adp = p.epoch_id.iloc[i], str(e.epoch_id.iloc[j]).removesuffix(".fits")
        if pd.isna(timestamp(raw)) or not ADP.fullmatch(adp):
            continue
        rows.append(dict(bjd=left[i],eso_bjd=right[j],era=p.era.iloc[i],raw_id=raw,adp=adp,
                         rv_ms=p.rv.iloc[i],rv_error_ms=p.err.iloc[i],drs_absolute_rv_ms=e.rv.iloc[j],
                         night=int(np.floor(left[i]))))
    if not rows:
        return pd.DataFrame()
    nights = (pd.DataFrame(rows).sort_values(["night","rv_error_ms","bjd"])
              .drop_duplicates("night").sort_values("bjd"))
    groups = [frame for _,frame in nights.groupby("era",sort=True)]
    if len(groups)==2 and min(map(len,groups))>=10:
        split = min(12,min(map(len,groups)))
        small = min(range(2),key=lambda i:len(groups[i]))
        quotas = [MAX_NIGHTS-split]*2
        quotas[small] = split
    else:
        groups = [max(groups,key=len)]
        quotas = [MAX_NIGHTS]
    return pd.concat([g.iloc[np.linspace(0,len(g)-1,min(len(g),n)).round().astype(int)]
                      for g,n in zip(groups,quotas)]).sort_values("bjd")


def fetch(session,adp,deadline,budget,rec):
    url = "https://dataportal.eso.org/dataPortal/file/"+quote(adp,safe=":")
    rec.update(url=url,retrieved_utc=datetime.now(timezone.utc).isoformat(),bytes=0)
    with session.get(url,stream=True,timeout=(6,12)) as resp:
        rec["http_status"] = resp.status_code
        resp.raise_for_status()
        chunks = []
        for block in resp.iter_content(65536):
            rec["bytes"] += len(block)
            budget[0] += len(block)
            if rec["bytes"]>MAX_FILE or budget[0]>MAX_TOTAL or time.monotonic()>deadline:
                raise OSError("download budget exceeded")
            chunks.append(block)
    data = b"".join(chunks)
    rec["sha256"] = hashlib.sha256(data).hexdigest()
    return data


def measure(data,obs,rec):
    with fits.open(io.BytesIO(data),memmap=False) as hdus:
        header = hdus[0].header
        ext = hdus[1]
        if str(header.get("INSTRUME","")).strip().upper()!="HARPS":
            raise ValueError("instrument is not southern HARPS")
        observed = timestamp(header.get("DATE-OBS",""))
        if pd.isna(observed):
            observed = timestamp(header.get("ORIGFILE",""))
        diff = (observed-timestamp(obs.raw_id)).total_seconds()
        rec.update(utc_difference_s=diff,header_object=header.get("OBJECT",""),
                   drift_quality=header.get("ESO DRS DRIFT QC",""),
                   pipeline=header.get("ESO DRS VERSION",header.get("PROCSOFT","")))
        if not np.isfinite(diff) or abs(diff)>1.0:
            raise ValueError("FITS exposure start fails UTC cross-check")
        if str(rec["drift_quality"]).strip().upper() in {"FAILED","FAIL","FALSE","F"}:
            raise ValueError("pipeline drift quality failure")
        if "ESO DRS BJD" in header and abs(float(header["ESO DRS BJD"])-obs.bjd)>2/86400:
            raise ValueError("FITS and RVBank BJD differ by >2 seconds")
        frame = str(ext.header.get("SPECSYS",header.get("SPECSYS",""))).upper()
        if frame and "BARY" not in frame:
            raise ValueError("spectral frame conflicts with HARPS barycentric s1d")
        if len(ext.data)!=1 or not {"WAVE","FLUX"}.issubset(ext.columns.names):
            raise ValueError("unexpected ESO 1D binary table")
        cols = ext.columns
        wave = (np.asarray(ext.data["WAVE"][0],float)*u.Unit(cols["WAVE"].unit)).to_value(u.AA)
        flux = np.asarray(ext.data["FLUX"][0],float)
        # The classic ESO DRS s1d product is in barycentric air wavelengths.
        # Where a product explicitly identifies vacuum wavelengths, require a
        # different reduction rather than mixing wavelengths between epochs.
        medium = " ".join(str(ext.header.get(k,header.get(k,""))) for k in ("AIRORVAC","WAVEMED"))
        if "VAC" in medium.upper():
            raise ValueError("explicit vacuum wavelength scale")
        beta = obs.drs_absolute_rv_ms/C_MS
        if not np.isfinite(beta) or abs(beta)>.005:
            raise ValueError("absolute stellar RV unavailable or out of range")
        wave = wave/np.sqrt((1+beta)/(1-beta))
        rec["wave_frame"] = frame or "documented HARPS s1d barycentric"
        for band,(center,width) in BANDS.items():
            mask = abs(wave-center)<=width
            w,y = wave[mask],flux[mask]
            if len(w)<20 or np.any(~np.isfinite(y)) or np.any(np.diff(w)<=0):
                raise ValueError("missing or nonfinite "+band+" samples")
            if ((w[-1]-w[0])/(2*width)<.98 or np.diff(w).max()>.075 or
                    np.mean(y>0)<.98):
                raise ValueError("incomplete, gapped or nonpositive "+band+" flux")
            wt = np.maximum(0.,1-abs(w-center)/width) if band in ("h","k") else np.ones(len(w))
            rec["f_"+band] = np.trapezoid(y*wt,w)/np.trapezoid(wt,w)
            if band in ("r","v"):
                diff = np.diff(y)
                noise = 1.4826*np.median(abs(diff-np.median(diff)))/np.sqrt(2.)
                rec["sn_proxy_"+band] = np.median(y)/noise if noise>0 else np.nan
                if not np.isfinite(rec["sn_proxy_"+band]) or rec["sn_proxy_"+band]<8:
                    raise ValueError("low blue-continuum quality: "+band)
        rec["hk_ratio"] = (rec["f_h"]+rec["f_k"])/(rec["f_r"]+rec["f_v"])
        rec["v_to_r"] = rec["f_v"]/rec["f_r"]
        if not np.isfinite(rec["hk_ratio"]) or rec["hk_ratio"]<=0:
            raise ValueError("invalid H&K flux ratio")


def compare(usable,sig,era):
    row = dict(star=sig.star,signal_index=sig.signal_index,period_d=sig.period_d,era=era,
               n_nights=len(usable),interpretation="insufficient spectral coverage")
    if len(usable)<10:
        return row
    t = usable.bjd.to_numpy()
    span = np.ptp(t)
    phase = ((t-sig.tperi_bjd)/sig.period_d)%1
    phase_gap = np.diff(np.sort(np.r_[phase,phase.min()+1])).max()
    row.update(span_d=span,phase_gap=phase_gap)
    if span<sig.period_d or phase_gap>.4:
        return row
    a = np.column_stack([np.ones(len(t)),(t-t.mean())/365.25])
    def resid(x,design=a):
        return x-design@np.linalg.lstsq(design,x,rcond=None)[0]
    s = resid(usable.hk_ratio.to_numpy())
    rv = resid(usable.rv_ms.to_numpy())
    rho, pcorr = spearmanr(s,rv)
    row.update(rv_spearman_rho=rho,rv_spearman_p=pcorr,
               colour_spearman_rho=spearmanr(s,resid(usable.v_to_r.to_numpy())).statistic)
    vals = []
    for harmonic,label in [(1,"orbital"),(.5,"half_frequency"),(2,"double_frequency")]:
        angle = 2*np.pi*(t-t.min())*harmonic/sig.period_d
        b = np.column_stack([a,np.sin(angle),np.cos(angle)])
        r = resid(usable.hk_ratio.to_numpy(),b)
        if np.linalg.cond(b)>1e7 or np.sum(r*r)<=0:
            row["interpretation"]="degenerate diagnostic"
            return row
        dof = len(t)-b.shape[1]
        prob = f_dist.sf(max(0.,(s@s-r@r)/(r@r)*dof/2),2,dof)
        row[label+"_p"] = prob
        vals.append(prob)
    if not np.isfinite([rho,pcorr,*vals]).all():
        row["interpretation"]="degenerate diagnostic"
    elif (abs(rho)>=.4 and pcorr<.01/7) or min(vals)<.01/21:
        row["interpretation"]="association with instrumental H&K proxy"
    else:
        row["interpretation"]="no association resolved in this subset"
    return row


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("work","out","eso-epochs"):
        p.add_argument("--"+name,type=Path,required=True)
    args = p.parse_args()
    require_survey_workflow()
    peaks = select_peaks(args.out)
    queue = peaks.star.drop_duplicates().tolist()
    pri = read(args.work/"epochs.csv.gz")
    eso = read(args.eso_epochs)
    names = defaultdict(set)
    for name in eso.star.unique():
        names[key(name)].add(name)
    primary = pri.groupby("star",sort=False)
    archive = eso.groupby("star",sort=False)
    deadline, budget = time.monotonic()+MAX_SECONDS,[0]
    epoch_rows, trials = [],[]
    outages = 0
    with requests.Session() as session:
        session.headers["User-Agent"]="Open-Space HARPS spectrum activity check"
        for i,star in enumerate(queue):
            fits_for_star = peaks[peaks.star.eq(star)]
            first = fits_for_star.iloc[0]
            issue = ""
            if i>=MAX_STARS:
                issue="outside three-target budget"
            choices = {key(v) for v in [star,first.identity,*first.aliases.split("|")] if str(v).strip()}
            ids = set().union(*(names.get(v,set()) for v in choices))
            if len(ids)!=1:
                issue = issue or "no unique component-preserving ESO identifier"
            if issue:
                for sig in fits_for_star.itertuples():
                    trials.append(dict(star=star,signal_index=sig.signal_index,period_d=sig.period_d,
                                       interpretation=issue,n_nights=0))
                continue
            name = next(iter(ids))
            selected = match_exposures(primary.get_group(star),archive.get_group(name))
            observed=[]
            for obs in selected.itertuples(index=False):
                row = dict(star=star,eso_id=name,**obs._asdict(),status="pending")
                epoch_rows.append(row)
                if outages>=3 or time.monotonic()>deadline or budget[0]>=MAX_TOTAL:
                    row["status"]="network/time/byte budget exhausted"
                    continue
                try:
                    data = fetch(session,obs.adp,deadline,budget,row)
                    outages=0
                except (OSError,requests.RequestException) as exc:
                    outages+=1
                    row.update(status="download failed",detail=str(exc)[:250])
                    continue
                try:
                    measure(data,obs,row)
                except (ValueError,KeyError,IndexError,OSError,TypeError) as exc:
                    row.update(status="spectrum unusable",detail=str(exc)[:250])
                    continue
                row["status"]="measured"
                observed.append(row)
            measured = pd.DataFrame(observed)
            for sig in fits_for_star.itertuples():
                if measured.empty:
                    trials.append(dict(star=star,signal_index=sig.signal_index,period_d=sig.period_d,
                                       interpretation="no usable measured nights",n_nights=0,
                                       n_selected=len(selected),eso_id=name))
                else:
                    for era,group in measured.groupby("era",sort=True):
                        trials.append(dict(**compare(group,sig,era),n_selected=len(selected),eso_id=name))
    measurements = pd.DataFrame(epoch_rows)
    fields = ["star","eso_id","raw_id","adp","bjd","era","rv_ms","hk_ratio","status","url","sha256"]
    measurements.reindex(columns=fields+[c for c in measurements if c not in fields]).to_csv(
        args.out/"spectral-activity-epochs.csv",index=False)
    result = pd.DataFrame(trials)
    fields = ["star","signal_index","period_d","era","n_nights","interpretation"]
    result.reindex(columns=fields+[c for c in result if c not in fields]).to_csv(
        args.out/"spectral-activity-targets.csv",index=False)
    print(f"H&K follow-up: {len(queue)} eligible target names, "
          f"{sum(r['status']=='measured' for r in epoch_rows)} usable spectra",flush=True)


if __name__=="__main__":
    main()
