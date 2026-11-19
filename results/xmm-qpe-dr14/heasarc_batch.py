#!/usr/bin/env python3
"""Frozen 4XMM-DR14 host cohort: bounded HEASARC PN PPS screen.

Only PN merged, position-verified, HDUCLAS2=NET RATE1+RATE2+RATE3
(0.2-2 keV) is considered screened. Download failures are never negatives.
FBKTSR is a source-masked, independent-field particle-flare diagnostic,
NOT a validated off-source 0.2-2 keV background spectrum or pileup test.
"""
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
import re
import time
from pathlib import Path

import numpy as np
import requests
from astropy.io import fits

import control_validation as cv
from heasarc_screen import HOST, listing

BASE = Path(__file__).resolve().parent
OUT = BASE / "outputs"
MANIFEST = OUT / "archive_provisional_nuclei_sdss_gaia.csv"
FREEZE = OUT / "archive_final_host_cohort.json"
SESSION = requests.Session()
SESSION.headers["User-Agent"] = "Open-Space-XMM-DR14-HEASARC-PN-screen/1.0"
MAX_FILE = 70 * 1024**2
CUTOFF = 39 * 60

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def save(name, obj):
    tmp = OUT / (name + ".tmp")
    tmp.write_text(json.dumps(cv.safe(obj), sort_keys=True, indent=2, allow_nan=False) + "\n")
    tmp.replace(OUT / name)

def manifest():
    fr = json.loads(FREEZE.read_text())
    assert cv.checksum(MANIFEST) == fr["final_csv_sha256"], "Frozen host cohort checksum changed"
    rows = list(csv.DictReader(MANIFEST.open(newline="")))
    assert len(rows) == 967 and len({r["detid"] for r in rows}) == 967
    obs = {}
    for r in rows:
        assert re.fullmatch(r"\d{10}", r["obs_id"])
        obs.setdefault(r["obs_id"], []).append(r)
    assert len(obs) == fr["post_gaia_provisional_unique_obsids"] == 670
    return obs, fr

def get_file(obsid, name, receipts, attempts=3):
    # Literal names only from directory listing: no redirects to another host.
    assert re.fullmatch(r"P" + obsid + r"PN[A-Z0-9_.+-]{8,90}\.FTZ", name, re.I), name
    url = HOST + obsid + "/PPS/" + name
    for attempt in range(1, attempts + 1):
        t = time.monotonic()
        rec = {"obs_id": obsid, "filename": name, "url": url,
               "attempt": attempt, "time_utc": utc(), "limit_bytes": MAX_FILE}
        try:
            with SESSION.get(url, timeout=(12, 45), stream=True) as response:
                rec["http_status"] = response.status_code
                rec["content_length"] = response.headers.get("content-length")
                rec["final_url"] = response.url
                if not response.url.startswith(HOST + obsid + "/PPS/"):
                    raise ValueError("Refused redirect outside NASA HEASARC observation PPS")
                response.raise_for_status()
                if int(response.headers.get("Content-Length") or "0") > MAX_FILE:
                    raise ValueError("File exceeds 70MiB compressed safety ceiling")
                body = bytearray()
                sha = hashlib.sha256()
                for block in response.iter_content(2**18):
                    if not block: continue
                    body.extend(block)
                    sha.update(block)
                    if len(body) > MAX_FILE:
                        raise ValueError("File streamed past 70MiB compressed limit")
                rec["bytes"] = len(body)
                rec["compressed_sha256"] = sha.hexdigest()
                if len(body) < 100:
                    raise ValueError("Truncated compressed product")
                raw = gzip.decompress(body) if name.upper().endswith(".FTZ") else bytes(body)
                rec["uncompressed_bytes"] = len(raw)
                rec["fits_sha256"] = hashlib.sha256(raw).hexdigest()
                if len(raw) < 2880 or not raw.startswith(b"SIMPLE"):
                    raise ValueError("Not FITS")
                rec["status"] = "downloaded_and_gzip_verified"
                return raw
        except Exception as exc:
            rec["status"] = "failed"
            rec["error"] = repr(exc)[:300]
            if attempt < attempts: time.sleep(attempt)
        finally:
            rec["elapsed_s"] = round(time.monotonic() - t, 2)
            receipts.append(rec)
    return None

def net_merged(data):
    with fits.open(io.BytesIO(data), memmap=False) as hd:
        rates = next((x for x in hd[1:] if x.data is not None and
             hasattr(x.data, "names") and all(f"RATE{i}" in x.data.names for i in (1,2,3))), None)
        if rates is None:
            raise ValueError("No 0.2-2 keV merged RATE1-3; broad-band RATE is insufficient")
        if str(rates.header.get("HDUCLAS2", "")).strip().upper() != "NET":
            raise ValueError("Source light curve is not certified NET")
        if str(rates.header.get("TIMESYS", "TT")).upper() != "TT":
            raise ValueError("Unexpected time system; require TT")
        h = rates.header
        if "SRC_RA" not in h or "SRC_DEC" not in h:
            raise ValueError("No explicit source coordinates in PN FITS")
        position = (float(h["SRC_RA"]), float(h["SRC_DEC"]))
        assert 0 <= position[0] < 360 and -90 <= position[1] <= 90
        t, y, e, bg, frac, hard, td = cv.combine_energies(rates)
        ref = cv.mjd_reference(h)
        curve = cv.rebin(t, y, e, bg, frac, hard, td)
        hdr = {"hduclas2": "NET", "mjdref_tt": ref,
               "source_position_deg": list(position),
               "columns": list(rates.data.names), "n_native_bins": len(t),
               "source_gti_hdu_present": "SRC_GTIS" in hd, "background_gti_hdu_present": "BKG_GTIS" in hd}
    return position, curve, ref, hdr

def flare_background(raw, name):
    try:
        with fits.open(io.BytesIO(raw), memmap=False) as hd:
            h = next((x for x in hd[1:] if x.data is not None and
                      hasattr(x.data, "names") and
                      {"TIME","RATE","ERROR","FRACEXP"}.issubset(x.data.names)), None)
            if h is None: return {"status":"no_time_rate_columns"}
            arr=h.data
            t=np.asarray(arr["TIME"],float)
            y=np.asarray(arr["RATE"],float)
            e=np.asarray(arr["ERROR"],float)
            f=np.asarray(arr["FRACEXP"],float)
            good=np.isfinite(t)&np.isfinite(y)&np.isfinite(e)&(f>=.5)&(e>=0)
            if np.count_nonzero(good)<15:
                return {"status":"insufficient_offsource_bins","n_good_bins":int(sum(good))}
            quiet=float(np.quantile(y[good],.3))
            mad=float(1.4826*np.median(np.abs(y[good]-np.median(y[good]))))
            threshold=max(3*max(quiet,0),quiet+5*mad)
            flare=good & (y>=threshold)
            rec={"status":"source_masked_particle_flare_control", "filename":name,
                 "n_good_bins":int(np.count_nonzero(good)),
                 "n_high_particle_bins":int(np.count_nonzero(flare)),
                 "flare_rate_threshold_cps":threshold,"baseline_rate_cps":quiet,
                 "mjdref_tt":cv.mjd_reference(h.header),
                 "hduclas2":h.header.get("HDUCLAS2"),
                 "energy_band":"pipeline particle-flare background, not independent 0.2-2keV source aperture",
                 "source_masked_regions":len(hd["REGION"].data) if "REGION" in hd and hd["REGION"].data is not None else None}
        return rec, (t[good],y[good],flare[good],rec["mjdref_tt"])
    except Exception as e:
        return {"status":"background_fits_parse_error","error":repr(e)[:250]}, None

def background_peak_comparison(analysis, bg):
    if bg is None:
        return {"status":"no_independent_particle_curve","all_peaks_clear":None}
    t,y,flare,ref=bg
    pairs=[]
    max_spacing=max(float(np.nanmedian(np.diff(np.sort(t))))*1.6,180)
    for peak in analysis.get("peaks",[]):
        at=(peak["mjd_tt"]-ref)*86400
        nearest=int(np.argmin(abs(t-at)))
        close=bool(abs(t[nearest]-at)<=max_spacing)
        pairs.append({"peak_mjd_tt":peak["mjd_tt"],"nearest_bg_separation_s":float(abs(t[nearest]-at)),
                      "bg_cps":float(y[nearest]) if close else None,
                      "source_masked_particle_flare_coincident":bool(flare[nearest]) if close else None})
    return {"status":"compared_source_masked_particle_flare_field",
            "n_source_peaks_tested":len(pairs),"all_peaks_clear":all(
                x["source_masked_particle_flare_coincident"] is False for x in pairs) if pairs else None,
            "peak_checks":pairs}

def background_only_finder(curve):
    # Source BACK1-3 is extracted/scaled from an off-source aperture inside
    # the source product. It is not the independent field particle curve.
    d=dict(curve)
    d["rate"]=np.asarray(curve["background"],float).copy()
    d["background"]=np.zeros_like(d["rate"])
    d["error"]=np.maximum(np.asarray(curve["error"],float),1e-4)
    result=cv.find_eruptions(d,50814.,(.5,24.))
    return {"method":"frozen finder on scaled source-specific off-source background aperture",
            "independent_from_source_counts":True,
            "independent_field_particle_GTIs":False,
            "recurrence_unvetted":bool(result.get("range_match") and result.get("n_peaks",0)>=2),
            "n_peaks":result.get("n_peaks",0)}

def poisson_injection(curve, detid):
    # Count-domain surrogate (independently draw gross and scaled background
    # counts). No event positions/epatplot, hence not event-list completeness.
    seed=int(hashlib.sha256(("dr14-pps-inject-v1:"+detid).encode()).hexdigest()[:16],16)
    rng=np.random.default_rng(seed)
    valid=np.asarray(curve["good"])&np.isfinite(curve["rate"])&np.isfinite(curve["error"])&np.isfinite(curve["background"])
    ix=np.flatnonzero(valid)
    if len(ix)<35:return {"status":"insufficient_good_bins","seed":seed}
    centers=np.array([i for i in ix if i>=5 and i+5<len(valid) and valid[i-4:i+5].all()],int)
    if len(centers)<2:return {"status":"no_two_exposed_centers","seed":seed}
    i=int(rng.choice(centers))
    delta=abs(curve["t_s"][centers]-curve["t_s"][i])
    others=centers[(delta>=2*3600)&(delta<=22*3600)]
    if not len(others):return {"status":"no_exposed_2_to_22h_pair","seed":seed}
    j=int(rng.choice(others))
    t=np.asarray(curve["t_s"],float)
    exp=np.maximum(np.asarray(curve["exposure"],float),1)
    bkg=np.maximum(np.nan_to_num(curve["background"],nan=0),0)
    quiet=np.asarray(curve["rate"],float)[ix]
    quiet=quiet[quiet<=np.quantile(quiet,.55)]
    base=np.full(len(t),0.,float)
    base[ix]=np.maximum(rng.choice(quiet,size=len(ix)),0)
    errs=np.asarray(curve["error"],float)[ix]
    trials=[]
    for ampl_sigma in (8.,16.):
        amp=ampl_sigma*max(float(np.median(errs)),.001)
        pulse=amp*(np.exp(-.5*((t-t[i])/850)**2)+np.exp(-.5*((t-t[j])/850)**2))
        pulse[~valid]=0
        gross=rng.poisson(np.clip((base+bkg+pulse)*exp,0,1e9))
        background=rng.poisson(np.clip(bkg*exp,0,1e9))
        injected=(gross-background)/exp
        error=np.sqrt(gross+background+1)/exp
        work=dict(curve)
        work["rate"]=np.where(valid,injected,np.nan)
        work["error"]=np.where(valid,error,np.nan)
        out=cv.find_eruptions(work,50814.,(.5,24.))
        detected=[(p["mjd_tt"]-50814.)*86400 for p in out.get("peaks",[])]
        m1=any(abs(x-t[i])<=900 for x in detected)
        m2=any(abs(x-t[j])<=900 for x in detected)
        trials.append({"amplitude_sigma":ampl_sigma,"width_s":850,
                       "input_peak_count":2,"matched_peak_count":int(m1)+int(m2),
                       "recovered_train":bool(m1 and m2 and out.get("range_match"))})
    return {"status":"count_domain_surrogate_complete","seed":seed,
            "n_trials":len(trials),"n_recovered":sum(x["recovered_train"] for x in trials),
            "limitations":"Not a Poisson event-list spatial/background-area injection or independent survey sensitivity",
            "trials":trials}

def scan_observation(obsid, rows):
    record={"obs_id":obsid,"n_frozen_detections":len(rows),"started_utc":utc(),
            "source_results":[{"detid":r["detid"],"obs_id":obsid,"srcid":r["srcid"],
               "ra":r["ra"],"dec":r["dec"],"pps_srcnum":r["pps_srcnum"],
               "pn_pileup_diagnostic":r["pn_pileup_diagnostic"],
               "status":"not_screened"} for r in rows],
            "file_receipts":[]}
    attempts=listing(obsid,retries=3)
    record["directory_receipts"]=attempts
    if "http_status" not in attempts[-1] or attempts[-1].get("entries",0)==0:
        record["status"]="mirror_PPS_directory_unavailable"
        return record
    # Parse only names whose exact observation prefix and FITS extension are verified.
    # listing() reports samples, so fetch directory list anew for the full filenames.
    import html
    from urllib.parse import urljoin, unquote, urlparse
    url=HOST+obsid+"/PPS/"
    try:
        with SESSION.get(url,timeout=(12,30)) as response:
            response.raise_for_status()
            rawlist=response.text
        matches=re.findall(r'href\s*=\s*["\x27]([^"\x27]+)',rawlist,flags=re.I)
        names=[]
        for href in matches:
            full=urljoin(url,html.unescape(href))
            part=urlparse(full)
            if part.netloc.lower()!="heasarc.gsfc.nasa.gov" or not part.path.startswith("/FTP/xmm/data/rev0/"+obsid+"/PPS/"):
                continue
            name=unquote(part.path.rsplit("/",1)[-1])
            if re.fullmatch("P"+obsid+r"PN[A-Z0-9_.+-]{8,90}\.FTZ",name,re.I):
                names.append(name)
        names=sorted(set(names))
    except Exception as exc:
        record["status"]="mirror_PPS_listing_parse_failure"
        record["listing_error"]=repr(exc)[:260]
        return record
    pn=[x for x in names if "SRCTSR8" in x.upper() and
        x.upper().startswith("P"+obsid+"PNX000")]
    bgfiles=[x for x in names if "FBKTSR" in x.upper()]
    record["pn_merged_srctsr8_available"]=len(pn)
    record["pn_fbktsr_available"]=len(bgfiles)
    record["pn_per_exposure_src_tsr_only"]=not bool(pn)
    bgdata=None
    if bgfiles:
        chosen=sorted(bgfiles,key=lambda x:("PNX000" not in x,x))[0]
        raw=get_file(obsid,chosen,record["file_receipts"])
        if raw is not None:
            b=flare_background(raw,chosen)
            if isinstance(b,tuple):
                record["field_particle_background"],bgdata=b
            else:
                record["field_particle_background"]=b
        else:
            record["field_particle_background"]={"status":"download_failed"}
    else:
        record["field_particle_background"]={"status":"no_source_masked_pn_FBKTSR"}
    if not pn:
        record["status"]="no_merged_soft_band_pn_products"
        return record
    # Original PPS index is an optimisation, NOT identity evidence. Every
    # downloaded light curve must pass the measured FITS SRC_RA/DEC match.
    expected={f"{int(r['pps_srcnum']):03X}" for r in rows if r.get("pps_srcnum","").isdigit()}
    preferred=[n for n in pn if re.search(r"SRCTSR8([0-9A-F]{3})\.FTZ$",n,re.I) and
               re.search(r"SRCTSR8([0-9A-F]{3})\.FTZ$",n,re.I).group(1).upper() in expected]
    alternatives=[n for n in pn if n not in preferred]
    pick=preferred+alternatives[:max(0, min(25,len(rows)*8)-len(preferred))]
    record["n_source_products_scheduled"]=len(pick)
    record["possible_uninspected_source_products"]=len(pn)-len(pick)
    for name in pick:
        if all(r["status"]=="position_verified_net_pn_soft_screened" for r in record["source_results"]):
            break
        raw=get_file(obsid,name,record["file_receipts"],attempts=3)
        if raw is None: continue
        try:
            position,curve,ref,header=net_merged(raw)
            distances=sorted((cv.sep_arcsec(position,(float(r["ra"]),float(r["dec"]))),idx)
                             for idx,r in enumerate(rows))
            nearest,idx=distances[0]
            if nearest > 15 or (len(distances)>1 and distances[1][0]-nearest<2):
                record.setdefault("unmatched_fits_sources",[]).append({
                    "filename":name,"source_position_deg":list(position),
                    "nearest_sep_arcsec":round(nearest,3),"ambiguous":bool(
                        len(distances)>1 and distances[1][0]-nearest<2)})
                continue
            row=record["source_results"][idx]
            if row["status"]=="position_verified_net_pn_soft_screened": continue
            finder=cv.find_eruptions(curve,ref,(0.5,24.0))
            bgcomp=background_peak_comparison(finder,bgdata)
            srcbg=background_only_finder(curve)
            row.update(status="position_verified_net_pn_soft_screened",
                 source_sep_arcsec=round(nearest,3),pn_product=name,
                 merged_pn_fits_sha256=hashlib.sha256(raw).hexdigest(),
                 source_fits_header=header,
                 good_300s_bins=int(np.sum(curve["good"])),
                 finder=finder,
                 unvetted_repeated_peak_trigger=bool(
                     finder.get("n_peaks",0)>=2 and finder.get("range_match")),
                 independent_field_particle_background=bgcomp,
                 offsource_scaled_background_control=srcbg,
                 poisson_count_domain_injection=poisson_injection(curve,row["detid"]),
                 pileup_SAS_epatplot="NOT PERFORMED",
                 independent_soft_background_region_event_extraction="NOT PERFORMED",
                 stellar_blazar_multiepoch_novelty_vetting="NOT PERFORMED")
        except Exception as exc:
            record.setdefault("product_parse_failures",[]).append({
                "filename":name,"error":repr(exc)[:300]})
    n=sum(x["status"]=="position_verified_net_pn_soft_screened" for x in record["source_results"])
    record["pn_detections_screened"]=n
    record["status"]="has_verified_pn_curves" if n else "no_position_verified_soft_band_pn_source_curve"
    record["finished_utc"]=utc()
    return record

def main():
    config=json.loads((BASE/"control_batch.json").read_text())
    assert config["stage"]=="heasarc_batch"
    batch=config["batch_id"]
    assert re.fullmatch(r"[a-zA-Z0-9_]+",batch)
    start=int(config["start"])
    count=int(config["count"])
    assert 0<=start<670 and 1<=count<=120
    cohort,fr=manifest()
    # The three pre-registered pilot targets precede chronological ordering;
    # all remaining observations are deterministic, no new-data reprioritization.
    pilot=["0872392901","0306630101","0200730401"]
    ordered=[x for x in pilot if x in cohort]+sorted(set(cohort)-set(pilot))
    chosen=ordered[start:start+count]
    receipt={"as_of_utc":utc(),"batch_id":batch,"start":start,"count_requested":count,
             "count_selected":len(chosen),"cohort_manifest_sha256":fr["final_csv_sha256"],
             "cohort_obsids":670,"cohort_provisional_detections":967,
             "selection":"three pre-registered pilots, then ascending 10-digit OBSID",
             "thresholds_sha256":cv.checksum(BASE/"control_thresholds.json"),
             "method":"unchanged control_validation.find_eruptions on PN merged NET RATE1+2+3 at 300s",
             "mirror":HOST,"observations":[],"status":"RUNNING"}
    filename=f"heasarc_batch_{batch}.json"
    t0=time.monotonic()
    for obsid in chosen:
        if time.monotonic()-t0>CUTOFF:
            receipt["status"]="TIME_BUDGET_PAUSE_NOT_COMPLETE"
            break
        try:
            row=scan_observation(obsid,cohort[obsid])
        except Exception as exc:
            row={"obs_id":obsid,"status":"uncaught_observation_error","error":repr(exc)[:400],
                 "source_results":[{"detid":r["detid"],"status":"not_screened"} for r in cohort[obsid]]}
        receipt["observations"].append(row)
        valid=[s for x in receipt["observations"] for s in x.get("source_results",[])
               if s.get("status")=="position_verified_net_pn_soft_screened"]
        flagged=[s for s in valid if s.get("unvetted_repeated_peak_trigger")]
        receipt.update(obsids_attempted=len(receipt["observations"]),
                       obsids_with_screened_pn=len({s["obs_id"] for s in valid if "obs_id" in s}
                           or {x["obs_id"] for x in receipt["observations"] if any(
                               s.get("status")=="position_verified_net_pn_soft_screened"
                               for s in x.get("source_results",[]))}),
                       pn_detections_screened=len(valid),
                       unvetted_recurrent_triggers=len(flagged),
                       candidates_validated=0,
                       independent_soft_region_and_event_pileup_gate=False,
                       archive_670_observation_complete=False,
                       ready_for_review=False)
        save(filename,receipt)
        print("HEASARC_BATCH",batch,obsid,row["status"],
              row.get("pn_detections_screened",0),"cumulative",len(valid),
              "flags",len(flagged),flush=True)
    receipt["status"]="BATCH_FINISHED_UNVETTED" if receipt["status"]=="RUNNING" else receipt["status"]
    receipt["finished_utc"]=utc()
    receipt["elapsed_s"]=round(time.monotonic()-t0,2)
    receipt["limitations"]=[
        "HEASARC public mirror covers a US archive subset; missing PPS cannot be treated as a negative",
        "PPS merged PN 0.2-2 keV FITS-only screen: lack of merged RATE1-3 is not a non-eruption",
        "Source-masked off-source FBKTSR is a hard particle-flare screen, not a matched 0.2-2 keV source-aperture background",
        "Poisson rate/count-domain surrogates are not SAS event-list injections or independent completeness",
        "No SAS epatplot pile-up, energy spectral temperature evolution, MOS corroboration or trial-adjusted archive FPR",
        "Stellar/SIMBAD/4LAC and eRASS/Swift/Chandra/TNS/ATel novelty vetting deferred for every new trigger",
        "SDSS-footprint provisional hosts are not complete all-sky galaxy-nucleus selection"]
    save(filename,receipt)
    print(json.dumps({x:receipt.get(x) for x in
        ("status","obsids_attempted","pn_detections_screened","unvetted_recurrent_triggers","candidates_validated")},indent=2))

if __name__=="__main__": main()
