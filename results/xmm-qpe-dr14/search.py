#!/usr/bin/env python3
"""4XMM-DR14 reproducible detection screen. No discoveries without host/flare validation."""
import csv
import datetime as dt
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import shutil
import tarfile
import tempfile
import traceback
from urllib.parse import urlencode
import numpy as np
import requests
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.ndimage import median_filter
from scipy.signal import find_peaks

ROOT = Path(__file__).parent
OUT = ROOT / "outputs"
OUT.mkdir(exist_ok=True)
MAX_PRODUCTS = 12
CAT_MD5 = "47c6bd113cb8510db418665472c81b48"
CAT_URLS = [
    "https://heasarc.gsfc.nasa.gov/FTP/xmm/data/catalogues/4XMM_DR14cat_v1.0.fits.gz",
    "https://xmmssc.irap.omp.eu/Catalogue/4XMM-DR14/4XMM_DR14cat_v1.0.fits.gz",
]
# J2000 positions, from public SIMBAD / QPE discovery papers; no position invented.
CONTROLS = {
    "GSN 069": (19.786098, -34.191813),
    "RX J1301.9+2747": (195.50054, 27.78267),
    "eRO-QPE1": (37.946917, -10.336197),
    "eRO-QPE2": (38.704042, -44.325458),
    "eRO-QPE3": (210.2222, -28.7665),
    "AT2019qiz": (71.657833, -10.226367),
}
RECEIPTS = []
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "Open-Space/4XMM-DR14-archive-research"})
SUMMARY = {"as_of_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "dataset": "4XMM-DR14 v1.0 (2024-07-09)",
           "result": "inconclusive: no validated new QPE",
           "limits": [], "checks": {}, "recovery_fraction": None,
           "false_positive_rate": None, "validated_new_qpe": 0}
def write_json(name, obj):
    (OUT / name).write_text(json.dumps(obj, indent=2, sort_keys=True, default=str) + "\n")
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()
def download(urls, path, max_bytes, expected_md5=None):
    for url in urls:
        rec = {"url": url, "time_utc": dt.datetime.now(dt.timezone.utc).isoformat()}
        try:
            with SESSION.get(url, stream=True, timeout=(20, 90)) as r:
                rec.update(http_status=r.status_code, final_url=r.url,
                           content_type=r.headers.get("Content-Type"),
                           content_length=r.headers.get("Content-Length"))
                r.raise_for_status()
                n, h, m = 0, hashlib.sha256(), hashlib.md5()
                if int(r.headers.get("Content-Length", "0")) > max_bytes:
                    raise ValueError("file exceeds configured byte budget")
                with path.open("wb") as f:
                    for block in r.iter_content(1024 * 1024):
                        if not block:
                            continue
                        n += len(block)
                        if n > max_bytes:
                            raise ValueError("stream exceeds configured byte budget")
                        h.update(block)
                        m.update(block)
                        f.write(block)
                rec.update(bytes=n, sha256=h.hexdigest(), md5=m.hexdigest())
                if expected_md5 and m.hexdigest() != expected_md5:
                    raise ValueError("catalogue MD5 != published checksum")
                if n < 1024:
                    raise ValueError("implausibly short data file")
                rec["verified"] = bool(expected_md5)
                RECEIPTS.append(rec)
                return path
        except Exception as exc:
            rec["error"] = f"{type(exc).__name__}: {str(exc)[:320]}"
            RECEIPTS.append(rec)
            path.unlink(missing_ok=True)
    return None
def finite(v, default=np.nan):
    try:
        x = float(v)
        return x if math.isfinite(x) else default
    except (TypeError, ValueError):
        return default
def array(d, field, default=np.nan):
    return np.asarray(d[field]) if field in d.names else np.full(len(d), default)
def txt(v):
    return v.decode(errors="replace").strip() if isinstance(v, bytes) else str(v).strip()
def controls_match(d):
    ra, dec = np.asarray(array(d,"RA"),float), np.asarray(array(d,"DEC"),float)
    result = {}
    for name,(cra,cdec) in CONTROLS.items():
        rough = np.where((np.abs((ra-cra+180)%360-180)<0.03) & (np.abs(dec-cdec)<0.03))[0]
        if len(rough):
            true_sep=SkyCoord(ra[rough]*u.deg,dec[rough]*u.deg).separation(
                    SkyCoord(cra*u.deg,cdec*u.deg)).arcsec
            result[name]=[int(i) for i in rough[true_sep<=12]]
        else:
            result[name]=[]
    return result
def peaks(t, y, err):
    good=np.isfinite(t)&np.isfinite(y)&np.isfinite(err)&(err>0)
    t,y,err=t[good],y[good],err[good]
    if len(t)<18:
        return {"status":"insufficient_good_bins","bins":len(t)}
    order=np.argsort(t)
    t,y,err=t[order],y[order],err[order]
    span=float((t[-1]-t[0])/3600)
    cadence=float(np.nanmedian(np.diff(t)))
    if cadence<=0 or span<3:
        return {"status":"short_or_invalid","span_h":span}
    quiet=float(np.nanpercentile(y,25))
    noise=max(float(np.nanmedian(err)),float(1.4826*np.nanmedian(
              np.abs(y-np.nanmedian(y)))),1e-5)
    z=median_filter(y,size=3,mode="nearest")
    threshold=max(4*max(quiet,0),quiet+5*noise)
    ix,_=find_peaks(z,height=threshold,prominence=max(3*noise,.5*threshold),distance=3)
    ix=[int(p) for p in ix if p>=2 and p<len(t)-2 and
        not np.any(np.diff(t[p-2:p+3])>2.5*cadence) and
        (y[p]-quiet)/max(err[p],1e-5)>=5]
    intervals=(np.diff(t[ix])/3600).tolist() if len(ix)>1 else []
    repeat=len(ix)>=2 and any(.5<=x<=24 for x in intervals)
    return {"status":"unvetted_repeated_peak_flag" if repeat else "no_repeated_peak_flag",
            "bins":len(t),"span_h":round(span,3),"cadence_s":round(cadence,1),
            "baseline_rate":round(quiet,5),"n_peaks":len(ix),
            "peak_times_h":[round(float((t[p]-t[0])/3600),3) for p in ix[:15]],
            "peak_rates":[round(float(y[p]),5) for p in ix[:15]],
            "intervals_h":[round(x,3) for x in intervals[:14]],
            "peak_to_quiet":round(float(np.max(y[ix])/max(quiet,1e-4)),2) if ix else None}
def parse_lc(data):
    try:
        with fits.open(io.BytesIO(data),memmap=False) as hd:
            for ext in hd[1:]:
                if ext.data is None or not hasattr(ext.data,"names"):
                    continue
                cols={x.upper():x for x in ext.data.names}
                if "TIME" not in cols or "RATE" not in cols:
                    continue
                d=ext.data
                t=np.asarray(d[cols["TIME"]],float)
                y=np.asarray(d[cols["RATE"]],float)
                err=np.asarray(d[cols["ERROR"]],float) if "ERROR" in cols else np.sqrt(np.maximum(y,0))
                if "FRACEXP" in cols:
                    ok=np.asarray(d[cols["FRACEXP"]],float)>=.7
                    t,y,err=t[ok],y[ok],err[ok]
                background={}
                for key in ("BACKV","BACK_RATE","BKG_RATE","BG_RATE"):
                    if key in cols:
                        background["background_column"]=key
                        background["background_median"]=finite(np.nanmedian(d[cols[key]]))
                        break
                return peaks(t,y,err),background
    except Exception as exc:
        return {"status":"parse_error","detail":str(exc)[:180]},{}
    return {"status":"no_time_rate_columns"},{}
def fetch_pps(r):
    params={"obsno":r["obs_id"],"level":"PPS","extension":"FTZ",
            "instname":"PN","name":"SRCTSR"}
    url="https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?"+urlencode(params)
    with tempfile.TemporaryDirectory() as td:
        tarpath=Path(td)/"archive.tar"
        if download([url],tarpath,120*1024*1024) is None:
            return {"status":"unavailable"}
        products=[]
        try:
            with tarfile.open(tarpath,"r:*") as tar:
                for member in tar:
                    if not member.isfile():
                        continue
                    nm=Path(member.name).name.upper()
                    if "SRCTSR" not in nm or "PN" not in nm[10:15]:
                        continue
                    if not nm.endswith((".FTZ",".FITS",".FIT")) or member.size>20*1024*1024:
                        continue
                    # The source-specific PPS product index is hexadecimal (not catalog SRC_NUM).
                    if f"{int(r['pps_srcnum']):03X}" not in nm[-7:]:
                        continue
                    data=tar.extractfile(member).read()
                    if nm.endswith(".FTZ"):
                        data=gzip.decompress(data)
                    flagged,bg=parse_lc(data)
                    products.append({"name":nm,"sha256":hashlib.sha256(data).hexdigest(),
                                     "lightcurve":flagged,**bg})
                    if len(products)>=12:
                        break
        except Exception as exc:
            return {"status":"tar_or_product_error","detail":str(exc)[:200]}
        return {"status":"products_parsed" if products else "no_matching_pn_srctsr",
                "products":products}
def catalog():
    with tempfile.TemporaryDirectory() as td:
        gz=Path(td)/"4xmm.fits.gz"
        raw=Path(td)/"4xmm.fits"
        if download(CAT_URLS,gz,1250*1024*1024,CAT_MD5) is None:
            SUMMARY["limits"].append("Full release data not retrieved/verified: no observed detection denominator")
            return
        with gzip.open(gz,"rb") as source,raw.open("wb") as out:
            shutil.copyfileobj(source,out,1024*1024)
        with fits.open(raw,memmap=True) as hd:
            tbl=next(x for x in hd[1:] if isinstance(x,fits.BinTableHDU))
            d=tbl.data
            n=len(d)
            names=list(d.names)
            write_json("schema.json",{"n_detections":n,"columns":names,
                                       "uncompressed_bytes":raw.stat().st_size})
            ts=np.asarray(array(d,"TSERIES",0),bool)
            flag=np.asarray(array(d,"SUM_FLAG",9),float)
            cts=np.asarray(array(d,"EP_8_CTS"),float)
            exp=np.nanmax(np.stack([np.asarray(array(d,f"{c}_ONTIME"),float)
                  for c in ("PN","M1","M2")]),axis=0)
            hr3=np.asarray(array(d,"EP_HR3"),float)
            hr4=np.asarray(array(d,"EP_HR4"),float)
            prob=np.asarray(array(d,"EP_CHI2PROB"),float)
            var=np.asarray(array(d,"EP_FVAR"),float)
            good=ts&(flag<=1)&(cts>=250)&(exp>=21600)&(hr3<-.3)&(hr4<-.4)
            indices=np.flatnonzero(good)
            SUMMARY["checks"]["full_catalogue_rows"]=n
            SUMMARY["checks"]["tseries_rows"]=int(ts.sum())
            SUMMARY["checks"]["eligible_soft_long_detections"]=len(indices)
            SUMMARY["checks"]["catalogue_verified_md5"]=CAT_MD5
            def value(i,key):
                return txt(d[key][i]) if key in names else ""
            def record(i):
                p=finite(prob[i],1.)
                v=finite(var[i],0.)
                return {"detid":value(i,"DETID"),"srcid":value(i,"SRCID"),
                        "obs_id":value(i,"OBS_ID").zfill(10),
                        "src_num":value(i,"SRC_NUM"),
                        "pps_srcnum":value(i,"PPS_SRCNUM") or value(i,"SRC_NUM"),
                        "iauname":value(i,"IAUNAME"),
                        "ra":round(finite(d["RA"][i]),6),
                        "dec":round(finite(d["DEC"][i]),6),
                        "cts":round(finite(cts[i]),1),
                        "ontime_s":round(finite(exp[i]),1),
                        "hr3":round(finite(hr3[i]),3),
                        "hr4":round(finite(hr4[i]),3),
                        "chi2prob":p,"fvar":v,
                        "pn_pileup":value(i,"PN_PILEUP"),
                        "score":round(float(-math.log10(max(p,1e-15))+2*min(max(v,0),5)),3)}
            eligible=[record(int(i)) for i in indices]
            eligible.sort(key=lambda r:(r["detid"],r["obs_id"]))
            fieldnames=list(eligible[0]) if eligible else ["detid","obs_id"]
            with (OUT/"frozen_detections.csv").open("w",newline="") as f:
                writer=csv.DictWriter(f,fieldnames=fieldnames)
                writer.writeheader()
                writer.writerows(eligible)
            top=sorted(eligible,key=lambda r:(-r["score"],r["detid"]))[:8]
            matched=controls_match(d)
            ctrls={}
            for name,ids in matched.items():
                rows=[record(i) for i in ids if ts[i]]
                ctrls[name]=sorted(rows,key=lambda r:(-r["ontime_s"],r["detid"]))[:2]
            requested=top.copy()
            for name,rows in ctrls.items():
                if rows:
                    requested.append({"control":name,**rows[0]})
            seen=set()
            unique=[]
            for r in requested:
                key=(r["obs_id"],r["pps_srcnum"])
                if key not in seen:
                    unique.append(r)
                    seen.add(key)
            unique=unique[:MAX_PRODUCTS]
            frozen={"release":"4XMM-DR14 v1.0","published_md5":CAT_MD5,
                    "full_catalogue_detections":n,"eligible_count":len(eligible),
                    "full_manifest_sha256":sha256(OUT/"frozen_detections.csv"),
                    "screen_cuts":"TSERIES; SUM_FLAG<=1; EP_8_CTS>=250; camera ontime>=21600s; EP_HR3<-.3; EP_HR4<-.4",
                    "ranking":"-log10(max(EP_CHI2PROB,1e-15)) + 2*clip(EP_FVAR,0,5)",
                    "ranked_first_8":top,"known_qpe_matches":ctrls,
                    "pps_requests":unique,"max_pps_requests":MAX_PRODUCTS,
                    "script_sha256":sha256(Path(__file__)),
                    "note":"X-ray prefilter only; galaxy nuclei NOT classified"}
            write_json("frozen_selection.json",frozen)
            SUMMARY["checks"]["frozen_manifest_sha256"]=frozen["full_manifest_sha256"]
            SUMMARY["checks"]["control_catalogue_matches"]={name:len(ids) for name,ids in matched.items()}
        write_json("receipts.json",RECEIPTS)
        examined=[]
        for r in unique:
            examined.append({**r,"pps":fetch_pps(r)})
            write_json("lightcurve_checks.json",examined)
            write_json("receipts.json",RECEIPTS)
        flagged=[r for r in examined if any(
            p["lightcurve"].get("status")=="unvetted_repeated_peak_flag"
            for p in r["pps"].get("products",[]))]
        SUMMARY["checks"]["pps_requests"]=len(examined)
        SUMMARY["checks"]["pps_detections_with_parsed_lightcurve"]=sum(
            bool(r["pps"].get("products")) for r in examined)
        SUMMARY["checks"]["unvetted_repeat_flags"]=len(flagged)
        SUMMARY["checks"]["known_qpe_flare_recovery"]="NOT VALIDATED"
        write_json("unvetted_flags.json",flagged)
        SUMMARY["limits"]+= [
            "This is only an X-ray prefilter: no Legacy Surveys/SDSS/2MASS galaxy-centroid or Gaia stellar rejection",
            "At most 12 detections were requested from PPS; incomplete sampling cannot yield an archive-wide null",
            "Full source-vs-background light curves, soft-band recurrence, source pile-up, spectral temperature evolution not validated",
            "No actual control-QPE flare recovery or injection-based completeness; no empirical false-positive rate",
            "No eRASS1/Swift/Chandra epoch check or ATel/TNS/literature clearance of any unvetted flags",
        ]
def main():
    try:
        catalog()
    except Exception as exc:
        SUMMARY["limits"].append("Pipeline exception: "+repr(exc)[:250])
        write_json("error.json",{"message":repr(exc),"traceback":traceback.format_exc()[-5000:]})
    finally:
        if "pps_requests" not in SUMMARY["checks"]:
            SUMMARY["limits"].append("No time-series analysis completed; null is NOT supported")
        write_json("receipts.json",RECEIPTS)
        write_json("summary.json",SUMMARY)
        lines=["# No validated new QPE (bounded archive screen)","",
               "Run UTC: "+SUMMARY["as_of_utc"],"",
               "## Executed checks","",
               "~~~json",json.dumps(SUMMARY["checks"],indent=2),"~~~","",
               "## Unresolved limitations",""]
        lines+=["- "+x for x in SUMMARY["limits"]]
        lines+=["","The source list is frozen before looking at PPS light curves.",
                "No discovery label or remarkable-findings entry is justified.",""]
        (OUT/"RUN_REPORT.md").write_text("\n".join(lines))
        print(json.dumps({"result":SUMMARY["result"],
                          "checks":SUMMARY["checks"],
                          "limits":SUMMARY["limits"]},indent=2))
if __name__=="__main__":
    # Explicitly separate no-network control replay from control downloads.
    # An archive stage is intentionally NOT exposed by this dispatch.
    from pathlib import Path as _Path
    import json as _json
    _mode = _json.loads((_Path(__file__).resolve().parent / "control_batch.json").read_text()).get("stage")
    if _mode == "replay":
        from replay_controls import main as _replay_main
        _replay_main()
    elif _mode == "archive_probe":
        from archive_probe import main as _probe_main
        _probe_main()
    elif _mode == "archive_hosts":
        from archive_hosts import main as _hosts_main
        _hosts_main()
    elif _mode == "archive_freeze":
        from archive_freeze import main as _freeze_main
        _freeze_main()
    elif _mode == "archive_diagnostics":
        from archive_diagnostics import main as _diag_main
        _diag_main()
    elif _mode == "archive_pilot":
        from archive_pilot import main as _pilot_main
        _pilot_main()
    else:
        from control_validation import main as _controls_main
        _controls_main()
