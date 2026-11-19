#!/usr/bin/env python3
"""One manually dispatched DJA v4.4 z>9 data-gated run, no discovery claim.

Results are unvetted 1D model proxies; downloaded 2D images are coadds,
not disjoint original exposures. This script does not make NIRCam or
literature novelty claims.
"""
import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.request import Request, urlopen

import numpy as np
from astropy.io import fits
import select_targets as selector
import download_targets as targeted

HERE = Path(__file__).resolve().parent
FROZEN = HERE.parent / "jwst-dja-nirspec-highz-screen"
SOURCES = json.loads((FROZEN / "frozen_sources.json").read_text())
CONTROLS = json.loads((FROZEN / "literature_baseline.json").read_text())["published_positive_controls"]
CANUCS = {"abell370-v4","macs0416-v4","macs0417-v4","macs1149-v4","macs1423-v4"}
FITCOL = ("catalog_row_index root file ra dec grating catalog_z_best catalog_z_phot "
          "published_controls download_sha256 download_bytes valid_bins two_sided "
          "z16_high_conditional z50_high_conditional z84_high_conditional "
          "p_z_gt12_proxy p_z_gt14_proxy p_lowz_gal_proxy p_dwarf_proxy "
          "p_artifact_proxy priority_class figure two_d_status status error").split()

def now():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temp=path.with_name(path.name+".partial")
    temp.write_text(json.dumps(data, indent=2, sort_keys=True)+"\n")
    temp.replace(path)

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1<<20),b""):
            h.update(chunk)
    return h.hexdigest()

def item_for(row):
    root, name, url = targeted.safe(row)
    if root in CANUCS:
        url=f"https://s3.amazonaws.com/grizli-canucs/nirspec/{root}/{name}"
    return dict(root=root,file=name,url=url)

def get(url, path, max_mib):
    path.parent.mkdir(parents=True,exist_ok=True)
    part=path.with_name(path.name+".partial")
    part.unlink(missing_ok=True)
    size=0
    try:
        req=Request(url,headers={"User-Agent":"Open-Space-DJA-v4.4-audit/1.0",
                                 "Accept":"application/octet-stream"})
        with urlopen(req,timeout=45) as response, part.open("wb") as out:
            status=response.status
            final_url=response.url
            if status!=200:
                raise ValueError(f"HTTP status {status}")
            while True:
                block=response.read(1<<20)
                if not block: break
                size+=len(block)
                if size>max_mib*(1<<20):
                    raise ValueError(f"File exceeded cap of {max_mib} MiB")
                out.write(block)
        if size<2880:
            raise ValueError("Empty or too-short response body")
        part.replace(path)
        return dict(http_status=status,final_url=final_url,bytes=size,sha256=sha(path))
    finally:
        part.unlink(missing_ok=True)

def summary_file(work, receipt):
    item=SOURCES["dja"]["catalog"]
    path=work/item["filename"]
    urls=[item["url"],
          "https://s3.amazonaws.com/msaexp-nirspec/extractions/"+item["filename"]]
    attempts=[]
    for url in urls:
        try:
            info=get(url,path,210)
            sha256,md5=selector.hashes(path)
            if md5.lower()!=item["publisher_md5"].lower():
                raise ValueError(f"Publisher MD5 mismatch: {md5} != {item['publisher_md5']}")
            attempts.append(dict(url=url,**info,observed_md5=md5,
                                 expected_md5=item["publisher_md5"],verified=True))
            receipt.update(summary_attempts=attempts,summary_sha256=sha256,
                           summary_md5_verified=True)
            return path
        except Exception as exc:
            attempts.append(dict(url=url,utc=now(),error=type(exc).__name__+": "+str(exc)[:350]))
            path.unlink(missing_ok=True)
    receipt.update(summary_attempts=attempts,summary_md5_verified=False)
    raise RuntimeError("All DJA v4.4 catalog GETs failed or returned nonmatching MD5")

def load_model():
    spec=importlib.util.spec_from_file_location("dja_pr195_screen",FROZEN/"screen.py")
    model=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    model.HI=np.round(np.arange(9.,18.001,.1),2)
    return model

def arrays_1d(hdul):
    if "SPEC1D" not in hdul:
        raise ValueError("No SPEC1D FITS extension")
    table=hdul["SPEC1D"].data
    names={name.lower():name for name in hdul["SPEC1D"].columns.names}
    for col in ("wave","flux"):
        if col not in names:
            raise ValueError(f"No SPEC1D.{col}")
    wave=np.asarray(table[names["wave"]],float).ravel()
    flux=np.asarray(table[names["flux"]],float).ravel()
    if "full_err" in names:
        err=np.asarray(table[names["full_err"]],float).ravel()
    elif "err" in names:
        raw=np.asarray(table[names["err"]],float).ravel()
        err=np.hypot(raw,.02*flux)
    else:
        raise ValueError("No FITS error column")
    if len(wave)!=len(flux) or len(err)!=len(wave) or len(wave)<20:
        raise ValueError("Malformed 1D array lengths")
    valid=np.isfinite(wave)&np.isfinite(flux)&np.isfinite(err)&(err>0)&(wave>=.55)&(wave<=5.6)
    if "valid" in names:
        valid &= np.asarray(table[names["valid"]],bool).ravel()
    if "npix" in names:
        valid &= np.asarray(table[names["npix"]],float).ravel()>0
    if valid.sum()<20:
        raise ValueError(f"Insufficient valid 1D samples: {valid.sum()}")
    if not np.all(np.diff(wave[valid])>0):
        raise ValueError("Nonmonotonic wavelength")
    return wave,flux,err,valid

def fit(hdul,row,model):
    wave,flux,err,valid=arrays_1d(hdul)
    bank=model.bank(wave,"prism" if row.get("grating","").lower()=="prism" else "grating")
    y=np.where(valid,flux,0.)[:,None]
    sigma=np.where(valid,err,1.)[:,None]
    result=model.fit(wave,y,sigma,valid[:,None],bank)
    metrics=dict(valid_bins=int(valid.sum()),two_sided=bool(result["two_sided"][0]))
    for key,original in [
      ("z16_high_conditional","q16_high"),("z50_high_conditional","q50_high"),
      ("z84_high_conditional","q84_high"),("p_z_gt12_proxy","p_gt12_proxy"),
      ("p_z_gt14_proxy","p_gt14_proxy"),("p_lowz_gal_proxy","p_lowz_gal_proxy"),
      ("p_dwarf_proxy","p_dwarf_proxy"),("p_artifact_proxy","p_artifact_proxy")]:
        metrics[key]=float(result[original][0])
        if not np.isfinite(metrics[key]):
            raise ValueError("Nonfinite numerical model value "+key)
    return metrics,(wave,flux,err,valid)

def plot_2d(hdul,arrays,z,filename,label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    wave,flux,err,valid=arrays
    sci=np.asarray(hdul["SCI"].data,float) if "SCI" in hdul and hdul["SCI"].data is not None else None
    if sci is None or sci.ndim!=2 or min(sci.shape)<2:
        return "NATIVE_2D_SCI_MISSING"
    wht=np.asarray(hdul["WHT"].data,float) if "WHT" in hdul and hdul["WHT"].data is not None else None
    sn=(sci*np.sqrt(np.maximum(wht,0))) if wht is not None and wht.shape==sci.shape else None
    # Interpolated wavelengths approximate the detector x-coordinate only;
    # retain this limitation on each saved figure.
    left,right=float(wave[valid].min()),float(wave[valid].max())
    edge=.121567*(1+z)
    fig,ax=plt.subplots(3,1,figsize=(9,7),layout="constrained")
    for panel,image,title in [(ax[0],sci,"DJA SCI coadd (not independent exposures)"),
                              (ax[1],sn,"DJA SCI × sqrt(WHT); approximate S/N")]:
        if image is None:
            panel.text(.5,.5,"WHT unavailable",transform=panel.transAxes,ha="center")
        else:
            ok=image[np.isfinite(image)]
            vmax=max(float(np.nanpercentile(np.abs(ok),97)),1e-15) if ok.size else 1.
            panel.imshow(image,origin="lower",aspect="auto",interpolation="nearest",
                         cmap="coolwarm",vmin=-vmax,vmax=vmax,
                         extent=(left,right,0,sci.shape[0]))
        panel.axvline(edge,color="gold",ls="--")
        panel.set_ylabel("spatial pixel")
        panel.set_title(title,fontsize=9)
    ax[2].errorbar(wave[valid],flux[valid],yerr=err[valid],
                   fmt=".",alpha=.6,ms=2,lw=.5)
    ax[2].axvline(edge,color="gold",ls="--")
    ax[2].set_ylabel("1D flux (µJy)")
    ax[2].set_xlabel("Observed wavelength (µm; 2D x-axis approximate)")
    fig.suptitle(f"{label}: z={z:.3f} — UNVETTED, NO NIRCam/independent DQ",fontsize=10)
    filename.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(filename,dpi=80)
    plt.close(fig)
    return "DJA_COADD_2D_PLOTTED_NOT_INDEPENDENT_VISIT_REVIEW"

def download(item,work):
    record=dict(url=item["url"],root=item["root"],file=item["file"],utc=now())
    dest=work/"fits"/item["root"]/item["file"]
    try:
        info=get(item["url"],dest,512)
        with fits.open(dest,memmap=False) as hd:
            hd.verify("exception")
            if "SPEC1D" not in hd: raise ValueError("FITS has no SPEC1D")
        return dict(record,**info,status="VALID_FITS_SPEC1D_DOWNLOADED")
    except Exception as exc:
        dest.unlink(missing_ok=True)
        return dict(record,status="FAILED_GET_OR_FITS",error=type(exc).__name__+": "+str(exc)[:350])

def write_csv(path,records,columns):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=columns,extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

def compact_audit(source,dest):
    columns=["dja_catalog_row_index","root","file","ra","dec","grating",
             "z_best","z_phot","targeted_selected","highz_hints","published_controls"]
    with source.open(newline="") as inp,dest.open("w",newline="") as out:
        writer=csv.DictWriter(out,fieldnames=columns,extrasaction="ignore")
        writer.writeheader()
        count=0
        for row in csv.DictReader(inp):
            writer.writerow(row);count+=1
    if count!=80367:
        raise ValueError(f"Audit denominator {count} != 80367")

def run(work,out):
    work.mkdir(parents=True,exist_ok=True)
    out.mkdir(parents=True,exist_ok=True)
    receipt={"started_at_utc":now(),"status":"STARTED_NO_SCIENCE_CLAIM","stage":"catalog",
             "release":"DJA v4.4, Zenodo 15472354",
             "model":"PR #195 exploratory Lyman+wing vs dusty Balmer, line, dwarf, artifact; uncalibrated",
             "nircam_independently_validated":0,"literature_independently_reviewed":0,
             "fully_validated_discoveries":0}
    write_json(out/"run_receipt.json",receipt)
    try:
        catalog=summary_file(work,receipt)
        receipt["stage"]="selection"
        write_json(out/"run_receipt.json",receipt)
        selection=selector.select(catalog,work/"selection")
        compact_audit(work/"selection"/"catalog_audit.csv",out/"compact_catalog_index.csv")
        shutil.copy2(work/"selection"/"selection_receipt.json",out/"selection_receipt.json")
        shutil.copy2(work/"selection"/"selected_spectra.csv",out/"selected_spectra.csv")
        receipt.update(selected_rows=selection["selected_spectrum_rows"],
                       unmatched_published_controls=selection["published_controls_unmatched"],
                       stage="individual_fits_downloads")
        write_json(out/"run_receipt.json",receipt)
        with (work/"selection"/"selected_spectra.csv").open(newline="") as f:
            selected=list(csv.DictReader(f))
        unique={}
        for row in selected:
            item=item_for(row)
            unique.setdefault((item["root"],item["file"]),item)
        items=list(unique.values())
        receipt["distinct_fits_urls"]=len(items)
        records=[]
        with ThreadPoolExecutor(max_workers=6) as pool:
            futures=[pool.submit(download,item,work) for item in items]
            for future in as_completed(futures):
                records.append(future.result())
        records.sort(key=lambda r:(r["root"],r["file"]))
        with (out/"download_receipts.jsonl").open("w") as f:
            for record in records:
                f.write(json.dumps(record,sort_keys=True)+"\n")
        good={(r["root"],r["file"]):r for r in records
              if r["status"]=="VALID_FITS_SPEC1D_DOWNLOADED"}
        receipt.update(downloaded_fits=len(good),download_failed=len(items)-len(good),
                       stage="1d_model_and_2d_figures")
        write_json(out/"run_receipt.json",receipt)
        model=load_model() if good else None
        fits_rows=[];review=[];made_figures=0;missing_2d=0;fitted=0
        controls_by_name={x["name"]:x for x in CONTROLS}
        for row in selected:
            key=(row["root"],row["file"])
            record={col:"" for col in FITCOL}
            record.update(catalog_row_index=row.get("dja_catalog_row_index",""),root=row["root"],
                          file=row["file"],ra=row.get("ra",""),dec=row.get("dec",""),
                          grating=row.get("grating",""),catalog_z_best=row.get("z_best",""),
                          catalog_z_phot=row.get("z_phot",""),
                          published_controls=row.get("published_controls",""),
                          priority_class=("KNOWN_PUBLISHED_CONTROL" if row.get("published_controls")
                                          else "UNREVIEWED_NOT_NEW"))
            data=good.get(key)
            if data is None:
                record["status"]="FITS_MISSING_OR_FAILED"
                fits_rows.append(record)
                continue
            record.update(download_sha256=data["sha256"],download_bytes=data["bytes"])
            try:
                with fits.open(work/"fits"/row["root"]/row["file"],memmap=False) as hd:
                    metrics,arrays=fit(hd,row,model)
                    record.update(metrics)
                    fitted+=1
                    try: catz=float(row.get("z_best",""))
                    except (TypeError,ValueError): catz=float("nan")
                    ctrl=bool(row.get("published_controls"))
                    # Union of published/catalog high-z and proxy z>12 support;
                    # NOT a detection or a calibrated P(z>12).
                    quick=(metrics["z50_high_conditional"]>12 and
                           metrics["p_z_gt12_proxy"]>=.05)
                    if ctrl or (np.isfinite(catz) and catz>12) or quick:
                        z=catz if np.isfinite(catz) and catz>12 else metrics["z50_high_conditional"]
                        if ctrl and (not np.isfinite(catz) or catz<=12):
                            name=row["published_controls"].split(";")[0].split(":")[0]
                            z=float(controls_by_name[name]["z_spec"])
                        label=row["root"]+"/"+row["file"]
                        fig_id=hashlib.sha256(label.encode()).hexdigest()[:18]
                        fig=out/"figures"/(fig_id+".png")
                        two_status=plot_2d(hd,arrays,z,fig,row["file"])
                        record["two_d_status"]=two_status
                        if fig.exists():
                            record["figure"]=str(fig.relative_to(out))
                            made_figures+=1
                        else:
                            missing_2d+=1
                        review.append({"root":row["root"],"file":row["file"],
                                       "ra":row.get("ra"),"dec":row.get("dec"),
                                       "dja_z_best":row.get("z_best"),
                                       "fit_highz_conditional_median":metrics["z50_high_conditional"],
                                       "p_gt12_proxy_uncalibrated":metrics["p_z_gt12_proxy"],
                                       "published_control_match":row.get("published_controls"),
                                       "priority":record["priority_class"],
                                       "status_2d":two_status,
                                       "nircam_status":"UNVALIDATED",
                                       "priority_search_status":"PENDING ADS/arXiv/DJA/JADES/CEERS/UNCOVER/MoM/GLASS"})
                    record["status"]="UNVETTED_1D_MODEL_ONLY"
            except Exception as exc:
                record["status"]="FITTING_OR_PLOTTING_FAILED"
                record["error"]=type(exc).__name__+": "+str(exc)[:350]
            fits_rows.append(record)
        write_csv(out/"fit_results.csv",fits_rows,FITCOL)
        write_json(out/"zgt12_literature_review_queue.json",review)
        receipt.update(fit_succeeded=fitted,fit_failed_or_missing=len(selected)-fitted,
                       z_gt12_review_rows=len(review),two_d_figures=made_figures,
                       missing_2d_sci=missing_2d,finished_at_utc=now(),stage="finished",
                       independently_literature_reviewed=0,independently_nircam_validated=0)
        if not good:
            receipt["status"]="BLOCKED_NO_INDIVIDUAL_FITS"
        elif len(good)!=len(items) or fitted!=len(selected):
            receipt["status"]="PARTIAL_DOWNLOAD_OR_FIT_FAILURE_NO_NULL"
        else:
            receipt["status"]="ONE_D_FITS_COMPLETED_VALIDATION_PENDING_NO_DISCOVERY"
        write_json(out/"run_receipt.json",receipt)
        return 0 if receipt["status"]=="ONE_D_FITS_COMPLETED_VALIDATION_PENDING_NO_DISCOVERY" else 2
    except Exception as exc:
        receipt.update(status="BLOCKED_OR_FAILED_NO_NULL",failure=type(exc).__name__+": "+str(exc)[:1000],
                       traceback=traceback.format_exc(limit=5)[-2500:],
                       finished_at_utc=now())
        write_json(out/"run_receipt.json",receipt)
        print(receipt["status"],receipt["failure"],file=sys.stderr)
        return 2

def self_test():
    from astropy.io.fits import PrimaryHDU,BinTableHDU,ImageHDU,HDUList,Column
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        w=np.linspace(.6,5.2,180)
        f=np.where(w<.121567*14.,0.,.7)
        tab=BinTableHDU.from_columns([
            Column(name="wave",format="D",array=w),
            Column(name="flux",format="D",array=f),
            Column(name="err",format="D",array=np.ones_like(f)*.1)],name="SPEC1D")
        img=np.tile(f,(12,1))
        path=Path(tmp)/"synthetic.fits"
        HDUList([PrimaryHDU(),tab,ImageHDU(img,name="SCI"),
                 ImageHDU(np.ones_like(img)*100,name="WHT")]).writeto(path)
        with fits.open(path) as hd:
            metrics,arrays=fit(hd,{"grating":"PRISM"},load_model())
            assert np.isfinite(metrics["p_z_gt12_proxy"])
            p=Path(tmp)/"figure.png"
            assert plot_2d(hd,arrays,13.,p,"synthetic").startswith("DJA_COADD")
            assert p.is_file()
    print("PASS fake FITS read, exploratory fit and 2D figure (not an observation)")

if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",type=Path);p.add_argument("--out",type=Path)
    p.add_argument("--self-test",action="store_true")
    args=p.parse_args()
    if args.self_test:self_test()
    else:
        if not args.work or not args.out:p.error("--work and --out needed")
        sys.exit(run(args.work,args.out))
