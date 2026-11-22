"""Audit independently corroborated Ly-alpha spectra with native red-side imaging."""
import argparse
from collections import Counter
from copy import deepcopy
import csv
import gzip
from hashlib import file_digest
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from acquisition import Spectra, prepare
from evidence_summary import adjudications, frozen_digest, ledger, stable_key
from legacy_bricks import legacy_brick_cones
from forced_legacy import forced_photometry
from photometry import bound, continuum_upper, imaging_ew, LEGACY_R, NMAGGY_FNU, C_A_PER_S, response
from vetting import (Cache,Ancillary,atom,number,validate_local,enrich,cube_test,tap_batch,
                     external_check,ordinary_agn_disfavored,AGN_GAP,AGN_PASSED,PAIR_EXTERNAL,PAIR_REPEAT)
import spectral as spec

Z_FILTERS=("MzLS-z","decam2014-z","decamDR1-z")
EW_GAP="redward imaging does not independently constrain EW above 240 A"
DIRECT_GAP="direct imaging checks unavailable"


def write(handle,obj):
    handle.write(json.dumps(atom(obj),allow_nan=False,separators=(",",":"))+"\n")


def digest(path):
    with path.open("rb") as source:
        return file_digest(source,"sha256").hexdigest()


def metric(x):
    value=number(x)
    return value if math.isfinite(value) else None


def support(p,v):
    passed=set(v.get("passed") or ())
    z=p.get("z_lya")
    if z is None:
        return False
    prior=v.get("DESI") or {}
    vac=(number(prior.get("dex_flag"))==1 and number(prior.get("target_sep_arcsec"),math.inf)<=1.25
         and number(prior.get("vi_quality"))>=3 and abs(number(prior.get("vi_z"))-z)<.008)
    external=(v.get("SPARCL") or {}).get("matches",[])
    pair=bool(passed & {PAIR_EXTERNAL,PAIR_REPEAT})
    anchor=vac or pair or any(number(m.get("redshift_warning"))==0 and
                               abs(number(m.get("redshift"))-z)<.006 for m in external)
    instrument=[m.get("diagnostics",{}) for m in external]
    if number(prior.get("dex_flag"))==1 and number(prior.get("target_sep_arcsec"),math.inf)<=1.25:
        instrument.append(v.get("DESI_lines") or {})
    line=any(number(d.get("lines",{}).get("primary",{}).get("snr"))>=4 for d in instrument)
    repeat=any(number(r.get("snr"))>=4 and str(r.get("shotid"))!=str(p.get("shotid"))
               for r in v.get("repeats",[]))
    # Historical external adjudications may predate a propagated VAC field.
    # Respect their already recorded confirmations for eligibility; recompute
    # current VAC association and all other gates before any promotion.
    anchor |= bool(passed & {"independent DESI visual redshift","independent survey redshift agrees"})
    line |= "independent instrument confirms line" in passed
    repeat |= "independent HETDEX shot" in passed
    return bool(anchor and (line or repeat or pair))


def bands_for(cap,band):
    if band=="r":
        return ("BASS-r",) if cap=="north" else ("decam2014-r","decamDR1-r") if cap=="south" else LEGACY_R
    return ("MzLS-z",) if cap=="north" else ("decam2014-z","decamDR1-z") if cap=="south" else Z_FILTERS


def interval(p,origin,filters,flux,noise,lf,le,aperture,transmission=1.,upper_line=None,psf_uncertainty=0.):
    """Positive lower and conditional upper bounds; units are physical cgs and AB nmgy."""
    lam=number(p.get("wave_vac"),p["wave"])
    try:
        upper_ph=(max(0.,flux)+3*noise)*(1+psf_uncertainty)/transmission
        b=bound(p,origin,filters,upper_ph,lf,le,aperture)
        coeff,_=continuum_upper(lam,1.,filters)
        if not b or not b["eligible"] or coeff is None:
            return None
        b.update(flux_nanomaggy=flux,noise_nanomaggy=noise,transmission=transmission,
                 continuum_noise=coeff*noise*(1+psf_uncertainty)/transmission,ew0_upper_A=None,
                 aperture_fractional_uncertainty=psf_uncertainty)
        if upper_line is not None and flux>3*noise:
            # Maximize predicted photon counts per unit line-centre continuum
            # over the allowed continuum slopes, even in the blue filter tail.
            limits=[]
            for name in filters:
                w,r,_=response(name)
                envelope=np.maximum(1.,(w/lam)**-3)
                limits.append(.95*max(0.,1-psf_uncertainty)*(flux-3*noise)*NMAGGY_FNU*C_A_PER_S*
                              trapezoid(r/w,w)/trapezoid(r*w*envelope,w))
            c_lo=min(limits)
            if c_lo>0:
                b["ew0_upper_A"]=float(upper_line/(c_lo*(1+p["z_lya"])))
                b["continuum_lower_cgs_per_A"]=float(c_lo)
        return b
    except (ValueError,KeyError,OSError) as exc:
        return dict(eligible=False,origin=origin,error=str(exc))


def measure(p,v,cone,native):
    cube=v.get("cube") or {}
    cen=cube.get("central_aperture") or {}
    fc=number(cen.get("value"),0.)*1e-17
    ec=number(cen.get("scatter"),0.)*1e-17
    lam=number(p.get("wave_vac"),p["wave"])
    # The total PSF line is a deliberately generous numerator for upper-EW
    # diagnostics; a finite narrowband or continuum aperture does not bound
    # arbitrary low-surface-brightness emission outside that aperture.
    sigma=number(p.get("sigma_A"),math.inf)
    window=math.erf(6./(math.sqrt(2)*sigma)) if sigma>0 else 0
    upper=(1.1*max((fc+3*ec)/window,p["flux"]+3*p["flux_err"])
           if window>.5 and fc>0 and ec>0 else None)
    tests=[]
    cat_tests=[]
    ap_tests=[]
    matches=cone.get("matches",[])
    cap=cone.get("selected_cap")
    if len(matches)==1 and number(matches[0].get("sep_arcsec"),math.inf)<=1.25:
        obj=matches[0]
        for band in ("r","z"):
            noise=number(obj.get("flux_err_"+band))
            flux=metric(obj.get("flux_"+band))
            trans=number(obj.get("mw_transmission_"+band))
            good=(obj.get("type")!="DUP" and obj.get("maskbits")==0 and obj.get("allmask_"+band)==0 and
                  number(obj.get("nobs_"+band))>0 and 0<=number(obj.get("fracmasked_"+band))<.15 and
                  0<=number(obj.get("fracflux_"+band))<.5 and 0<trans<=1.001 and noise>0 and flux is not None)
            if good:
                b=interval(p,"native Tractor "+band,bands_for(cap,band),flux,noise,p["flux"],p["flux_err"],
                           "total optical / PSF-integrated line",transmission=trans)
                if b:
                    b.update(band=band,tractor_id=obj.get("tractor_id"))
                    tests.append(b)
                    if b["eligible"]:
                        cat_tests.append(b)
    for band,m in native.items():
        if band not in ("r","z") or not isinstance(m,dict) or not m.get("eligible"):
            continue
        flux=metric(m.get("flux"))
        err=number(m.get("flux_err"))
        sep=3600*math.hypot(((number(m.get("ra"),math.inf)-number(cube.get("ra"),0)+180)%360-180)*
                            math.cos(math.radians(p["dec"])),number(m.get("dec"),math.inf)-number(cube.get("dec"),0))
        if (not math.isfinite(sep) or sep>.3 or flux is None or err<=0 or fc<=0 or ec<=0 or
            abs(number(m.get("aperture_radius_arcsec"),0)-number(cube.get("aperture_radius_arcsec"),0))>.05 or
            number(m.get("blank_apertures"))<16 or number(cen.get("n"))<16):
            continue
        frac=number((m.get("psf") or {}).get("aperture_fractional_uncertainty"),math.inf)
        if not math.isfinite(frac) or frac<0:
            continue
        b=interval(p,"native coadd "+band,bands_for(cap,band),flux,err,fc,ec,"PSF-matched 1.5 arcsec radius",
                   upper_line=upper,psf_uncertainty=frac)
        if b:
            b.update(band=band,center_offset_arcsec=sep)
            tests.append(b)
            if b["eligible"]:
                ap_tests.append(b)
    old=(v.get("imaging_EW") or {}).get("bounds",[])
    other=[]
    for t in old:
        if not t.get("eligible") or metric(t.get("ew0_lower_A")) is None:
            continue
        if "Legacy Tractor" in t.get("origin","") and cat_tests:
            continue
        if t.get("error_r_nanomaggy") and t.get("continuum_upper_cgs_per_A") and t.get("r_upper_nanomaggy_before_cal_margin"):
            t=dict(t,continuum_noise=t["error_r_nanomaggy"]*t["continuum_upper_cgs_per_A"]/
                   (1.05*t["r_upper_nanomaggy_before_cal_margin"]))
            ap_tests.append(t)
        else:
            other.append(t)
    # If the native r image could establish 240 A at zero target signal, it
    # is the primary native band regardless of the actual target measurement.
    # A z exposure then checks only the two-band conditional upper bound.
    deep_r=any(t.get("origin")=="native coadd r" and
               t["line_flux_lower_cgs"]>240*3*1.05*t["continuum_noise"]*(1+p["z_lya"]) for t in ap_tests)
    ap_eligible=[t for t in ap_tests if not (deep_r and t.get("origin")=="native coadd z")]
    chosen=[min(group,key=lambda b:b["continuum_noise"]) for group in (cat_tests,ap_eligible) if group]+other
    lower=min((t["ew0_lower_A"] for t in chosen),default=None)
    # Only count a model exclusion when both clean coadd bands independently
    # exclude 240, preventing one-line contamination of a single broadband filter.
    uppers={t["band"]:t.get("ew0_upper_A") for t in ap_tests if t.get("band") in ("r","z")}
    upper_two=max(uppers.values()) if len(uppers)==2 and all(x is not None for x in uppers.values()) else None
    return dict(ew0_lower_A=lower,ew0_upper_A_two_band_model=upper_two,bounds=tests,deep_native_r=deep_r,
                selected=[t["origin"] for t in chosen],previous_bounds=old,
                rule="Beta -3 to 0, redward AB throughput >=99%, photometry ±3σ with 5% calibration margin; "
                     "line ±3σ with 10% margin. Select by continuum noise before flux within total and aperture "
                     "routes; retain the smallest eligible lower bound across routes. Upper is conditional on "
                     "isolated continuum-dominated bands, compact line emission, and the reported apertures.")


def can_promote(p,v,img):
    return (number(img.get("ew0_lower_A"))>240 and support(p,v) and
            (p.get("stream")=="imaging_rescue" or number(p.get("ew0_lower"))>240))


def full_check(p,v,cone,img,cache,cat,aux):
    """Re-evaluate all source, literature, independent line and active-nucleus gates."""
    m=cat.metadata(p["row_index"])
    rec={**m,**p,"catalog":m}
    if not rec.get("shuffled_detail_complete") and not rec.get("lya_family"):
        rec.update(enrich(rec,cat))
    local=validate_local(rec,cat,aux)
    if local['disposition']!="unresolved":
        return local
    if aux.errors.get("LAN"):
        local["LAN_catalog_unavailable"]=True
    local["Legacy"]=cone
    for name in ("NED","SIMBAD"):
        attempt=tap_batch(name,[rec])[0]
        local[name]=attempt if attempt.get("status")=="queried" else v.get(name,attempt)
    rec["rescue_audit"]={k:v[k] for k in ("cube","HSC","SPARCL","Legacy_forced") if v.get(k)}
    checked=external_check(rec,local,cache)
    checked["imaging_EW"]=img
    if checked["disposition"]=="excluded":
        return checked
    remaining=set(checked["reasons"])
    if number(img.get("ew0_lower_A"))>240:
        remaining.discard(EW_GAP)
    if cone.get("status")=="queried" or any(t.get("eligible") for t in img["bounds"]):
        remaining.discard(DIRECT_GAP)
    if not remaining-{AGN_GAP}:
        if ordinary_agn_disfavored(rec,checked,p["z_lya"]):
            remaining.discard(AGN_GAP)
            checked["passed"]=list(dict.fromkeys([*checked["passed"],AGN_PASSED,"native red-side EW bound exceeds 240 A"]))
        else:
            remaining.add(AGN_GAP)
    checked.update(disposition="unresolved" if remaining else "candidate",reasons=sorted(remaining))
    return checked


def null_rebuild(cat,p):
    if p.get("stream")=="nuisance_shuffle":
        from nuisance_control import _draw,_model
        from repeats import spectra_for_rows
        design=p["mock_design"]
        f,e=spectra_for_rows(cat,[p["row_index"],*design["donor_rows"]])
        models=[_model(y,se) for y,se in zip(f,e)]
        y,se=_draw(models[0],models[1:],e[0],design["draw_seed"])
    else:
        from scan import BATCH
        i=p["row_index"]
        start=i//BATCH*BATCH
        f,e=cat.read(start,min(start+BATCH,cat.size))
        y,se=spec.shuffle_fibres(f,e,seed=410087+start)
        y,se=y[i-start],se[i-start]
    return spec.detail(cat.wave,y,se,p,density_unit=1e-17)


def clone_docket(n,old,path,line):
    return dict(first_pass_line=n,previous_ledger=path,previous_ledger_line=line,
                previous_disposition=old["vet"]["disposition"],
                previous_snapshot_sha256=frozen_digest(old["peak"],old["vet"]).hex(),
                key=list(stable_key(old["peak"])),peak=old["peak"],vet=old["vet"])


def existing(cat,indir,source,counts):
    if source=="observed":
        candidates=adjudications(indir,include_ew_imaging=False)
    else:
        candidates=((n,r,"vetted.jsonl.gz",r["vet"]["disposition"]) for n,r in ledger(indir/source/"vetted.jsonl.gz"))
    seen=set()
    out=[]
    for n,(line,entry,path,_initial) in enumerate(candidates,1):
        p,v=entry["peak"],entry["vet"]
        if p["channel"]!="ew":
            continue
        counts["ew_rows_in_register"]+=1
        seen.add((p["row_index"],p["wave"]))
        if v["disposition"]!="unresolved" or not support(p,v):
            continue
        counts["independently_corroborated_eligible"]+=1
        d=clone_docket(n,entry,path,line)
        if source!="observed":
            d["mock_details"]=null_rebuild(cat,p)
        out.append(d)
    return out,seen


def beyond(cat,aux,indir,seen,counts):
    """Inspect every frozen imaging-prescreen survivor outside the vetting ledger."""
    selected=[]
    for n,row in ledger(indir/"ew-rescue-spectra.jsonl.gz"):
        if row.get("status")!="science_peak_for_imaging":
            continue
        counts["scout_prescreen_science_peaks"]+=1
        if (row["row_index"],row["wave"]) in seen:
            continue
        counts["outside_existing_register"]+=1
        i=row["row_index"]
        m=cat.metadata(i)
        p=dict(row,ra=m.get("ra_det") if m.get("ra_det") is not None else row["ra"],
               dec=m.get("dec_det") if m.get("dec_det") is not None else row["dec"],
               channel="ew",stream="imaging_rescue",release=cat.release)
        # A distinct visit is required by the repeat prior; a verified DESI
        # crossmatch can also confirm the primary with one HETDEX observation.
        rec={**m,**p,"catalog":m}
        rec.update(enrich(rec,cat))
        v=validate_local(rec,cat,aux)
        if v["disposition"]!="unresolved" or not support(rec,v):
            counts["outside_failing_local_or_corroboration"]+=1
            continue
        counts["outside_scout_with_independent_support"]+=1
        selected.append(dict(source_ledger="ew-rescue-spectra.jsonl.gz",source_line=n,
                             key=list(stable_key(p)),peak=p,vet=v))
    return selected


def trim(cache):
    """Remove only downloaded native imaging files from this temporary cache."""
    root=cache.root/"legacy-tractor"
    if not root.exists():
        return
    files=sorted((p for p in root.rglob("*") if p.is_file() and ".fits" in p.name and not p.name.endswith(".json")),
                 key=lambda p:p.stat().st_mtime)
    used=sum(p.stat().st_size for p in files)
    for p in files:
        if used<1_500_000_000:
            break
        if "survey-bricks" in p.name:
            continue
        sz=p.stat().st_size
        p.unlink()
        used-=sz


def decision(rec,v,cone,cache,cat,aux,counts,cached_bands=()):
    """Apply the identical native continuum and external gates to data and controls."""
    if "cube" not in v or not v["cube"].get("central_aperture"):
        meta=cat.metadata(rec["row_index"])
        try:
            v["cube"]=cube_test({**meta,**rec},cache.fetch)
        except Exception as exc:
            v["cube"]={"status":"unavailable","error":str(exc)}
    def get_band(band):
        cube=v["cube"]
        for line,measurement in cached_bands:
            item=measurement.get(band) or {}
            psf=item.get("psf") or {}
            if (item.get("eligible") and cube.get("central_aperture") and
                (cube.get("ra"),cube.get("dec"))==(rec["ra"],rec["dec"]) and
                (item.get("ra"),item.get("dec"))==(rec["ra"],rec["dec"]) and
                item.get("aperture_radius_arcsec")==cube.get("aperture_radius_arcsec")==1.5 and
                psf.get("cube_seeing_arcsec")==cube.get("seeing_arcsec") and
                psf.get("cube_moffat_beta")==cube.get("two_psf_beta")):
                counts["archived_band_measurements_reused"]+=1
                return {band:dict(item,reused_from_ew_imaging_native_line=line)}
        return forced_photometry(rec,cone,cache,cube=cube,bands=(band,))
    native=get_band("r")
    img=measure(rec,v,cone,native)
    if not (img["deep_native_r"] and number(img["ew0_lower_A"])>240):
        native.update(get_band("z"))
        img=measure(rec,v,cone,native)
    lower=img["ew0_lower_A"]
    upper=img["ew0_upper_A_two_band_model"]
    counts["native_cones_queried"]+=cone.get("status")=="queried"
    counts["aperture_r_eligible"]+=bool((native.get("r") or {}).get("eligible"))
    counts["aperture_z_eligible"]+=bool((native.get("z") or {}).get("eligible"))
    counts["with_EW_lower"]+=lower is not None
    counts["lower_above_240"]+=number(lower)>240
    counts["two_band_model_upper_below_240"]+=upper is not None and upper<240
    if can_promote(rec,v,img):
        v=full_check(rec,v,cone,img,cache,cat,aux)
        counts["full_literature_agn_rechecks"]+=1
    else:
        v["imaging_EW"]=img
        reasons=set(v.get("reasons") or [])
        if not reasons:
            reasons.add(EW_GAP)
        if cone.get("status")=="queried":
            reasons.discard(DIRECT_GAP)
        if upper is not None and upper<240:
            reasons.add("both matched red bands put the compact continuum-model rest EW below 240 A")
        if number(lower)<=240:
            reasons.add(EW_GAP)
        v.update(disposition="unresolved",reasons=sorted(reasons))
    v["imaging_EW"]=img
    return v,native,img


def run(indir,out,release,cache_root):
    out.mkdir(parents=True,exist_ok=True)
    cf=("sample","source_register","source_line","row_index","source_id","detectid","ra","dec","z_lya","wave",
        "flux","flux_err","rest_EW_spectral_lower_A","rest_EW_imaging_lower_A","rest_EW_conditional_upper_A",
        "locus","legacy_outcome","forced_r","forced_z","disposition","reasons")
    summary={"release":release,"input_sha256":digest(indir/"vetted.jsonl.gz"),
             "additional_input_sha256":digest(indir/"ew-rescue-spectra.jsonl.gz"),"counts":{}}
    path=prepare(release,cache_root,out)
    cache=Cache(cache_root)
    aux=Ancillary(cache)
    try:
      with Spectra(path,release) as cat, (out/"outcomes.csv").open("w",newline="") as table, \
           gzip.open(out/"native.jsonl.gz","wt",compresslevel=3) as native_log:
        writer=csv.DictWriter(table,cf)
        writer.writeheader()
        native_line=0
        for sample in ("observed","matched-shuffle-vet","shuffle-vet"):
            counts=Counter()
            regular,seen=existing(cat,indir,sample,counts)
            adds=beyond(cat,aux,indir,seen,counts) if sample=="observed" else []
            if sample!="observed":
                summary.setdefault("control_sha256",{})[sample]=digest(indir/sample/"vetted.jsonl.gz")
            rows=regular+adds
            dest=out if sample=="observed" else out/sample
            dest.mkdir(exist_ok=True)
            rows.sort(key=lambda x:(round(x["peak"]["ra"],1),round(x["peak"]["dec"],1),x["peak"]["ra"]))
            with gzip.open(dest/"vetted.jsonl.gz","wt",compresslevel=3) as vetted, \
                 gzip.open(dest/"additional.jsonl.gz","wt",compresslevel=3) as extra:
              for start in range(0,len(rows),4):
                batch=rows[start:start+4]
                coords=sorted({(d["peak"]["ra"],d["peak"]["dec"]) for d in batch})
                cones=dict(zip(coords,legacy_brick_cones(coords,cache),strict=True))
                for docket in batch:
                    p=docket["peak"]
                    v=deepcopy(docket["vet"])
                    rec={**p,**docket.pop("mock_details",{})}
                    if sample!="observed":
                        rec["shuffled_detail_complete"]=True
                    cone=cones[p["ra"],p["dec"]]
                    v,native,img=decision(rec,v,cone,cache,cat,aux,counts)
                    lower=img["ew0_lower_A"]
                    upper=img["ew0_upper_A_two_band_model"]
                    # Permanent input references are written once in this ledger
                    # per assessed acquisition/PSF, never as lower-limit proxies.
                    native_line+=1
                    write(native_log,dict(sample=sample,row_index=p["row_index"],wave=p["wave"],
                                          cone=cone,bands=native,imaging=img))
                    v["native_ew_photometry"]=dict(ledger="native.jsonl.gz",line=native_line,
                        observed_footprint=cone.get("outcome"),two_band_model_upper_A=upper)
                    docket.update(vet=v)
                    destination=extra if "source_ledger" in docket else vetted
                    write(destination,docket)
                    counts["additional" if destination is extra else "registered"]+=1
                    counts["disposition_"+v["disposition"]]+=1
                    if destination is extra:
                        counts["additional_"+v["disposition"]]+=1
                    writer.writerow(dict(sample=sample,source_register=docket.get("source_ledger","vetted.jsonl.gz"),
                        source_line=docket.get("source_line",docket.get("first_pass_line")),row_index=p["row_index"],
                        source_id=p.get("source_id"),detectid=p.get("detectid"),ra=p["ra"],dec=p["dec"],z_lya=p["z_lya"],
                        wave=p["wave"],flux=p["flux"],flux_err=p["flux_err"],rest_EW_spectral_lower_A=p.get("ew0_lower"),
                        rest_EW_imaging_lower_A=lower,rest_EW_conditional_upper_A=upper,
                        locus="outside scout" if destination is extra else "existing ledger",
                        legacy_outcome=cone.get("outcome"),forced_r=(native.get("r") or {}).get("status"),
                        forced_z=(native.get("z") or {}).get("status"),disposition=v["disposition"],
                        reasons=" | ".join(v.get("reasons",[]))))
                trim(cache)
            summary["counts"][sample]=dict(counts)
            print(release,sample,dict(counts),flush=True)
    finally:
        aux.close()
    with (out/"outcomes.csv").open() as infile:
        rows=list(csv.DictReader(infile))
    observed=[r for r in rows if r["sample"]=="observed"]
    informative=[r for r in observed if metric(r["rest_EW_imaging_lower_A"]) is not None]
    upper_rows=[r for r in observed if metric(r["rest_EW_conditional_upper_A"]) is not None]
    best={"passed":[r for r in observed if r["disposition"]=="candidate"],
          "largest_measured_lower_bounds":sorted(informative,key=lambda r:-float(r["rest_EW_imaging_lower_A"]))[:15],
          "smallest_two_band_conditional_upper_bounds":sorted(upper_rows,key=lambda r:float(r["rest_EW_conditional_upper_A"]))[:15]}
    (out/"best.json").write_text(json.dumps(best,indent=2,allow_nan=False)+"\n")
    summary["interpretation"]="Independent identification and line corroboration are required before the native imaging test. Two positive band-specific EW upper limits below 240 exclude the compact/aperture broadband-continuum model only. Unknown foreground or extended emission is retained as unresolved."
    (out/"summary.json").write_text(json.dumps(atom(summary),indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("input",type=Path)
    ap.add_argument("output",type=Path)
    ap.add_argument("--release",required=True,choices=("sc1","sc2"))
    ap.add_argument("--cache",required=True,type=Path)
    ap.add_argument("--scout-control-only",action="store_true")
    args=ap.parse_args()
    if args.scout_control_only:
        from ew_scout_control import run_control
        run_control(args.input,args.output,args.release,args.cache)
    else:
        run(args.input,args.output,args.release,args.cache)
