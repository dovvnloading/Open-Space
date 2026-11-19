#!/usr/bin/env python3
"""Batch-complete HI-beam Euclid DR1-F follow-up, network only on Actions.

Direct IRSA SIA -> vetted public Euclid S3 MER FITS byte ranges, plus IRSA
HATS cone searches, DESI DR1 spectra and Legacy DR10 photometric redshifts.
Every unobserved/failed service is an UNKNOWN, never a nondetection.
"""
from __future__ import annotations
import argparse, csv, gzip, hashlib, io, json, math, os, re, time
from pathlib import Path
from urllib.parse import urlsplit
import numpy as np
import pandas as pd
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.ndimage import gaussian_filter, binary_dilation

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
BATCH=OUT/"beam_batches"
SIA="https://irsa.ipac.caltech.edu/SIA"
TAP="https://datalab.noirlab.edu/tap/sync"
BUCKET="nasa-irsa-euclid-dr1"
HATS_PREFIX="s3://nasa-irsa-euclid-dr1/contributed/dr1/hats/reformatted_catalogs/"
BANDS=("VIS","Y","J","H")
HEADERS={"User-Agent":"Open-Space-Euclid-HI-all-beam-followup/2026.11"}
C=299792.458
RE=(2.,5.,10.,20.)
KNOWN_AGC={"229383","229384","229385","322753"}
def dig(b):return hashlib.sha256(b).hexdigest()
def fnum(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None
def jsonwrite(p,obj):
    p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_suffix(p.suffix+".tmp")
    tmp.write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
    tmp.replace(p)
def sep(ra,dec,ra2,dec2):
    a,b,c,d=map(math.radians,[ra,dec,ra2,dec2])
    z=math.sin((b-d)/2)**2+math.cos(b)*math.cos(d)*math.sin((a-c)/2)**2
    return 7200*math.degrees(math.asin(min(1.,math.sqrt(max(0.,z)))))
def posbox(ra,dec,radius,racol="ra",deccol="dec"):
    rr=radius/3600
    half=rr/max(.01,math.cos(math.radians(dec)))
    lo=(ra-half)%360;hi=(ra+half)%360
    ras=(f"{racol} BETWEEN {lo:.8f} AND {hi:.8f}" if lo<=hi
        else f"({racol}>={lo:.8f} OR {racol}<={hi:.8f})")
    return f"({deccol} BETWEEN {max(-90,dec-rr):.8f} AND {min(90,dec+rr):.8f}) AND ({ras})"
def tap(sql,limit,ra,dec,radius):
    receipt={"endpoint":TAP,"adql_sha256":dig(sql.encode()),"limit":limit,
             "status":"QUERY_ERROR","row_cap_reached":None}
    try:
        response=requests.post(TAP,data={"REQUEST":"doQuery","LANG":"ADQL",
                    "FORMAT":"csv","QUERY":sql},headers=HEADERS,timeout=(12,52))
        receipt.update(http=response.status_code,bytes=len(response.content),
                       response_sha256=dig(response.content))
        response.raise_for_status()
        raw=response.content.decode("utf-8-sig",errors="replace")
        if raw.lstrip().startswith("<"):raise ValueError("TAP XML error/unsupported CSV: "+raw[:220])
        rows=[{k.lower().strip():v for k,v in row.items()}
               for row in csv.DictReader(io.StringIO(raw))]
        if not rows and not raw.strip():raise ValueError("empty response")
        if rows and ("ra" not in rows[0] or "dec" not in rows[0]):
            raise ValueError("missing coordinate columns")
        receipt.update(status="OK",returned_rows=len(rows),
                       row_cap_reached=len(rows)>=limit)
        inside=[]
        for r in rows:
            a=fnum(r.get("ra"));b=fnum(r.get("dec"))
            if a is None or b is None:continue
            s=sep(ra,dec,a,b)
            if s<=radius+.2:inside.append((s,r))
        inside.sort(key=lambda q:q[0])
        receipt["within_beam"]=len(inside)
        return receipt,inside
    except Exception as e:
        receipt["error"]=type(e).__name__+": "+str(e)[:280]
        return receipt,[]
def spectroscopy(row):
    ra,dec,radius=float(row.ra_deg),float(row.dec_deg),float(row.beam_fwhm_arcsec)/2
    query=(f"SELECT TOP 200 targetid,mean_fiber_ra AS ra,mean_fiber_dec AS dec,"
           "z,zwarn,deltachi2,spectype,zcat_primary FROM desi_dr1.zpix WHERE "+
           posbox(ra,dec,radius,"mean_fiber_ra","mean_fiber_dec"))
    rec,hits=tap(query,200,ra,dec,radius)
    rv=fnum(row.vopt_kms_or_derived)
    if rv is None:return {"receipt":rec,"redshift_matches":[],"no_radio_redshift":True}
    matching=[]
    for s,r in hits:
        z=fnum(r.get("z"))
        if z is None or z<=-1 or 1+rv/C<=0:continue
        dv=C*math.log((1+z)/(1+rv/C))
        if abs(dv)<=350 and fnum(r.get("zwarn"))==0 and (fnum(r.get("deltachi2")) or 0)>=25 and str(r.get("spectype","")).upper()=="GALAXY":
            matching.append({"targetid":r.get("targetid"),"sep_arcsec":round(s,2),
               "z":z,"delta_v_kms":round(dv,1),"deltachi2":fnum(r.get("deltachi2"))})
    return {"receipt":rec,"redshift_matches":matching[:30],
            "matching_count":len(matching),"no_radio_redshift":False,
            "interpretation":"beam-level velocity association; HI host unproven"}
def photoz(row):
    ra,dec,radius=float(row.ra_deg),float(row.dec_deg),float(row.beam_fwhm_arcsec)/2
    box=posbox(ra,dec,radius,"t.ra","t.dec")
    # The photo_z table is keyed to tractor by RELEASE, BRICKID and OBJID.
    # Query as a spatially restricted join rather than a capped DR10-only list.
    query=(f"SELECT TOP 600 t.ra,t.dec,t.release,t.brickid,t.objid,"
           "t.type,t.flux_r,t.flux_z,t.maskbits,p.z_phot_median,"
           "p.z_phot_l95,p.z_phot_u95 "
           "FROM ls_dr10.tractor AS t JOIN ls_dr10.photo_z AS p "
           "ON t.release=p.release AND t.brickid=p.brickid AND t.objid=p.objid "
           f"WHERE {box}")
    rec,rows=tap(query,600,ra,dec,radius)
    # The fallback is a *catalogue morphology* receipt, never a photo-z.
    fallback=None
    if rec["status"]!="OK":
        alt=(f"SELECT TOP 1000 ra,dec,ls_id,type,flux_r,flux_z,maskbits "
             "FROM ls_dr10.tractor WHERE "+posbox(ra,dec,radius))
        fallback,objects=tap(alt,1000,ra,dec,radius)
        fallback["morphological_nonstellarity_not_redshift"]=True
    else:objects=rows
    zhi=fnum(row.vopt_kms_or_derived)
    zhi=(zhi/C if zhi is not None and zhi>0 else None)
    matches=[]
    if rec["status"]=="OK" and zhi is not None:
        for s,r in rows:
            median,lo,hi=map(fnum,[r.get("z_phot_median"),r.get("z_phot_l95"),r.get("z_phot_u95")])
            if None in (median,lo,hi):continue
            if lo<=zhi<=hi and abs(median-zhi)<=max(.05,.2*(1+zhi)):
                matches.append({"sep_arcsec":round(s,2),"z_hi":round(zhi,6),
                 "z_phot_median":median,"z_phot_l95":lo,"z_phot_u95":hi,
                 "type":r.get("type"),"release":r.get("release"),
                 "brickid":r.get("brickid"),"objid":r.get("objid")})
    return {"receipt":rec,"fallback_tractor":fallback,"photoz_matches":matches[:30],
       "matching_count":len(matches),"observed_tractor_count":len(objects),
       "photoz_identity_proven":False,"photometric_redshift_is_not_HI_host_proof":True}
def hats_open():
    import lsdb
    cats={};err={}
    for mode in ("wide_survey","deep_survey"):
        prefix=HATS_PREFIX+"mer_catalogue_"+mode
        try:
            cats[mode]=lsdb.open_catalog(prefix,columns=[
                "object_id","ra","dec","flux_vis_psf","point_like_prob"])
        except Exception as e:
            err[mode]=type(e).__name__+": "+str(e)[:280]
    return cats,err
def hats_cone(row,cats,open_errors):
    radius=float(row.beam_fwhm_arcsec)/2
    result={"catalogs":{},"complete":False,"nearest":[]}
    for mode in ("wide_survey","deep_survey"):
        if mode not in cats:
            result["catalogs"][mode]={"status":"OPEN_ERROR","error":open_errors.get(mode)}
            continue
        try:
            df=cats[mode].cone_search(float(row.ra_deg),float(row.dec_deg),
                                     radius_arcsec=radius).compute()
            items=[]
            for _,rr in df.iterrows():
                a=fnum(rr.get("ra"));b=fnum(rr.get("dec"))
                if a is None or b is None:continue
                items.append({"mode":mode,"object_id":str(rr.get("object_id")),
                     "sep_arcsec":round(sep(float(row.ra_deg),float(row.dec_deg),a,b),2),
                     "flux_vis_psf":fnum(rr.get("flux_vis_psf")),
                     "point_like_prob":fnum(rr.get("point_like_prob"))})
            items.sort(key=lambda x:x["sep_arcsec"])
            result["catalogs"][mode]={"status":"OK","n_sources_in_beam":len(items),
                                      "nearest":items[:4]}
            result["nearest"].extend(items[:15])
        except Exception as e:
            result["catalogs"][mode]={"status":"QUERY_ERROR",
                    "error":type(e).__name__+": "+str(e)[:280]}
    result["complete"]=all(v.get("status")=="OK" for v in result["catalogs"].values())
    result["nearest"].sort(key=lambda q:q["sep_arcsec"])
    result["nearest"]=result["nearest"][:10]
    return result
def s3key(rec):
    # Reuse PR #211's verified IRSA SIA -> public MER FITS routing.
    link=(rec.get("access_url") or "").strip()
    url=urlsplit(link)
    if url.scheme!="https" or url.hostname not in (
            "irsa.ipac.caltech.edu","irsacloud.ipac.caltech.edu"):
        raise ValueError("non-IRSA SIA URL")
    base="/ibe/data/euclid/dr1/"
    if base not in url.path:raise ValueError("non-DR1 source")
    tail=url.path.split(base,1)[1]
    if not tail.startswith("MER/") or not tail.endswith(".fits"):
        raise ValueError("non-MER FITS science")
    return "dr1/"+tail
def sia(row):
    ra,dec,radius=float(row.ra_deg),float(row.dec_deg),float(row.beam_fwhm_arcsec)/2
    records=[];checks=[]
    for mode in ("wide","deep"):
        params={"COLLECTION":f"euclid_dr1_{mode}_DpdMerBksMosaic",
                "POS":f"CIRCLE {ra:.8f} {dec:.8f} {radius/3600+.003:.7f}",
                "FORMAT":"image/fits","MAXREC":200,"RESPONSEFORMAT":"CSV"}
        rec={"mode":mode,"status":"QUERY_ERROR","params":params}
        try:
            r=requests.get(SIA,params=params,headers=HEADERS,timeout=(12,45))
            rec.update(http=r.status_code,response_sha256=dig(r.content))
            r.raise_for_status()
            dat=r.content.decode("utf-8-sig","replace")
            if dat.lstrip().startswith("<"):raise ValueError("unexpected SIA XML response")
            q=list(csv.DictReader(io.StringIO(dat)))
            rec.update(status="OK",rows=len(q),capped=len(q)>=200)
            for o in q:
                if str(o.get("dataproduct_subtype","")).lower() not in ("","science"):
                    continue
                band=str(o.get("energy_bandpassname","")).upper()
                if band not in BANDS:continue
                try:
                    key=s3key(o)
                    if (band=="VIS" and "/VIS/" not in key) or (band!="VIS" and "/NISP/" not in key):
                        continue
                    records.append((mode,band,key))
                except ValueError:continue
        except Exception as e:rec["error"]=type(e).__name__+": "+str(e)[:280]
        checks.append(rec)
    unique=[]
    for item in records:
        if item not in unique:unique.append(item)
    return checks,unique
def photometric_injections(raw,scale,zero,band,hi_key):
    """Masked, fixed-position exponential injections on actual native MER pixels.

    Measured locally, not a blind all-beam recovery experiment or a validated
    stellar-population upper bound. Deterministic seeded injection positions.
    """
    k=max(1,round(1./scale))
    ny,nx=raw.shape;ny2=ny//k;nx2=nx//k
    if min(ny2,nx2)<48:return {"status":"TOO_SMALL_FOR_SKY_INJECTIONS"}
    a=raw[:ny2*k,:nx2*k]
    z=a.reshape(ny2,k,nx2,k)
    n=np.isfinite(z).sum(axis=(1,3))
    im=np.where(n>=.9*k*k,np.nansum(z,axis=(1,3)),np.nan)
    if np.isfinite(im).mean()<.70:return {"status":"INSUFFICIENT_FINITE_PIXELS"}
    median=float(np.nanmedian(im))
    sky=np.where(np.isfinite(im),im-median,0.)
    pix=im[np.isfinite(im)]
    noise=1.4826*float(np.median(np.abs(pix-np.median(pix))))
    if not np.isfinite(noise) or noise<=0:return {"status":"BAD_SKY_NOISE"}
    smooth=gaussian_filter(sky,sigma=2.)
    sm_med=float(np.median(smooth))
    sm_std=1.4826*float(np.median(np.abs(smooth-sm_med)))
    mask=binary_dilation((smooth-sm_med)>4*max(sm_std,1e-12),
                         iterations=3) | ~np.isfinite(im)
    rng=np.random.default_rng(int(dig((hi_key+band).encode())[:16],16))
    results=[]
    for reff in RE:
        psf=.35 if band=="VIS" else .65
        radius=int(math.ceil(4*reff/(k*scale)+4))
        if 2*radius+2>=min(ny2,nx2):
            results.append({"re_arcsec":reff,"status":"TEMPLATE_TOO_LARGE"});continue
        coords=[]
        for _ in range(300):
            cy=int(rng.integers(radius,ny2-radius))
            cx=int(rng.integers(radius,nx2-radius))
            if mask[cy,cx]:continue
            submask=mask[cy-radius:cy+radius+1,cx-radius:cx+radius+1]
            if submask.mean()>.07:continue
            coords.append((cy,cx))
            if len(coords)==28:break
        if len(coords)<16:
            results.append({"re_arcsec":reff,"status":"INSUFFICIENT_CLEAN_INJECTION_SITES",
                            "sites":len(coords)});continue
        yy,xx=np.mgrid[-radius:radius+1,-radius:radius+1]
        rr=np.hypot(xx,yy)*k*scale
        template=np.exp(-1.678*rr/reff)
        template=gaussian_filter(template,sigma=max(.2,psf/(k*scale)))
        template/=template.sum()
        weight=template/(template*template).sum()
        scores=np.array([np.sum(sky[y-radius:y+radius+1,x-radius:x+radius+1]*weight)
                         for y,x in coords])
        scatter=1.4826*np.median(np.abs(scores-np.median(scores)))
        if not np.isfinite(scatter) or scatter<=0:
            results.append({"re_arcsec":reff,"status":"NO_SKY_SCATTER"});continue
        threshold=max(5*scatter,float(np.quantile(scores,.99)))
        # Actual artificial light is added to the observed sky cutouts and
        # recovered by the same template score, not an idealized noise formula.
        flux_grid=np.geomspace(max(scatter*.4,1e-7),max(scatter*.4,1e-7)*300,14)
        rates=[]
        for flux in flux_grid:
            detected=0
            for y,x in coords:
                injected=sky[y-radius:y+radius+1,x-radius:x+radius+1]+flux*template
                snr=float(np.sum(injected*weight))
                if snr>threshold:detected+=1
            rates.append(detected/len(coords))
        idx=next((i for i,v in enumerate(rates) if v>=.90),None)
        if idx is None:
            results.append({"re_arcsec":reff,"status":"INJECTION_90PCT_NOT_REACHED",
                            "sites":len(coords)});continue
        limit=float(flux_grid[idx])
        m=zero-2.5*math.log10(limit)
        results.append({"re_arcsec":reff,"status":"LOCAL_INJECTION_90PCT",
            "sites":len(coords),"recovered_fraction":rates[idx],
            "false_positive_threshold_native_flux":float(threshold),
            "injected_total_flux_native":limit,"total_mag_AB":round(m,3),
            "mean_SB_within_re_mag_AB_arcsec2":round(m+2.5*math.log10(2*math.pi*reff**2),3),
            "source_profile":"exponential_n1_PSF_approx","position_known_injection":True,
            "blank_selection":"masked_local_sky_only"})
    return {"status":"LOCAL_INJECTION_EXPERIMENT","binned_scale_arcsec":round(k*scale,3),
        "masked_fraction":round(float(mask.mean()),4),"noise_per_bin":float(noise),
        "injections":results,"science_quality_flags_verified":False,
        "blind_beam_recovery_calibrated":False}
def read_stamp(key,row,band):
    ra,dec,radius=float(row.ra_deg),float(row.dec_deg),float(row.beam_fwhm_arcsec)/2
    record={"band":band,"s3_key":key,"status":"PIXEL_ERROR"}
    uri=f"s3://{BUCKET}/{key}"
    try:
        with fits.open(uri,use_fsspec=True,mode="readonly",
          fsspec_kwargs={"anon":True,"client_kwargs":{"region_name":"us-east-1"},
                         "config_kwargs":{"connect_timeout":10,"read_timeout":55,
                                          "retries":{"max_attempts":2}}}) as hdul:
            candidates=[h for h in hdul if h.header.get("NAXIS")==2 and
                        h.header.get("NAXIS1",0)>0 and h.header.get("NAXIS2",0)>0]
            if not candidates:raise ValueError("missing image HDU")
            h=candidates[0];hdr=h.header;wcs=WCS(hdr).celestial
            x,y=wcs.world_to_pixel(SkyCoord(ra*u.deg,dec*u.deg))
            scale=float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))*3600)
            if not (.07<scale<.13):raise ValueError("invalid MER pixel scale")
            if not(np.isfinite(x) and np.isfinite(y)):raise ValueError("missing valid WCS")
            pad=int(math.ceil(radius/scale))+3
            nx,ny=int(hdr["NAXIS1"]),int(hdr["NAXIS2"])
            x0=max(0,int(round(x))-pad);x1=min(nx,int(round(x))+pad+1)
            y0=max(0,int(round(y))-pad);y1=min(ny,int(round(y))+pad+1)
            if min(x1-x0,y1-y0)<30:raise ValueError("target outside science image")
            if (x1-x0)*(y1-y0)>7000000:raise ValueError("cutout exceeds pixel budget")
            pixels=np.asarray(h.section[y0:y1,x0:x1],dtype=np.float32)
            yy,xx=np.indices(pixels.shape,dtype=np.float32)
            inside=(xx+x0-x)**2+(yy+y0-y)**2<=(radius/scale)**2
            finite=float(np.isfinite(pixels[inside]).mean()) if inside.any() else 0.
            full=(x-pad>=0 and y-pad>=0 and x+pad<nx and y+pad<ny)
            source_unit=str(hdr.get("BUNIT","")).strip()
            zero=fnum(hdr.get("MAGZERO"))
            record.update(native_sha256=dig(np.ascontiguousarray(pixels,dtype="<f4").tobytes()),
                read_shape=[int(pixels.shape[0]),int(pixels.shape[1])],
                native_scale_arcsec=round(scale,4),magzero=zero,
                bunit=source_unit,beam_finite_fraction=round(finite,4),
                full_beam_in_single_mosaic=bool(full),
                image_data_release=str(hdr.get("DATASETR",""))[:60],
                pplanid=str(hdr.get("PPLANID",""))[:120])
            if zero is None or not source_unit.lower().endswith("/s"):
                record["status"]="UNVALIDATED_PHOTOMETRIC_CALIBRATION"
            elif not full or finite<.97:
                record["status"]="PARTIAL_BEAM_SCIENCE_PIXELS"
            else:
                # All four bands use the FITS MAGZERO and native flux units;
                # NISP 'ELECTRON/s' is not excluded as a spurious unit.
                record["status"]="FULL_BEAM_NATIVE_PIXELS"
                record["injections"]=photometric_injections(
                    pixels,scale,zero,band,str(row.survey)+":"+str(row.id))
    except Exception as e:
        record["error"]=type(e).__name__+": "+str(e)[:360]
    return record
def imaging(row):
    receipts,entries=sia(row)
    bands={}
    for band in BANDS:
        options=[(mode,key) for mode,b,key in entries if b==band]
        trials=[]
        for mode,key in options[:3]:
            r=read_stamp(key,row,band)
            r["mode"]=mode
            trials.append(r)
            if r["status"]=="FULL_BEAM_NATIVE_PIXELS":break
        if trials:
            full=next((r for r in trials if r["status"]=="FULL_BEAM_NATIVE_PIXELS"),None)
            bands[band]=full or max(trials,key=lambda r:r.get("beam_finite_fraction",0))
            bands[band]["attempted_tiles"]=len(trials)
        else:
            bands[band]={"band":band,"status":"NO_SIA_MER_COVERAGE",
                         "attempted_tiles":0}
    return {"sia_receipts":receipts,"bands":bands,
        "full_four_band_science_pixels":all(
            bands[b]["status"]=="FULL_BEAM_NATIVE_PIXELS" for b in BANDS)}
def conditional_mass(row,images,spec,photo):
    """Explicitly conditional old/dusty H model, never a validated Mstar limit."""
    h=images["bands"]["H"]
    if not images["full_four_band_science_pixels"]:return None
    if spec.get("matching_count") or photo.get("matching_count"):return None
    d=fnum(row.distance_mpc);loghi=fnum(row.logMHI)
    if d is None or d<=0 or loghi is None:return None
    experiments=h.get("injections",{}).get("injections",[])
    measurements=[q for q in experiments if q.get("status")=="LOCAL_INJECTION_90PCT"
                  and q.get("re_arcsec")==10.]
    if not measurements:return None
    mag=measurements[0]["total_mag_AB"]
    # Solar H AB~4.71; maximally conservative example M/L_H=5, A_H=1 mag,
    # absolute K-correction allowance=0.5 mag. Not a population/SED posterior.
    distance_modulus=5*math.log10(d*1e6/10)
    upper=5 * 10**(.4*(4.71-(mag-distance_modulus)+1.5))
    return {"status":"MODEL_CONDITIONAL_NOT_VALIDATED","assumed_distance_mpc":d,
            "assumed_ML_H":5,"assumed_A_H_mag":1.,"assumed_K_allowance_mag":.5,
            "assumed_solar_H_AB_mag":4.71,"assumed_re_arcsec":10.,
            "illustrative_Mstar_upper_msun":round(upper,2),
            "illustrative_MHI_over_Mstar_lower":round(10**loghi/upper,4),
            "valid_physical_upper_bound":False}
def examine(row,cats,errors):
    rec={"survey":str(row.survey),"id":str(row.id),
         "source_row":int(row.source_row),"ra_deg":float(row.ra_deg),
         "dec_deg":float(row.dec_deg),"beam_fwhm_arcsec":float(row.beam_fwhm_arcsec),
         "velocity_kms":fnum(row.vopt_kms_or_derived),
         "w50_kms":fnum(row.w50_kms),"logMHI":fnum(row.logMHI),
         "distance_mpc":fnum(row.distance_mpc),"coverage_code":int(row.euclid_coverage_code),
         "nearby_HI_20arcmin_dv500":int(row.nearby_HI_20arcmin_dv500),
         "validated_starless":False,"validated_Mstar_upper_msun":None,
         "validated_MHI_Mstar_lower":None}
    stages=[("desi_dr1",lambda:spectroscopy(row)),
        ("legacy_dr10_photoz",lambda:photoz(row)),
        ("euclid_hats",lambda:hats_cone(row,cats,errors)),
        ("euclid_mer_beam_pixels",lambda:imaging(row))]
    for name,fun in stages:
        try:rec[name]=fun()
        except Exception as e:
            rec[name]={"status":"UNEXPECTED_STAGE_ERROR","error":repr(e)[:500]}
    try:
        rec["conditional_stellar_mass"]=conditional_mass(
             row,rec["euclid_mer_beam_pixels"],rec["desi_dr1"],
             rec["legacy_dr10_photoz"])
    except Exception as e:
        rec["conditional_stellar_mass"]={"status":"MODEL_ERROR","error":str(e)[:250]}
    rec["hi_checks"]={
        "galactic_hvc_risk":(rec["velocity_kms"] is None or abs(rec["velocity_kms"])<600),
        "rfi_spectral_line_cube_review":"NOT_PERFORMED",
        "interferometric_position_confirmation":"NOT_PERFORMED",
        "tidal_debris_group_membership":"NOT_EXCLUDED",
        "beam_confusion":("NEARBY_HI_POSSIBLE" if rec["nearby_HI_20arcmin_dv500"] else "NOT_EXCLUDED"),
        "w50_not_dynamical_mass":True,
        "previous_dark_candidates_exact_AGC":str(row.id) in KNOWN_AGC if str(row.survey)=="ALFALFA_alpha100" else None,
        "literature_alias_reconciliation":"NOT_EXHAUSTIVE"}
    rec["referee_status"]="UNVALIDATED: HI cube/RFI/tidal/confusion and image-quality/SED checks remain"
    return rec
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--batch-size",type=int,default=8)
    args=ap.parse_args()
    if args.batch_size<1 or args.batch_size>32:raise ValueError("batch size must be 1..32")
    path=OUT/"footprint_hits.csv"
    if not path.exists():raise RuntimeError("Missing frozen all-1375 position source table")
    allrows=pd.read_csv(path,dtype={"id":str})
    allrows=allrows.sort_values(["survey","source_row"],kind="stable").reset_index(drop=True)
    BATCH.mkdir(parents=True,exist_ok=True)
    inp_hash=dig(path.read_bytes())
    expected=(len(allrows)+args.batch_size-1)//args.batch_size
    existing=sorted(BATCH.glob("batch_????_manifest.json"))
    done=set()
    for p in existing:
        meta=json.loads(p.read_text())
        if meta["input_sha256"]!=inp_hash or meta["batch_size"]!=args.batch_size:
            raise RuntimeError("Batch provenance/input mismatch: "+str(p))
        data=BATCH/p.name.replace("_manifest.json",".jsonl.gz")
        if not data.exists() or dig(data.read_bytes())!=meta["data_sha256"]:
            raise RuntimeError("Invalid/missing batch payload: "+str(data))
        done.add(meta["batch_index"])
    index=next((i for i in range(expected) if i not in done),None)
    if index is None:
        print("ALL_COMPLETE",len(allrows),expected,flush=True)
        return
    rows=allrows.iloc[index*args.batch_size:(index+1)*args.batch_size]
    print("BATCH",index,"/",expected,"rows",len(rows),flush=True)
    cats,errors=hats_open()
    output=[]
    for row in rows.itertuples(index=False):
        start=time.time()
        result=examine(row,cats,errors)
        result["elapsed_seconds"]=round(time.time()-start,2)
        output.append(result)
        print("POSITION",len(output),len(rows),result["survey"],result["id"],
              {b:d["status"] for b,d in result.get("euclid_mer_beam_pixels",{}).get("bands",{}).items()},
              "elapsed",result["elapsed_seconds"],flush=True)
    raw="\n".join(json.dumps(x,separators=(",",":"),allow_nan=False) for x in output)+"\n"
    file=BATCH/f"batch_{index:04d}.jsonl.gz"
    # Compress deterministically for manifest SHA verification on reruns.
    with file.open("wb") as f:
        with gzip.GzipFile(filename="",mode="wb",fileobj=f,mtime=0) as z:z.write(raw.encode())
    jsonwrite(BATCH/f"batch_{index:04d}_manifest.json",{
       "batch_index":index,"batch_size":args.batch_size,
       "input_sha256":inp_hash,"data_sha256":dig(file.read_bytes()),
       "n_positions":len(output),"first_row":index*args.batch_size,
       "last_row_exclusive":index*args.batch_size+len(rows),
       "full_four_band_pixel_beams":sum(x.get("euclid_mer_beam_pixels",{}).get(
             "full_four_band_science_pixels",False) for x in output),
       "validated_starless":0})
    jsonwrite(OUT/"beam_followup_status.json",{
       "input_sha256":inp_hash,"input_positions":len(allrows),
       "completed_batches":len(done)+1,"total_batches":expected,
       "completed_positions":min(len(allrows),(len(done)+1)*args.batch_size),
       "complete":len(done)+1==expected,
       "status":"OBSERVATIONAL_FOLLOWUP_NOT_STARLESS_CERTIFICATION"})
if __name__=="__main__":main()
