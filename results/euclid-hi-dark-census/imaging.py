#!/usr/bin/env python3
"""Bounded native Euclid DR1 VIS/NISP pixel sensitivity audit on frozen HI positions.

No source is called starless. The reported aperture sensitivities (when present)
are descriptive/local and not a completeness-calibrated stellar mass census.
"""
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astroquery.ipac.irsa import Irsa
import s3fs

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
S3=s3fs.S3FileSystem(anon=True,config_kwargs={"connect_timeout":8,"read_timeout":20})

def utc():return dt.datetime.now(dt.timezone.utc).isoformat()
def sha(x):return hashlib.sha256(x).hexdigest()
def finite(x):
    try:
        v=float(x);return v if math.isfinite(v) else None
    except (ValueError,TypeError):return None
def sky_limit(pixels,zp,ap_radius_arcsec=3.0,scale=0.1):
    """Empirical blank-aperture scatter for a fixed top-hat model; no PSF correction."""
    n,m=pixels.shape
    yy,xx=np.indices((n,m),dtype=float)
    x0=(m-1)/2;y0=(n-1)/2
    rad=ap_radius_arcsec/scale
    center=(xx-x0)**2+(yy-y0)**2<=rad**2
    finite_mask=np.isfinite(pixels)
    # No calibrated source/bright-star mask is present: be conservative and
    # mark every result preliminary; blank sky apertures can contain galaxies.
    sky=finite_mask&(((xx-x0)**2+(yy-y0)**2)>=(11/scale)**2)&(
        ((xx-x0)**2+(yy-y0)**2)<=(20/scale)**2)
    if np.mean(finite_mask)<.97 or np.sum(sky)<3000 or not np.all(finite_mask[center]):
        return {"status":"INVALID_PIXEL_COVERAGE","finite_fraction":float(np.mean(finite_mask))}
    baseline=float(np.median(pixels[sky]))
    target=float(np.sum(pixels[center]-baseline))
    samples=[]
    for deg in range(0,360,30):
        a=math.radians(deg)
        dx=(15/scale)*math.cos(a);dy=(15/scale)*math.sin(a)
        test=(xx-(x0+dx))**2+(yy-(y0+dy))**2<=rad**2
        if np.all(finite_mask[test]):samples.append(float(np.sum(pixels[test]-baseline)))
    if len(samples)<8:return {"status":"INSUFFICIENT_BACKGROUND_APERTURES","n":len(samples)}
    med=float(np.median(samples))
    scatter=1.4826*float(np.median(np.abs(np.array(samples)-med)))
    if not (math.isfinite(scatter) and scatter>0):
        return {"status":"UNMEASURABLE_EMPIRICAL_NOISE","scatter":scatter}
    # A range of correlated apertures and unsubtracted contaminants precludes
    # treating this as a validated 5-sigma completeness threshold.
    sensitivity=5*scatter
    aperture_area=np.count_nonzero(center)*scale**2
    avg_sb=zp-2.5*math.log10(sensitivity/aperture_area) if sensitivity>0 else None
    flux_limit_mag=zp-2.5*math.log10(sensitivity) if sensitivity>0 else None
    return {"status":"EXPLORATORY_NOT_VALIDATED",
      "center_flux_adu_per_sec":target,"blank_median_adu_per_sec":med,
      "blank_mad_sigma_adu_per_sec":scatter,"central_SNR_rough":(target-med)/scatter,
      "aperture_radius_arcsec":ap_radius_arcsec,"aperture_area_arcsec2":aperture_area,
      "n_blank_apertures":len(samples),"empirical_5sigma_mag_AB":flux_limit_mag,
      "empirical_5sigma_mean_SB_mag_AB_arcsec2":avg_sb,
      "flag_mask_verified":False,"completeness_calibrated":False,
      "interpretation":"Local background-aperture diagnostic only; unmasked sources/correlated samples."}

def scan_one(hi):
    coord=SkyCoord(ra=float(hi["ra_deg"])*u.deg,dec=float(hi["dec_deg"])*u.deg,frame="icrs")
    entry={"id":hi["hi_id"],"ra_deg":hi["ra_deg"],"dec_deg":hi["dec_deg"],
      "timestamp":utc(),"image_sensitivity":{}, "image_access":[],"validated_Mstar_limit":None,
      "validated_MHI_over_Mstar_bound":None}
    for mode in ("deep","wide"):
        collection=f"euclid_dr1_{mode}_DpdMerBksMosaic"
        try:
            images=Irsa.query_sia(pos=(coord,7*u.arcsec),collection=collection)
            entry["image_access"].append({"collection":collection,"rows":len(images)})
            img=images[(images["dataproduct_subtype"]=="science") &
                       (images["facility_name"]=="Euclid")]
            if not len(img):continue
            for band in ("VIS","Y","J","H"):
                found=img[img["energy_bandpassname"]==band]
                if not len(found):continue
                for row in found:
                    try:
                        cloud=json.loads(str(row["cloud_access"]))
                        key=cloud["aws"]["key"];bucket=cloud["aws"]["bucket_name"]
                        fpath=f"{bucket}/{key}"
                        # Avoid extra deep single-visit wide mosaics.
                        if not S3.exists(fpath):continue
                        with fits.open("s3://"+fpath,fsspec_kwargs={"anon":True},do_not_scale_image_data=False) as hdul:
                            head=hdul[0].header
                            if mode=="wide" and "-DATA_RELEASE-" not in str(head.get("PPLANID","")):
                                continue
                            unit=str(head.get("BUNIT","")).strip()
                            if not (unit=="ADU/s" or
                                    (band!="VIS" and unit.lower() in ("electrons","electron","e-"))):
                                raise ValueError("Uncalibrated/nonstandard BUNIT: "+unit)
                            zp=finite(head.get("MAGZERO"))
                            if zp is None:raise ValueError("Missing AB zero point")
                            wcs=WCS(head)
                            subset=Cutout2D(hdul[0].section,position=coord,
                                            size=44*u.arcsec,wcs=wcs,mode="partial",fill_value=np.nan)
                            pix=np.asarray(subset.data,dtype=float)
                            # Native MER pixel scale 0.1"/pix, validate rather than assume.
                            cd=abs(float(wcs.proj_plane_pixel_scales()[0].to_value(u.deg)))*3600
                            if not .095<cd<.105:raise ValueError("Unexpected pixel scale")
                            lim=sky_limit(pix,zp,scale=cd)
                            lim.update({"tile_id":str(row["obs_id"])[:9],
                                "s3_key":key,"s3_bucket":bucket,"source_MAGZERO_AB":zp,
                                "source_BUNIT":str(head.get("BUNIT")),
                                "source_DATASER":str(head.get("DATASETR","")),
                                "source_PPLANID":str(head.get("PPLANID",""))[:190]})
                            entry["image_sensitivity"][mode+"_"+band]=lim
                        break
                    except Exception as e:
                        entry["image_access"].append({"collection":collection,"band":band,
                            "status":"IMAGE_ERROR","error":str(e)[:200]})
        except Exception as e:
            entry["image_access"].append({"collection":collection,"status":"SIA_ERROR",
                                         "error":str(e)[:280]})
    entry["valid_local_pixel_estimates"]=sum(x.get("status")=="EXPLORATORY_NOT_VALIDATED"
                                      for x in entry["image_sensitivity"].values())
    return entry
def main():
    data=json.loads((OUT/"followup.json").read_text())
    options=[x for x in data["results"] if x.get("disposition") not in (
        "SPECTROSCOPIC_GALAXY_REDSHIFT_ASSOCIATION","PROCESSING_ERROR") and
        x.get("footprint")=="CENTER_IN_MOC"]
    # Strictly frozen case selection: first two by survey/source ID,
    # no post-photometry ranking or peeking.
    options.sort(key=lambda x:(x["survey"],x["hi_id"]))
    options=options[:2]
    payload={"utc":utc(),"policy":"First two non-spectroscopic source positions sorted survey/id from deterministic catalog-referee sample",
             "tested":len(options),"per_position":[]}
    for o in options:
        try:
            rec=scan_one(o)
        except Exception as exc:
            rec={"id":o["hi_id"],"error":str(exc)[:300],
                 "validated_Mstar_limit":None,"validated_MHI_over_Mstar_bound":None}
        payload["per_position"].append(rec)
        print("Euclid pixels:",rec["id"],rec.get("valid_local_pixel_estimates",0),flush=True)
    payload["validated_stellar_upper_limits_count"]=0
    payload["reason"]="Not validated: background apertures lack exposure/segmentation-quality masks, injection-completeness calibration, and robust source SED/distance association. Even a numeric local limit is NOT an astrophysical Mstar bound."
    (OUT/"native_pixel_checks.json").write_text(json.dumps(payload,indent=2,allow_nan=False)+"\n")
if __name__=="__main__":main()
