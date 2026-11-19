#!/usr/bin/env python3
"""Position-forced, failure-aware audit of the 35 unassociated top-60 eRASS1 sources.

The prior catalogue's *lack of a secure association* is NEVER substituted for
image coverage, a measured upper limit, a spectral fit, or archival absence.
"""
from __future__ import annotations
import csv, hashlib, io, json, math, os, re, tarfile, tempfile, time, threading
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode, quote
import numpy as np
import requests
from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
PREV=ROOT.parent/"erass1-xdins-search"/"output"
OUT.mkdir(exist_ok=True)
BASE="https://www.legacysurvey.org/viewer/fits-cutout"
KNOWN={
 "110607500001":"PSR J0659+1414 (B0656+14; known rotation-powered INS)",
 "111112000001":"RX J0720.4-3125",
 "112113200001":"RX J0806.4-4123",
 "119606900001":"RX J1308.6+2127 (RBS 1223)",
 "128212900001":"RX J1856.5-3754",
 "116414100004":"PSR J1057-5226 (B1055-52; known rotation-powered INS)",
}
PS1="https://ps1images.stsci.edu/cgi-bin/"
NU_G=6.30e14
AB_NMGY_FNU=3.631e-29 # erg cm^-2 s^-1 Hz^-1; AB 22.5 = 1 nMgy
HEAD={"User-Agent":"Open-Space-XDINS-35-referee/1.0 (public research reproducibility)"}
LS_LOCK=threading.Lock()
LS_NEXT=[0.0]

def readcsv(p):
    with open(p,newline="") as f: return list(csv.DictReader(f))
def writejson(name,obj):
    (OUT/name).write_text(json.dumps(obj,indent=2,sort_keys=True,allow_nan=False)+"\n")
def sha(p):
    h=hashlib.sha256()
    with open(p,"rb") as f:
        for b in iter(lambda:f.read(1<<20),b""):h.update(b)
    return h.hexdigest()
def finite(v):
    try:
        a=float(v)
        return a if math.isfinite(a) else None
    except (ValueError,TypeError):return None
def fetch(s,url,params=None,limit=12_000_000,timeout=(10,28)):
    r=s.get(url,params=params,timeout=timeout,stream=True)
    uri=r.url
    r.raise_for_status()
    buf=io.BytesIO()
    for b in r.iter_content(131072):
        buf.write(b)
        if buf.tell()>limit:raise RuntimeError("response_exceeds_size_cap_"+str(limit))
    b=buf.getvalue()
    if not b:raise RuntimeError("empty_response")
    return b,{"url":uri,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),"http":r.status_code}
def post_fetch(s,url,form,limit=4_000_000,timeout=(12,65)):
    r=s.post(url,data=form,timeout=timeout,stream=True)
    r.raise_for_status()
    buf=io.BytesIO()
    for b in r.iter_content(131072):
        buf.write(b)
        if buf.tell()>limit:raise RuntimeError("response_exceeds_size_cap")
    b=buf.getvalue()
    if not b:raise RuntimeError("empty_response")
    return b,{"url":r.url,"bytes":len(b),"sha256":hashlib.sha256(b).hexdigest(),"http":r.status_code,
              "request_method":"POST","ADQL_sha256":hashlib.sha256(form["QUERY"].encode()).hexdigest()}
def img_measure(img, header, ra, dec, pixscale, band, survey):
    img=np.asarray(img,dtype=float)
    if img.ndim!=2 or min(img.shape)<25:raise ValueError("not_2D_or_too_small")
    w=WCS(header,relax=True)
    if w.has_celestial:
        x,y=w.celestial.world_to_pixel_values(ra,dec)
    else:
        x,y=(img.shape[1]-1)/2,(img.shape[0]-1)/2
    if not np.isfinite(x+y) or x<12 or y<12 or x>img.shape[1]-13 or y>img.shape[0]-13:
        raise ValueError("target_off_image")
    yy,xx=np.indices(img.shape)
    d=np.hypot(xx-x,yy-y)
    radius=max(2.0,1.25/pixscale) if survey!="unWISE" else max(2.0,5.0/pixscale)
    radius=min(radius,6)
    aperture=(d<=radius)&np.isfinite(img)
    ring=(d>radius*2)&(d<min(11.5,radius*4))&np.isfinite(img)
    if aperture.sum()<10 or ring.sum()<30:raise ValueError("insufficient_valid_pixels")
    _,bg,sd=sigma_clipped_stats(img[ring],sigma=2.8,maxiters=4)
    flux=float(np.sum(img[aperture]-bg))
    # Correlated-noise and source-confusion conservative empirical aperture scatter.
    off=[]
    rr=radius
    for cy in np.arange(8,img.shape[0]-8,4.0):
        for cx in np.arange(8,img.shape[1]-8,4.0):
            if np.hypot(cx-x,cy-y)<3.2*rr:continue
            ap=np.hypot(xx-cx,yy-cy)<=rr
            if ap.sum()>0.9*aperture.sum() and np.all(np.isfinite(img[ap])):
                off.append(float(np.sum(img[ap]-bg)))
    empirical=None
    if len(off)>=12:
        z=np.asarray(off)
        empirical=float(1.4826*np.median(np.abs(z-np.median(z))))
    simple=float(sd*np.sqrt(aperture.sum()))
    sig=max(simple,empirical or 0)
    if not np.isfinite(sig) or sig<=0:raise ValueError("nonpositive_aperture_noise")
    # Aperture correction 1/0.7 makes a conservative FLUX upper limit.
    limit_nmg=5*sig/0.7
    snr=flux/sig
    result={"band":band,"survey":survey,"forced_nmg_ab":round(flux,5),
            "sigma_aperture_nmg_ab":round(sig,5),"five_sigma_nmg_ab":round(limit_nmg,5),
            "snr":round(float(snr),3),"position_pixel":[round(float(x),2),round(float(y),2)],
            "pixscale_arcsec":pixscale,"aperture_radius_arcsec":round(radius*pixscale,3),
            "offsource_aperture_samples":len(off),"detected_at_position_5sigma":bool(snr>=5),
            "crowding_warning":bool(empirical is not None and empirical>2*simple),
            "status":"measured"}
    if limit_nmg>0:result["five_sigma_ab_mag"]=round(22.5-2.5*math.log10(limit_nmg),3)
    return result
def apply_units(band_data,header,survey,band):
    img=np.asarray(band_data,dtype=float)
    if survey=="PS1":
        exptime=finite(header.get("EXPTIME"))
        if not exptime or exptime<=0:raise ValueError("PS1_EXPTIME_missing_no_calibration")
        zp=25+2.5*math.log10(exptime)
        img*=10**((22.5-zp)/2.5)
    elif survey=="unWISE":
        img*=10**(-0.4*({"1":2.699,"2":3.339}[band]))
    # LS10 viewer image coadds are nominally calibrated in AB nanomaggies.
    return img
def phot_fits(b,ra,dec,survey,band,pixscale):
    with fits.open(io.BytesIO(b),memmap=False,ignore_missing_simple=True) as hd:
        hdu=next((h for h in hd if isinstance(h,(fits.PrimaryHDU,fits.ImageHDU,fits.CompImageHDU))
                 and h.data is not None and np.ndim(h.data)>=2),None)
        if hdu is None:raise ValueError("no_image_hdu")
        d=np.asarray(hdu.data)
        hdr=hdu.header
        if d.ndim==3:
            if survey=="LS10":
                order=str(hdr.get("BANDS","grz")).strip()
                index=order.index(band) if band in order else "grz".index(band)
            else:index=0
            d=d[index]
        if d.ndim>3:raise ValueError("4D_image_unexpected")
        return img_measure(apply_units(d,hdr,survey,band),hdr,ra,dec,pixscale,band,survey)
def attempt(name,func):
    try:return func()
    except Exception as e:return {"status":"unavailable_or_uninterpretable","error":type(e).__name__+": "+str(e)[:220]}
def legacy(s,row):
    ra,dec=float(row["ra"]),float(row["dec"])
    # A global host pace avoids turning rate limits into false sky non-coverage.
    for trial in range(3):
        with LS_LOCK:
            time.sleep(max(0,LS_NEXT[0]-time.monotonic()))
            LS_NEXT[0]=time.monotonic()+5.0
            try:
                b,receipt=fetch(s,BASE,{"ra":ra,"dec":dec,"layer":"ls-dr10","size":80,"pixscale":0.262,"bands":"grz"},limit=8_000_000)
                break
            except requests.HTTPError as e:
                if e.response is None or e.response.status_code!=429:raise
                LS_NEXT[0]=time.monotonic()+18*(trial+1)
    else:
        raise RuntimeError("Legacy viewer still HTTP_429 after three spaced requests")
    ans={"receipt":receipt,"bands":{},"status":"downloaded"}
    for band in "grz":
        ans["bands"][band]=attempt(band,lambda band=band: phot_fits(b,ra,dec,"LS10",band,0.262))
    if not any(v.get("status")=="measured" for v in ans["bands"].values()):ans["status"]="uninterpretable"
    return ans
def ps1(s,row):
    ra,dec=float(row["ra"]),float(row["dec"])
    if dec<-30:return {"status":"outside_declination_footprint"}
    data,receipt=fetch(s,PS1+"ps1filenames.py",{"ra":ra,"dec":dec,"filters":"gr","type":"stack"},limit=500000)
    table=Table.read(io.BytesIO(data),format="ascii")
    cols={c.lower():c for c in table.colnames}
    if "filename" not in cols or "filter" not in cols:raise ValueError("PS1_missing_table_columns")
    out={"listing":receipt,"bands":{},"status":"queried_stacks"}
    for band in "gr":
        found=[t for t in table if str(t[cols["filter"]]).strip()==band]
        if not found:
            out["bands"][band]={"status":"no_stack_listed"}
            continue
        filename=str(found[0][cols["filename"]]).strip()
        def grab():
            b,r=fetch(s,PS1+"fitscut.cgi",{"ra":ra,"dec":dec,"size":96,"format":"fits","red":filename},limit=8_000_000)
            return {"receipt":r,**phot_fits(b,ra,dec,"PS1",band,0.25)}
        out["bands"][band]=attempt(band,grab)
    return out
def unwise(s,row):
    ra,dec=float(row["ra"]),float(row["dec"])
    # native unWISE coadds, rather than CatWISE catalogue nonmatches
    b,receipt=fetch(s,"https://unwise.me/cutout_fits",
       {"version":"neo7","ra":ra,"dec":dec,"size":64,"bands":"12",
        "file_img_m":"on","file_invvar_m":"on"},limit=20_000_000,timeout=(10,35))
    out={"receipt":receipt,"status":"downloaded","bands":{}}
    with tarfile.open(fileobj=io.BytesIO(b),mode="r:gz") as t:
        fs=[m for m in t.getmembers() if m.isfile() and m.name.endswith(".fits") and "-img-m" in m.name]
        for band in "12":
            matches=[m for m in fs if re.search(r"-w"+band+r"[-.]",m.name)]
            if not matches:
                out["bands"][band]={"status":"no_band_image_in_cutout"}
                continue
            def grab():
                f=t.extractfile(matches[0])
                return {"tile_file":matches[0].name,**phot_fits(f.read(),ra,dec,"unWISE",band,2.75)}
            out["bands"][band]=attempt(band,grab)
    return out
def des(s,row):
    ra,dec=float(row["ra"]),float(row["dec"])
    # SIA returns image locations, not a catalogue association.
    b,r=fetch(s,"https://datalab.noirlab.edu/sia/des_dr2",
              {"POS":f"{ra},{dec}","SIZE":"0.025"},limit=3_000_000,timeout=(8,22))
    tb=Table.read(io.BytesIO(b),format="votable")
    cols={c.lower():c for c in tb.colnames}
    urlcol=next((cols[k] for k in ("access_url","accessurl","access.reference") if k in cols),None)
    bandcol=next((cols[k] for k in ("bandpass","filter","band") if k in cols),None)
    out={"sia_receipt":r,"sia_rows":len(tb),"status":"queried_sia","bands":{}}
    if not urlcol:return dict(out,status="metadata_only_url_column_unavailable")
    for band in "gr":
        found=[x for x in tb if (not bandcol or band in str(x[bandcol]).lower())]
        if not found:
            out["bands"][band]={"status":"no_sia_image_found"}
            continue
        def grab():
            url=str(found[0][urlcol]).strip()
            if not url.startswith("https://"):raise ValueError("unsafe_or_invalid_sia_url")
            data,rec=fetch(s,url,limit=9_000_000,timeout=(8,30))
            # DES SIA cutouts may be survey-calibrated, but units must be proven.
            with fits.open(io.BytesIO(data)) as hd:
                hs=[h for h in hd if h.data is not None and getattr(h.data,"ndim",0)==2]
                if not hs:raise ValueError("SIA_no_2D_image")
                h=hs[0].header
                bunit=str(h.get("BUNIT","")).lower()
                if "nanomagg" not in bunit:raise ValueError("DES_flux_units_not_verified:"+str(h.get("BUNIT")))
            return {"receipt":rec,**phot_fits(data,ra,dec,"LS10",band,0.263)}
        out["bands"][band]=attempt(band,grab)
    return out
def dust(s,row):
    ra,dec=float(row["ra"]),float(row["dec"])
    b,rec=fetch(s,"https://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust",
                {"locstr":f"{ra} {dec}"},limit=100000,timeout=(8,18))
    t=b.decode(errors="replace")
    candidates=[r"<refPixelValueSandF>\\s*([\\d.eE+\\-]+)",
                r"<refPixelValueSFD>\\s*([\\d.eE+\\-]+)",
                r"<meanValueSandF>\\s*([\\d.eE+\\-]+)",
                r"<meanValueSFD>\\s*([\\d.eE+\\-]+)"]
    vals=[float(m.group(1)) for p in candidates if (m:=re.search(p,t))]
    if not vals:raise ValueError("IRSA_dust_XML_no_SFD_value; excerpt: "+t[:240].replace("\n"," "))
    ebv=vals[0]
    return {"status":"SandF_integrated_dust_proxy","E_BV_mag":round(ebv,4),
            "NH_proxy_cm2":round(8.8e21*ebv,0),"receipt":rec,
            "warning":"SFD full-column dust not a measured X-ray N_H; dust-to-gas varies"}
def one(row):
    s=requests.Session();s.headers.update(HEAD)
    out={"uid":row["uid"],"iauname":row["iauname"],"ra":float(row["ra"]),"dec":float(row["dec"])}
    coord=SkyCoord(out["ra"]*u.deg,out["dec"]*u.deg,frame="icrs")
    out["galactic_l_b_deg"]=[round(float(coord.galactic.l.deg),3),round(float(coord.galactic.b.deg),3)]
    for survey,fn in [("LS10",legacy),("PS1",ps1),("DES_DR2",des),("unWISE",unwise),("galactic_dust",dust)]:
        out[survey]=attempt(survey,lambda fn=fn:fn(s,row))
    s.close()
    return out
def heasarc_archives(rows):
    # A catalogue association is an archival detection, not a light curve or
    # a proof of long-term non-variability. Never interpret a failed query as absent.
    searches={}
    cones=" OR ".join("(1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',"+str(r["ra"])+","+str(r["dec"])+",0.006)))" for r in rows)
    for label,table in [("XMM","xmm4dr14"),("Chandra","csc"),("Swift","swift2sxps")]:
        query="SELECT TOP 1000 name,ra,dec FROM "+table+" WHERE "+cones
        s=requests.Session();s.headers.update(HEAD)
        try:
            b,receipt=post_fetch(s,"https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync",
                 {"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":query},limit=4_000_000,timeout=(12,65))
            if b.lstrip().startswith(b"<"):
                raw=b.decode(errors="replace")
                m=re.search(r'<INFO[^>]*name="QUERY_STATUS"[^>]*value="([^"]+)"[^>]*>(.*?)</INFO>',raw,re.I|re.S)
                if m and m.group(1).upper()!="OK":
                    raise ValueError("TAP_QUERY_STATUS_"+m.group(1)+": "+re.sub("<[^>]+>"," ",m.group(2))[:300])
                try:
                    tb=Table.read(io.BytesIO(b),format="votable")
                    rr=[{k:str(row[k]) for k in tb.colnames} for row in tb]
                    colnames={k.lower() for k in tb.colnames}
                except Exception as e:
                    raise ValueError("TAP_VOTABLE_parse_failed:"+str(e)[:200]+
                                     " QUERY_STATUS_excerpt:"+str(m.group(0)[:200] if m else raw[-400:]))
            elif b.lstrip().startswith(b"ERROR"):
                raise ValueError("TAP_plain_error:"+b[:350].decode(errors="replace"))
            else:
                rr=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig"))))
                colnames={k.lower() for k in rr[0]} if rr else set(b.decode().split("\n")[0].lower().split(","))
            if "ra" not in colnames or "dec" not in colnames:
                raise ValueError("TAP_no_parseable_RA_Dec_or_empty_columns:"+str(colnames))
            searches[label]={"status":"queried","receipt":receipt,"rows_returned":len(rr),"matches":{}}
            for x in rows:
                c=SkyCoord(float(x["ra"])*u.deg,float(x["dec"])*u.deg)
                found=[]
                for obj in rr:
                    try:
                        co=SkyCoord(float(obj["ra"])*u.deg,float(obj["dec"])*u.deg)
                        sep=c.separation(co).arcsec
                        if sep<max(12.0,2*float(x["poserr_arcsec"])): found.append({"name":obj.get("name"),"sep_arcsec":round(float(sep),2)})
                    except Exception:continue
                searches[label]["matches"][x["uid"]]=sorted(found,key=lambda z:z["sep_arcsec"])[:5]
        except Exception as e: searches[label]={"status":"query_failed_not_nondetection","error":str(e)[:250]}
        finally:s.close()
    return searches
def atnf(rows):
    s=requests.Session();s.headers.update(HEAD)
    url="https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz"
    try:
        b,rec=fetch(s,url,limit=10_000_000,timeout=(12,70))
        with tarfile.open(fileobj=io.BytesIO(b),mode="r:gz") as t:
            choices=[m for m in t if m.isfile() and m.name.endswith("psrcat.db")]
            if not choices:raise ValueError("psrcat.db missing")
            data=t.extractfile(choices[0]).read().decode(errors="replace")
        psr=[]
        for block in data.split("@"):
            fields={}
            for line in block.splitlines():
                line=line.strip()
                if not line or line.startswith("#"):continue
                spl=line.split()
                if len(spl)>=2 and spl[0] in ("RAJ","DECJ","PSRJ","PSRB"):
                    fields[spl[0]]=spl[1]
            if "RAJ" not in fields or "DECJ" not in fields:continue
            try:
                c=SkyCoord(fields["RAJ"],fields["DECJ"],unit=(u.hourangle,u.deg))
                psr.append((c,fields.get("PSRJ",fields.get("PSRB","?"))))
            except Exception:continue
        found={}
        for x in rows:
            c=SkyCoord(float(x["ra"])*u.deg,float(x["dec"])*u.deg)
            found[x["uid"]]=[{"pulsar":n,"sep_arcsec":round(float(c.separation(p).arcsec),2)}
              for p,n in psr if c.separation(p).arcsec<30][:12]
        return {"status":"queried","receipt":rec,"parsed_pulsars":len(psr),"matches":found}
    except Exception as e:return {"status":"query_failed_not_nondetection","error":str(e)[:260]}
    finally:s.close()
def kurpas(rows):
    s=requests.Session();s.headers.update(HEAD)
    url="https://vizier.cds.unistra.fr/viz-bin/asu-tsv"
    try:
        b,rec=fetch(s,url,{"-source":"J/A+A/687/A251","-out.max":"1000"},limit=2_000_000,timeout=(10,26))
        content=b.decode(errors="replace")
        # Presence of a catalogue page does not mean the candidates were parsed.
        return {"status":"unparsed_public_catalogue_needs_manual_novelty_verification",
                "receipt":rec,"header_excerpt":content[:400]}
    except Exception as e:return {"status":"unavailable_novelty_unverified","error":str(e)[:220]}
    finally:s.close()
def main():
    start=datetime.now(timezone.utc).isoformat()
    frozen=readcsv(PREV/"official_top60.csv")
    no=[r for r in frozen if r["status"]=="NO_SECURE_ASSOCIATION_LIMITS_UNVERIFIED"]
    original={r["uid"]:r for r in readcsv(PREV/"screened_main.csv")}
    assert len(frozen)==60 and len(no)==35 and len(set(x["uid"] for x in no))==35
    for x in no:assert x["uid"] in original
    for x in no:x.update({k:original[x["uid"]][k] for k in ("poserr_arcsec","soft_fraction","p1_counts","det_like","brighter_neighbors_60arcsec","ml_counts")})
    evidence={"started_utc":start,
       "frozen_input":{"official_top60_sha256":sha(PREV/"official_top60.csv"),
         "screened_main_sha256":sha(PREV/"screened_main.csv"),"count":len(no)},
       "warning":"Positions are eROSITA not optical; non-association is not an upper limit."}
    # Reuse previously hashed native FITS-derived photometry and only refresh
    # the dust sightlines. Avoid redownloading large coadds / rate-limited service.
    imgpath=OUT/"image_receipts.json"
    if imgpath.exists():
        old=json.loads(imgpath.read_text())
        if [v.get("uid") for v in old]==[v["uid"] for v in no]:
            images=old
            for row,im in zip(no,images):
                with requests.Session() as ss:
                    ss.headers.update(HEAD)
                    im["galactic_dust"]=attempt("galactic_dust",lambda: dust(ss,row))
        else:
            with ThreadPoolExecutor(max_workers=5) as pool:images=list(pool.map(one,no))
    else:
        with ThreadPoolExecutor(max_workers=5) as pool:images=list(pool.map(one,no))
    writejson("image_receipts.json",images)
    arch=heasarc_archives(no)
    writejson("archive_receipts.json",arch)
    at=atnf(no)
    writejson("atnf_receipts.json",at)
    lit=kurpas(no)
    writejson("literature_receipts.json",lit)
    out=[]
    for x,im in zip(no,images):
        assert x["uid"]==im["uid"]
        observed=[]
        for survey in ("LS10","PS1","DES_DR2"):
            for v in im.get(survey,{}).get("bands",{}).values():
                if v.get("status")=="measured" and v.get("band") in ("g","r"):
                    observed.append(v)
        # Conservative: only a true image-level <5sigma measurement yields a bound.
        nondet=[v for v in observed if not v["detected_at_position_5sigma"] and v["five_sigma_nmg_ab"]>0]
        detected=[v for v in observed if v["detected_at_position_5sigma"]]
        g=[v for v in nondet if v["band"]=="g"]
        chosen=max(g,key=lambda z:z["five_sigma_ab_mag"]) if g else None
        # f_X/(nu_g * f_nu(g 5sigma upper limit)); quoted observed-band, uncorrected
        # for Galactic extinction or absorbed-X-ray model. Never use g=V shortcut.
        ratio=(float(x["fX_0p2_2p3"])/(NU_G*AB_NMGY_FNU*chosen["five_sigma_nmg_ab"])
               if chosen else None)
        fx=float(x["fX_0p2_2p3"])
        lum={str(d):round(4*math.pi*(d*3.085677581e21)**2*fx,4) for d in (0.3,1,50)}
        known=KNOWN.get(x["uid"])
        classif="known INS" if known else "unclassified (INS claim not established)"
        out.append({"uid":x["uid"],"name":x["iauname"],"ra":float(x["ra"]),"dec":float(x["dec"]),
            "classification":classif,"known_identity":known,
            "new_INS_candidate_validated":False,
            "eROSITA":{"poserr_arcsec":finite(x["poserr_arcsec"]),"extent_arcsec":0,
                "extent_likelihood_max_6_by_selection":True,"source_counts":finite(x["ml_counts"]),
                "detection_likelihood":finite(x["det_like"]),"soft_fraction_P1P2_over_P123":finite(x["soft_fraction"]),
                "P1_counts":finite(x["p1_counts"]),"native_blackbody_vs_powerlaw_fit":"NOT PERFORMED",
                "variability":"NOT ESTABLISHED",
                "brighter_neighbour_count":int(x["brighter_neighbors_60arcsec"])},
            "optical":{"bands_with_position_forced_photometry":len(observed),
               "5sigma_detections_in_fixed_aperture":len(detected),
               "5sigma_nondetections_in_fixed_aperture":len(nondet),
               "best_g_5sigma_limit_AB":chosen["five_sigma_ab_mag"] if chosen else None,
               "g_observed_fX_over_nu_fnu_lower_bound":round(ratio,1) if ratio else None,
               "g_bound_exceeds_1e4":bool(ratio and ratio>1e4),
               "warning":"Aperture approximation, crowded-field confusion, Galactic extinction uncorrected; not a PSF/model fit"},
            "infrared":{"unWISE_bands":im.get("unWISE",{}).get("bands",{}),
                "status":im.get("unWISE",{}).get("status","unknown")},
            "galactic":{"coordinates_lb_deg":im["galactic_l_b_deg"],
                "dust":im.get("galactic_dust"),"Lx_erg_s_if_kpc":lum,
                "distance_warning":"50 kpc only illustrative; do not assign LMC membership from angular position"},
            "archival":{"detections":{n:v.get("matches",{}).get(x["uid"],[]) if v.get("status")=="queried" else None for n,v in arch.items()},
               "query_status":{n:v.get("status") for n,v in arch.items()},
               "variability":"Not measurable from positional matches alone"},
            "ATNF":{"status":at.get("status"),"matches":at.get("matches",{}).get(x["uid"]) if at.get("status")=="queried" else None},
            "literature":{"status":lit.get("status"),"known_M7":x["uid"] in ("111112000001","112113200001","119606900001","128212900001"),
               "Kurpas_2024_novelty":"NOT VERIFIED unless published named control"},
            "decision":"Known XDINS positive control, not discovery" if known
                      else "Do not claim AGN/star/cluster/spurious/novel INS from these data alone"})
    assert len(out)==35 and len([x for x in out if x["classification"]=="known INS"])==6
    assert not any(x["new_INS_candidate_validated"] for x in out)
    with open(OUT/"classifications.csv","w",newline="") as f:
        cols=["uid","name","classification","known_identity","eROSITA_poserr_arcsec","eROSITA_soft_fraction",
              "forced_optical_bands","g_5sigma_limit_AB","g_fX_over_nu_fnu_lower_bound","ratio_gt_1e4",
              "LS10_status","PS1_status","DES_status","unWISE_status","XMM_matches","Chandra_matches","Swift_matches","ATNF_matches","literature_status"]
        w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
        for x,im in zip(out,images):
            ar=x["archival"]["detections"]
            w.writerow(dict(uid=x["uid"],name=x["name"],classification=x["classification"],
                known_identity=x["known_identity"],eROSITA_poserr_arcsec=x["eROSITA"]["poserr_arcsec"],
                eROSITA_soft_fraction=x["eROSITA"]["soft_fraction_P1P2_over_P123"],
                forced_optical_bands=x["optical"]["bands_with_position_forced_photometry"],
                g_5sigma_limit_AB=x["optical"]["best_g_5sigma_limit_AB"],
                g_fX_over_nu_fnu_lower_bound=x["optical"]["g_observed_fX_over_nu_fnu_lower_bound"],
                ratio_gt_1e4=x["optical"]["g_bound_exceeds_1e4"],
                LS10_status=im["LS10"]["status"],PS1_status=im["PS1"]["status"],
                DES_status=im["DES_DR2"]["status"],unWISE_status=im["unWISE"]["status"],
                XMM_matches=None if ar["XMM"] is None else len(ar["XMM"]),
                Chandra_matches=None if ar["Chandra"] is None else len(ar["Chandra"]),
                Swift_matches=None if ar["Swift"] is None else len(ar["Swift"]),
                ATNF_matches=None if x["ATNF"]["matches"] is None else len(x["ATNF"]["matches"]),
                literature_status=x["literature"]["status"]))
    writejson("source_evidence.json",out)
    summary={"finished_utc":datetime.now(timezone.utc).isoformat(),
       "frozen_input":evidence["frozen_input"],"source_count":35,
       "known_INS":6,"known_M7":4,"known_radio_pulsars":2,"new_validated_INS":0,
       "still_unclassified":29,
       "with_position_forced_optical":sum(x["optical"]["bands_with_position_forced_photometry"]>0 for x in out),
       "with_g_lower_bound":sum(x["optical"]["g_observed_fX_over_nu_fnu_lower_bound"] is not None for x in out),
       "provisional_g_bounds_gt_1e4":sum(x["optical"]["g_bound_exceeds_1e4"] for x in out),
       "archival_query_status":{n:v["status"] for n,v in arch.items()},
       "ATNF_status":at["status"],"Kurpas_status":lit["status"],
       "native_X_ray_spectra_fit":False,"long_term_variability_established":False,
       "NH_spectral_fit":False,"novelty_audit_complete":False,
       "verdict":"NO NEW DISCOVERY CLAIM; incomplete independent multiwavelength vetting"}
    writejson("summary.json",summary)
    lines=["# Referee outcome: no independently demonstrated NEW XDINS","",
     "This file is generated from explicit source evidence, not catalogue blank fields. Two radio pulsars are correctly counted as known INS, not new XDINS.","",
     "| Check | Outcome |","|---|---|"]
    for k,v in summary.items():
        if k not in ("frozen_input","verdict"):lines.append("| "+k+" | "+str(v).replace("|","/")+" |")
    lines+=["","**Decision:** Four known XDINS plus two catalogued rotation-powered radio pulsars are not novel. The other 29 remain unclassified; no source may be promoted without calibrated optical depth and corrected flux ratio, native eROSITA/XMM/Chandra/Swift spectral response fits, variability check, Galactic N_H, and complete prior-publication comparison. Optical aperture measurements are provisional where blends, PSF and image coverage are unvalidated. Failed upstream requests are not negative detections.","",
    "Input receipts, per-source measurements and service failures: [source_evidence.json](source_evidence.json), [image_receipts.json](image_receipts.json), [archive_receipts.json](archive_receipts.json), [atnf_receipts.json](atnf_receipts.json), [literature_receipts.json](literature_receipts.json)."]
    (OUT/"REFEREE.md").write_text("\n".join(lines)+"\n")
    print(json.dumps(summary,sort_keys=True))
if __name__=="__main__":main()
