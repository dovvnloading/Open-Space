"""Fit trailed point sources to individual ZTF differences and check the science-image WCS against Gaia DR3."""
import csv, gzip, hashlib, io, json, math, sys
from datetime import datetime,timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request,urlopen
import numpy as np
from scipy.ndimage import convolve
from scipy.optimize import least_squares

HERE=Path(__file__).parent
OUT=HERE/"image-search"/"centroids"
OUT.mkdir(exist_ok=True)
sys.path.insert(0,str(HERE/"image-search"/"stack"))
import stack as base
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

catalog=json.loads((HERE/"image-search"/"archive-results.json").read_text())["ztf"]["ztf_nominal_footprints"]
meta=json.loads((HERE/"image-search"/"stack"/"stack-results.json").read_text())["frames"]
rows={r["filefracday"]:r for r in catalog}
receipts=[]
def get_fits(row,suffix):
    url=base.image_url(row,suffix)
    p=urlencode({"center":f"{row['pred_ra_deg']:.7f},{row['pred_dec_deg']:.7f}","size":"800arcsec","gzip":"false"})
    for url2 in (url+"?"+p,url):
        try:
            raw,rec=base.request(url2)
            receipts.append(rec)
            return base.get_array(raw)
        except Exception as e:
            receipts.append({"url":url2,"error":str(e)})
    raise ValueError("FITS unavailable")

def make_psf(x,y,xc,yc,sigma,trail):
    ts=np.linspace(-.5,.5,9)
    profiles=[np.exp(-((x-xc-t*trail[0])**2+(y-yc-t*trail[1])**2)/(2*sigma**2)) for t in ts]
    return np.mean(profiles,axis=0)/(2*np.pi*sigma**2)

def flux_map(data,sigma,trail):
    y,x=np.mgrid[-8:9,-8:9]
    p=make_psf(x,y,0,0,sigma,trail)
    p-=p.mean()
    kernel=p/(p*p).sum()
    valid=np.isfinite(data)
    a=convolve(np.where(valid,data,0),kernel,mode="constant",cval=0)
    from scipy.ndimage import minimum_filter
    a[~minimum_filter(valid,size=17,mode="constant",cval=0)]=np.nan
    return a

def centroid(data,x0,y0,sigma,trail,bound=4.0):
    x0,y0=float(x0),float(y0)
    xi,yi=round(x0),round(y0)
    if xi<8 or yi<8 or xi+9>data.shape[1] or yi+9>data.shape[0]:
        return None
    sub=data[yi-8:yi+9,xi-8:xi+9]
    if not np.isfinite(sub).all():return None
    yy,xx=np.mgrid[yi-8:yi+9,xi-8:xi+9]
    x,y=xx-x0,yy-y0
    background=float(np.median(sub))
    noise=float(1.4826*np.median(abs(sub-background)))
    if not noise>0:return None
    def calc(p):
        dx,dy,flux,b,gx,gy=p
        return flux*make_psf(x,y,dx,dy,sigma,trail)+b+gx*x+gy*y
    def residual(p):return ((calc(p)-sub)/noise).ravel()
    init=[0.,0.,float(np.max(sub)-background)*2*np.pi*sigma*sigma,background,0.,0.]
    fit=least_squares(residual,init,bounds=([-bound,-bound,0,-np.inf,-np.inf,-np.inf],
                                           [bound,bound,np.inf,np.inf,np.inf,np.inf]),max_nfev=120)
    dx,dy,flux=fit.x[:3]
    cov=np.linalg.pinv(fit.jac.T@fit.jac)
    rms=float(np.sqrt(np.mean(residual(fit.x)**2)))
    return {"x":float(x0+dx),"y":float(y0+dy),"dx":float(dx),"dy":float(dy),
            "formal_xy_sigmas":np.sqrt(np.maximum(0,np.diag(cov)[:2])).tolist(),
            "flux":float(flux),"formal_flux_snr":float(flux/max(1e-12,np.sqrt(cov[2,2]))),
            "noise_scaled_rms":rms,"boundary":max(abs(dx),abs(dy))>bound-.1,
            "stamp_peak":float(sub.max()),"sky_sigma":noise}

def gaia(date,ra,dec):
    params={"-source":"I/355/gaiadr3","-c":f"{ra} {dec}","-c.rs":"340",
            "-out":"Source,RA_ICRS,DE_ICRS,Gmag,pmRA,pmDE","-out.max":"4000"}
    uri="https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"+urlencode(params)
    for host in ("vizier.cds.unistra.fr","vizier.cfa.harvard.edu"):
        try:
            url=uri.replace("vizier.cds.unistra.fr",host)
            raw,rec=base.request(url,60)
            receipts.append(rec)
            (OUT/(date+"-Gaia.tsv.gz")).write_bytes(gzip.compress(raw,mtime=0))
            lines=[l for l in raw.decode().splitlines() if l.strip() and not l.startswith("#")]
            stars=[]
            for r0 in csv.DictReader(lines,delimiter="\t"):
                r={str(k).strip():str(v).strip() for k,v in r0.items()}
                try:
                    g=float(r["Gmag"])
                    if not 15<g<18.5:continue
                    ra,dec=float(r["RA_ICRS"]),float(r["DE_ICRS"])
                    if not r["pmRA"] or not r["pmDE"]:continue
                    stars.append({"ra":ra,"dec":dec,"g":g,"pmra":float(r["pmRA"]),"pmdec":float(r["pmDE"])})
                except (ValueError,KeyError):continue
            if not stars: raise ValueError("No usable Gaia rows")
            return sorted(stars,key=lambda s:s["g"])
        except Exception as e:
            receipts.append({"url":url,"error":str(e)})
    return []

def calibration(f,stars):
    try:
        dat,hdr=get_fits(f["row"],"sciimg.fits")
        w=WCS(hdr).celestial
        years=(f["jd"]-2457388.5)/365.25
        residuals=[]
        for s in stars:
            ra=s["ra"]+s["pmra"]*years/(3.6e6*np.cos(np.radians(s["dec"])))
            dec=s["dec"]+s["pmdec"]*years/3.6e6
            x,y=w.all_world2pix(ra,dec,0)
            fit=centroid(dat,x,y,f["sigma"],[0,0],2.5)
            if not fit or fit["boundary"] or fit["formal_flux_snr"]<10 or fit["noise_scaled_rms"]>4:continue
            a,b=w.all_pix2world(fit["x"],fit["y"],0)
            residuals.append([3600*(ra-a)*np.cos(np.radians(dec)),3600*(dec-b)])
            if len(residuals)>=75:break
        if len(residuals)<5:return {"n":len(residuals),"error":"Too few calibrated stars"}
        v=np.array(residuals);m=np.median(v,axis=0);mad=1.4826*np.median(abs(v-m),axis=0)
        return {"n":len(v),"observed_to_Gaia_EN_arcsec":m.tolist(),"single_star_MAD_arcsec":mad.tolist()}
    except Exception as e:return {"error":str(e)}

results={}
for date in ("20251027","20251029"):
    subset=sorted([m for m in meta if m["file"].startswith(date) and int(m["infobits"])<33554432],key=lambda m:m["jd"])
    stars=gaia(date,subset[0]["ra"],subset[0]["dec"])
    frames=[]
    for m in subset:
        try:
            row=rows[m["file"]]
            data,hdr=get_fits(row,"scimrefdiffimg.fits.fz")
            w=WCS(hdr).celestial
            x,y=map(float,w.all_world2pix(m["ra"],m["dec"],0))
            xa,ya=max(0,round(x)-402),max(0,round(y)-402)
            dat=data[ya:round(y)+403,xa:round(x)+403]
            scale=float(np.mean(proj_plane_pixel_scales(w))*3600)
            sig=m["seeing"]/scale/2.35482
            xp,yp=w.all_world2pix([subset[0]["ra"],subset[-1]["ra"]],[subset[0]["dec"],subset[-1]["dec"]],0)
            seconds=(subset[-1]["jd"]-subset[0]["jd"])*86400
            trail=[float((xp[1]-xp[0])*hdr["EXPTIME"]/seconds),float((yp[1]-yp[0])*hdr["EXPTIME"]/seconds)]
            f={**m,"row":row,"data":dat,"wcs":w,"x0":xa,"y0":ya,"scale":scale,"sigma":sig,"trail":trail,
               "map":flux_map(dat,sig,trail)}
            frames.append(f)
        except Exception as e:print("FRAME ISSUE",m["file"],str(e),flush=True)
    if len(frames)!=len(subset):
        results[date]={"failure":"Missing target frames"};continue
    joint=base.scan(frames)
    offset=joint["best_east_north_arcsec"]
    per=[]
    for f in frames:
        ra0=f["ra"]+offset[0]/(3600*np.cos(np.radians(f["dec"])))
        dec0=f["dec"]+offset[1]/3600
        x,y=f["wcs"].all_world2pix(ra0,dec0,0)
        fit=centroid(f["data"],x-f["x0"],y-f["y0"],f["sigma"],f["trail"])
        observed=calibration(f,stars)
        if fit:
            ra,dec=f["wcs"].all_pix2world(fit["x"]+f["x0"],fit["y"]+f["y0"],0)
            fit.update(uncorrected_fit_RA_DEC_deg=[float(ra),float(dec)],
                       measured_minus_ephemeris_EN_arcsec=[float((ra-f["ra"])*3600*np.cos(np.radians(dec))),float((dec-f["dec"])*3600)])
        per.append({"file":f["file"],"jd":f["jd"],"utc":datetime.fromtimestamp((f["jd"]-2440587.5)*86400,tz=timezone.utc).isoformat(),
                    "seeing":f["seeing"],"trail_xy_pixels":f["trail"],"scale_arcsec_pix":f["scale"],
                    "PSF_fit":fit,"local_gaia":observed})
    results[date]={"joint_trailed_PSF_filter":joint,"per_exposure_diagnostics":per,"gaia_stars_queried":len(stars)}
    print("CENTROIDS",date,json.dumps(results[date]),flush=True)
(OUT/"measurements.json").write_text(json.dumps({"dates":results,"receipts":receipts},indent=2,allow_nan=False)+"\n")
