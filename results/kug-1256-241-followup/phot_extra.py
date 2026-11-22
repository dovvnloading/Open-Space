"""Local-noise radio measurements and Swift UVOT aperture count rates for one nucleus."""
import csv
import gzip
import hashlib
import io
import math
import re
import zlib
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from urllib.parse import urljoin, urlsplit

import numpy as np
import requests
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from bs4 import BeautifulSoup
from scipy.stats import chi2

import common
from common import HOST_DEC, HOST_RA, OUT, event, fetch, status
from photometry import col, epoch, tap, num, off, HEASARC, IRSA

LOG = []
NRAO = "https://archive-new.nrao.edu/vlass/quicklook/"


def record(name, state, detail):
    LOG.append((name, state, str(detail)[:650]))
    status(name, state, str(detail)[:650])


def links(raw, base):
    return [urljoin(base, a["href"]) for a in BeautifulSoup(raw,"html.parser").find_all("a",href=True)
            if a["href"] not in ("../","./") and not a["href"].startswith(("?","#"))]


def basename(url):
    return urlsplit(url).path.rstrip("/").rsplit("/",1)[-1]


def bytes_range(key, url, lo, hi):
    if hi-lo > 12_000_000:
        raise ValueError("Requested stripe exceeds size bound")
    h = {"Range":f"bytes={lo}-{hi}","Accept-Encoding":"identity"}
    rec = {"key":key,"url":url,"method":"GET","headers":h,"utc":datetime.now(timezone.utc).isoformat()}
    try:
        with requests.get(url,headers=h,timeout=(12,40),stream=True) as r:
            rec.update(http=r.status_code,final_url=r.url)
            if r.status_code != 206 or not r.headers.get("Content-Range","").startswith(f"bytes {lo}-"):
                raise ValueError(f"range response unavailable ({r.status_code})")
            raw = r.raw.read(hi-lo+2)
            if len(raw) > hi-lo+1:
                raise ValueError("range exceeded")
            rec.update(bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest())
            return raw
    except Exception as exc:
        rec["error"] = str(exc)[:300]
        raise
    finally:
        with common.lock:
            common.receipts.append(rec)


def crop_fits(url, key):
    data = bytes_range(key+"_header",url,0,57599)
    end = next((i for i in range(0,len(data),80) if data[i:i+8] == b"END     "),None)
    if end is None:
        raise ValueError("Primary FITS header exceeds 57600 bytes")
    h = fits.Header.fromstring(data[:end+80].decode("ascii"),sep="")
    off = math.ceil((end+80)/2880)*2880
    size,ny,bit = int(h["NAXIS1"]),int(h["NAXIS2"]),int(h["BITPIX"])
    dt = {-32:">f4",-64:">f8",16:">i2",32:">i4"}[bit]
    wcs = WCS(h,naxis=2)
    x,y = map(float,wcs.world_to_pixel_values(HOST_RA,HOST_DEC))
    if not (45<=x<size-45 and 45<=y<ny-45):
        raise ValueError("nucleus outside radio tile's interior")
    scale = float(np.mean(proj_plane_pixel_scales(wcs)))*3600
    half = min(128,math.ceil(52/scale))
    x0,x1=max(0,round(x)-half),min(size,round(x)+half+1)
    y0,y1=max(0,round(y)-half),min(ny,round(y)+half+1)
    width = size*(abs(bit)//8)
    stripe = bytes_range(key+"_stripe",url,off+y0*width,off+y1*width-1)
    image = np.frombuffer(stripe,dtype=dt).reshape(y1-y0,size)[:,x0:x1].astype(float)
    image = image*float(h.get("BSCALE",1))+float(h.get("BZERO",0))
    return image,x-x0,y-y0,scale,h


def radio_point(name,a,x,y,scale,h,url,t_override="",end=""):
    yy,xx = np.indices(a.shape)
    radius = np.hypot(xx-x,yy-y)*scale
    beam = 3600*max(float(h.get("BMAJ",0)),float(h.get("BMIN",0)))
    sky = a[(radius>=max(12,3*beam)) & (radius<50) & np.isfinite(a)]
    if len(sky)<80 or not np.isfinite(a[round(y),round(x)]):
        raise ValueError("Radio cutout lacks finite central pixels/annulus")
    unit = str(h.get("BUNIT","")).lower().replace(" ","")
    if unit in ("jy/beam", "jybeam-1"):
        conv = 1000.
    elif unit in ("mjy/beam","mjybeam-1"):
        conv = 1.
    else:
        raise ValueError("Radio image units "+unit+" unsupported")
    bg = np.median(sky)
    sigma = conv*1.4826*np.median(np.abs(sky-bg))
    measured = conv*(float(a[round(y),round(x)])-bg)
    if not sigma>0:
        raise ValueError("Local radio noise unavailable")
    t = t_override or epoch(h.get("MJD-OBS",h.get("DATE-OBS","")))
    band = "3 GHz" if "VLASS" in name else "1.4 GHz"
    kind, v = ("detection", measured) if measured>=5*sigma else ("upper_limit", max(0,measured)+5*sigma)
    event(name,band,t,v,sigma,"mJy/beam",kind,source=url, end_mjd=end,
          note=f"Fixed nearest nucleus pixel {measured:.4g} mJy/beam after annulus median; annulus MAD sigma {sigma:.4g}. Limit=max(pixel,0)+5sigma. Beam major {beam:.2f} arcsec.")
    return f"{kind} {v:.4g} mJy/beam; local rms {sigma:.4g}, date from FITS observation keyword when present."


def tile_skycoord(url):
    m=re.search(r"J(\d\d)(\d\d)(\d\d)([+-])(\d\d)(\d\d)(\d\d)?",basename(url))
    if not m:
        return math.inf
    h, minute, sec, sign, deg, arcmin, arcsec = m.groups()
    ra=15*(int(h)+int(minute)/60+int(sec)/3600)
    dec=(int(deg)+int(arcmin)/60+int(arcsec or 0)/3600)*(1 if sign=="+" else -1)
    return math.hypot((ra-HOST_RA)*math.cos(math.radians(HOST_DEC)),dec-HOST_DEC)


def vlass():
    roots=links(fetch("nrao_vlass_quicklook_epochs",NRAO,timeout=25),NRAO)
    epochs={}
    for u in roots:
        name=basename(u)
        m=re.fullmatch(r"VLASS([1234]\.[12])(?:v(\d+))?",name,re.I)
        if m:
            k=m[1]
            if k not in epochs or int(m[2] or 1)>epochs[k][0]:
                epochs[k]=int(m[2] or 1),u.rstrip("/")+"/"
    if not epochs:
        raise ValueError("No VLASS quicklook observing epochs in NRAO listing")
    # Empirical footprints from the deposited first pass identify the precise VLASS tile.
    field_rows=[r for r in common.events if r["archive"]=="VLASS CADC fields"]
    def keyinfo(k):
        return [r for r in field_rows if str(r["source"]).startswith("VLASS"+k+".")]
    def half_epoch(pair):
        k, (_, root) = pair
        try:
            fields=keyinfo(k)
            if not fields:
                return
            tile_ids={r["source"].split(".")[2] for r in fields}
            tiles=[root+t+"/" for t in sorted(tile_ids)]
            if not tiles:
                return record("VLASS "+k,"no_matching_tile","No adjacent target tile IDs in quicklook listing")
            subt=[]
            for tile in tiles:
                tile=tile.rstrip("/")+"/"
                try:
                    subt += [(tile_skycoord(u),u) for u in links(fetch("vlass_"+k+"_"+basename(tile),tile,timeout=20),tile) if tile_skycoord(u)<.78]
                except Exception as e:
                    LOG.append(("VLASS "+k,"tile listing failed",str(e)[:200]))
            if not subt:
                return record("VLASS "+k,"no_matching_tile","No target-centred subtiles within 0.78 degree in adjacent NRAO tile folders checked.")
            for j, (_,u) in enumerate(sorted(set(subt))):
                try:
                    if u.lower().endswith(".fits"):
                        products=[u]
                    else:
                        u=u.rstrip("/")+"/"
                        products=links(fetch(f"vlass_files_{k}_{j}",u,timeout=20),u)
                    fitsurls=[p for p in products if re.search(r"\.I\..*image\.pbcor\.tt0(?:\.subim)?\.fits$",basename(p),re.I)]
                    if not fitsurls:
                        continue
                    url=sorted(fitsurls)[-1]
                    a,x,y,scale,h=crop_fits(url,f"vlass_{k}_{j}")
                    timed=next((r for r in fields if r["mjd"]!=""),{})
                    return record("VLASS "+k,"measured",radio_point("VLASS "+k,a,x,y,scale,h,url,timed.get("mjd",""),timed.get("end_mjd","")))
                except Exception as e:
                    LOG.append(("VLASS "+k,"subtile issue",str(e)[:220]))
            record("VLASS "+k,"unavailable",f"{len(subt)} nearby subtiles, none yielded a calibrated usable nucleus cutout")
        except Exception as e:
            record("VLASS "+k,"unavailable",e)
    with ThreadPoolExecutor(max_workers=6) as p:
        list(p.map(half_epoch,epochs.items()))


def first():
    url="https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
    params={"Position":f"{HOST_RA},{HOST_DEC}","Survey":"VLA FIRST (1.4 GHz)","Coordinates":"J2000","Projection":"Tan","Size":0.032,"Pixels":181,"Return":"FITS"}
    raw=fetch("first_skyview_cutout",url,params=params,timeout=55)
    if not raw.startswith(b"SIMPLE"):
        raise ValueError("SkyView did not return a FIRST FITS image")
    with fits.open(io.BytesIO(raw)) as hd:
        h=hd[0].header
        a=np.squeeze(hd[0].data)
        wcs=WCS(h,naxis=2)
        x,y=map(float,wcs.world_to_pixel_values(HOST_RA,HOST_DEC))
        scale=float(np.mean(proj_plane_pixel_scales(wcs)))*3600
        record("FIRST image","measured",radio_point("FIRST image",a,x,y,scale,h,url))


def observation_id(v):
    s = str(v).strip()
    return s.zfill(11) if re.fullmatch(r"\d{7,11}",s) else ""


def qfield():
    q=(f"SELECT * FROM swiftmastr WHERE CONTAINS(POINT('ICRS',ra,dec),"
       f"CIRCLE('ICRS',{HOST_RA},{HOST_DEC},0.21))=1")
    for i,server in enumerate([HEASARC,"https://heasarc.gsfc.nasa.gov/vo/tap/sync"]):
        try:
            return tap("swift_mission_fields_"+str(i),server,q)
        except Exception:
            if i:
                raise


def image_measure(img,exps,root,filt,obs):
    dates=0
    for ext in img:
        image=ext.data
        if image is None or np.ndim(image)!=2:
            continue
        h=ext.header
        name=h.get("EXTNAME")
        candidates=[u for u in exps if u.data is not None and np.ndim(u.data)==2 and
                    (name==u.header.get("EXTNAME") or (h.get("TSTART") is not None and abs(float(h["TSTART"])-float(u.header.get("TSTART",-1e10)))<1.))]
        if not candidates:
            continue
        if len(candidates)>1:
            raise ValueError("Ambiguous exposure map for sky image extension")
        em=candidates[0]
        wh=WCS(h,naxis=2)
        ew=WCS(em.header,naxis=2)
        x,y=map(float,wh.world_to_pixel_values(HOST_RA,HOST_DEC))
        if not math.isfinite(x+y):
            continue
        scale=float(np.mean(proj_plane_pixel_scales(wh)))*3600
        n=math.ceil(27/scale)
        x0,x1=round(x)-n,round(x)+n+1
        y0,y1=round(y)-n,round(y)+n+1
        if x0<0 or y0<0 or x1>image.shape[1] or y1>image.shape[0]:
            continue
        yy,xx=np.mgrid[y0:y1,x0:x1]
        ex,ey=ew.world_to_pixel_values(*wh.pixel_to_world_values(xx,yy))
        ex,ey=np.rint(ex).astype(int),np.rint(ey).astype(int)
        if np.any((ex<0)|(ey<0)|(ey>=em.data.shape[0])|(ex>=em.data.shape[1])):
            continue
        expo=np.asarray(em.data[ey,ex],dtype=float)
        counts=np.asarray(image[y0:y1,x0:x1],dtype=float)
        rr=np.hypot(xx-x,yy-y)*scale
        src=rr<=5
        ann=(rr>=15)&(rr<=25)
        good=(np.isfinite(counts)&np.isfinite(expo)&(counts>=0)&(expo>0))
        if not good[src].all() or good[ann].mean()<.95:
            continue
        ann &= good
        units=str(h.get("BUNIT","")).lower().strip()
        if units in ("count/s","counts/s","count/s/pixel","counts/s/pixel","count s-1"):
            counts=counts*expo
        elif units not in ("count","counts","count/pixel","counts/pixel"):
            raise ValueError("Unrecognized UVOT BUNIT "+repr(units))
        if np.min(expo[src]) < .85*np.max(expo[src]):
            continue
        rate=float(np.sum(counts[src]/expo[src]))
        bgrate=float(np.sum(counts[ann])/np.sum(expo[ann]))
        net=rate-src.sum()*bgrate
        err=math.sqrt(float(np.sum(counts[src]/expo[src]**2))+src.sum()**2*float(np.sum(counts[ann]))/np.sum(expo[ann])**2)
        # Report conservative one-sided Poisson 99.865% gross upper bound. No background subtraction for limits.
        sig=net/err if err else 0.
        kind,val=("detection",net) if sig>=5 else ("upper_limit",float(.5*chi2.ppf(.99865,2*(math.ceil(np.sum(counts[src]))+1))/np.min(expo[src])))
        t=epoch(h.get("DATE-OBS",img[0].header.get("DATE-OBS","")))
        end=epoch(h.get("DATE-END",""))
        event("Swift UVOT image",filt,t,val,err,"count/s",kind,source=obs+"/"+str(name),
              note=f"5 arcsec instrumental aperture, 15–25 arcsec local sky, gross-count conservative limit; net={net:.5g}. Coincidence loss and long-term sensitivity uncorrected; {root}.",end_mjd=end)
        dates+=1
    return dates


def swift():
    rr=qfield()
    visits={}
    for r in rr:
        oid=observation_id(col(r,"obsid","obs_id","observation_id","sequence_number"))
        t=epoch(col(r,"time","start_time","mjd_start","date_obs"))
        if oid and t:
            visits[oid]=Time(t,format="mjd").strftime("%Y_%m")
    if rr and not visits:
        raise ValueError("Could not resolve Swift observation IDs/dates from columns "+",".join(rr[0])[:400])
    n=0
    def visit(pair):
        oid,month=pair
        root=f"https://heasarc.gsfc.nasa.gov/FTP/swift/data/obs/{month}/{oid}/uvot/image/"
        try:
            listing=links(fetch("uvot_image_directory_"+oid,root,timeout=24),root)
        except Exception:
            listing=[]
        if not any(u.endswith(".img.gz") for u in listing):
            import xml.etree.ElementTree as ET
            bucket="https://nasa-heasarc.s3.amazonaws.com/"
            prefix=f"swift/data/obs/{month}/{oid}/uvot/image/"
            xml=fetch("swift_uvot_s3_files_"+oid,bucket,params={"list-type":"2","prefix":prefix},timeout=25)
            listing=[bucket+el.text for el in ET.fromstring(xml).findall(".//{*}Key") if el.text]
            root=bucket+prefix
        total=0
        for filt in ("uw2","um2","uw1"):
            imgs=[u for u in listing if basename(u).lower()==f"sw{oid}{filt}_sk.img.gz"]
            exp=[u for u in listing if basename(u).lower()==f"sw{oid}{filt}_ex.img.gz"]
            if imgs and exp:
                ib=fetch("swift_sky_"+oid+filt,imgs[0],timeout=65)
                eb=fetch("swift_exp_"+oid+filt,exp[0],timeout=65)
                with fits.open(io.BytesIO(gzip.decompress(ib)),memmap=False) as ih, fits.open(io.BytesIO(gzip.decompress(eb)),memmap=False) as eh:
                    total+=image_measure(ih,eh,root,{"uw2":"UVW2","um2":"UVM2","uw1":"UVW1"}[filt],oid)
        return total
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs={pool.submit(visit,x):x[0] for x in visits.items()}
        for future in as_completed(jobs):
            try:
                n+=future.result()
            except Exception as e:
                LOG.append(("Swift UVOT "+jobs[future],"file access/measurement failed",str(e)[:300]))
    record("Swift UVOT images","measured" if n else ("coverage_only" if visits else "no_match"),
           f"{len(rr)} Swift master fields with centres within 12.6 arcmin, {len(visits)} unique IDs, {n} image apertures; local annulus and valid exposure required.")


def galex_time(h):
    date=str(h.get("DATE-OBS","")).strip()
    if date and "T" not in date and h.get("TIME-OBS"):
        date+="T"+str(h["TIME-OBS"]).strip()
    for v in (h.get("MJD-OBS"),h.get("MJD_OBS"),date):
        t=epoch(v)
        if t:
            return t
    for k in ("EXPSTART","GPSSTART","TSTART"):
        try:
            val=float(h[k])
            if 6e8<val<1.3e9 and "GPS" in str(h.comments[k]).upper():
                return Time(val,format="gps").utc.mjd
        except (KeyError,ValueError,TypeError):
            continue
    return ""


def ga_products():
    urls=[r["final_url"] for r in common.receipts if str(r.get("key","")).startswith("galex_mcat_fits_") and "final_url" in r]
    out=[]
    for i,u in enumerate(dict.fromkeys(urls)):
        raw=fetch("galex_visit_reparse_"+str(i),u,timeout=45)
        with fits.open(io.BytesIO(gzip.decompress(raw))) as hd:
            for tab in hd:
                if not isinstance(tab,fits.BinTableHDU):
                    continue
                keys={k.lower():k for k in tab.data.names}
                if i==0:
                    LOG.append(("GALEX MCAT", "column names", ",".join(keys)))
                best=None
                for a,d in (("alpha_j2000","delta_j2000"),("ra","dec"),("alpha_j2000_nuv","delta_j2000_nuv"),("nuv_alpha_j2000","nuv_delta_j2000"),("ra_nuv","dec_nuv")):
                    if a not in keys or d not in keys:
                        continue
                    ra=np.array(tab.data[keys[a]],float)
                    dec=np.array(tab.data[keys[d]],float)
                    dd=3600*np.hypot((ra-HOST_RA)*math.cos(math.radians(HOST_DEC)),dec-HOST_DEC)
                    if np.any(np.isfinite(dd)):
                        idx=int(np.nanargmin(dd))
                        if best is None or dd[idx]<best[0]:
                            best=dd[idx],idx
                if best is None or best[0]>4:
                    LOG.append(("GALEX MCAT",basename(u),"No centroid within 4 arcsec" if best else "Unrecognised coordinate columns"))
                    continue
                r={k:str(tab.data[best[1]][v]) for k,v in keys.items() if np.ndim(tab.data[best[1]][v])==0}
                t=galex_time(tab.header) or galex_time(hd[0].header)
                coadd="/01-main/" in u
                for label,short in (("nuv","nd"),("fuv","fd")):
                    m=num(r,label+"_mag",label+"_mag_auto")
                    e=num(r,label+"_magerr",label+"_magerr_auto",label+"_mag_err")
                    flag=num(r,label+"_artifact",label+"_artifact_flag")
                    out.append({"file":basename(u),"band":label,"separation_arcsec":best[0],"mag":m,"error":e,"artifact":flag,"fields":str({k:v for k,v in r.items() if any(c in k for c in ("mag","flag","artifact","date","epoch"))})})
                    if not (0<m<35 and 0<e<.36) or (math.isfinite(flag) and flag>0):
                        continue
                    tm=t
                    if not tm and not coadd:
                        try:
                            cu=u.replace("-xd-mcat.fits.gz","-"+short+"-cnt.fits.gz")
                            b=bytes_range("galex_time_"+str(i)+label,cu,0,65535)
                            h0=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(b,57600)
                            en=next(y for y in range(0,len(h0),80) if h0[y:y+8]==b"END     ")
                            tm=galex_time(fits.Header.fromstring(h0[:en+80].decode(),sep=""))
                        except Exception as ex:
                            LOG.append(("GALEX visit date",basename(u),str(ex)[:200]))
                    bounds=[v for v in common.events if v["archive"]=="GALEX footprint" and v["band"]==label.upper()]
                    end=""
                    if coadd and bounds:
                        tm,end=bounds[0]["mjd"],bounds[0]["end_mjd"]
                    event("GALEX "+("coadd" if coadd else "visit"),label.upper(),tm,m,e,"mag AB","detection",u,
                          note=f"MCAT matched at {best[0]:.2f} arcsec; artifact={flag}, 4 arcsec match. Overlapping coadd/visits.",end_mjd=end)
    with (OUT/"galex-matched-bands.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=["file","band","separation_arcsec","mag","error","artifact","fields"])
        w.writeheader()
        w.writerows(out)
    n=sum(r["archive"] in ("GALEX coadd","GALEX visit") for r in common.events)
    record("GALEX matched products","measured" if n else "no_usable_points",f"{len(set(urls))} MCAT catalogs, {n} retained visit/coadd bands. Band limits require explicit sensitivity; a blank is unconstrained.")


def aperture_wise():
    # Independent fixed 8.25 arcsec (aperture 2) measurement, suitable for the extended host.
    allrows=[]
    for table in ("allsky_4band_p1bs_psd","allsky_3band_p1bs_psd","postcryo_p1bs_psd","neowiser_p1bs_psd"):
        try:
            q=(f"SELECT * FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{HOST_RA},{HOST_DEC},0.001111111))=1")
            rr=tap("aperture_"+table,IRSA,q)
            flags=Counter()
            for r in rr:
                t=epoch(col(r,"mjd"))
                if not t or not 55197<=t<60676 or any(num(r,k)<=0 for k in ("qual_frame","qi_fact")):
                    continue
                for i in (1,2):
                    field="w"+str(i)
                    cc=col(r,"cc_flags")[i-1:i]
                    moon=col(r,"moon_masked")[i-1:i]
                    fl=num(r,field+"flg_2")
                    if cc!="0" or moon not in ("","0") or fl!=0:
                        flags[f"{field}:cc={cc},moon={moon},apflag={fl}"]+=1
                        continue
                    m,e=num(r,field+"mag_2"),num(r,field+"sigm_2")
                    if 0<m<30 and 0<e<.36:
                        event("WISE fixed aperture",field.upper(),t,m,e,"mag Vega",source=table+"/"+col(r,"cntr"),note="8.25 arcsec aperture, cc_flags=0, aperture flag=0. Includes host and any blended source.")
                        allrows.append((field.upper(),t,m))
            LOG.append((table,"aperture selection",f"{len(rr)} rows, common rejection reasons: {dict(flags.most_common(8))}"))
        except Exception as ex:
            LOG.append((table,"aperture query failed",str(ex)[:220]))
    with (OUT/"wise-fixed-aperture.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["band","first_mjd","last_mjd","N","median_mag","mad_mag"])
        for band in ("W1","W2"):
            series=[(t,m) for b,t,m in allrows if b==band]
            for k in sorted({math.floor((t-55000)/180) for t,m in series}):
                sub=[(t,m) for t,m in series if math.floor((t-55000)/180)==k]
                ts,ms=zip(*sub)
                med=float(np.median(ms))
                w.writerow([band,min(ts),max(ts),len(sub),med,np.median(np.abs(np.array(ms)-med))])
    record("WISE fixed aperture","measured" if allrows else "no_usable_points",
           f"{len(allrows)} quality-filtered W1/W2 single-exposure aperture measurements; bands, sampling and scatter tabulated in wise-fixed-aperture.csv.")


def crts_retry():
    from photometry import rows_csv
    p={"RA":HOST_RA,"Dec":HOST_DEC,"Rad":0.05,"DB":"photcat","OUT":"csv"}
    for i,host in enumerate(("http://nesssi.cacr.caltech.edu","http://nunuku.caltech.edu","https://crts.iucaa.in","http://crts.iucaa.in")):
        try:
            b=fetch("crts_public_mirror_"+str(i),host+"/cgi-bin/getcssconedb_release_img.cgi",params=p,timeout=15)
            ss=BeautifulSoup(b,"html.parser").get_text("\n") if b"<" in b[:256] else b.decode(errors="replace")
            lines=[x.strip() for x in ss.splitlines() if "," in x]
            ix=next(j for j,x in enumerate(lines) if "mjd" in x.lower() and "mag" in x.lower())
            rr=rows_csv(("\n".join(lines[ix:])).encode())
            n=0
            for r in rr:
                r={k.strip():v for k,v in r.items()}
                t,m,e=epoch(col(r,"mjd")),num(r,"mag","v"),num(r,"magerr","mag_err")
                if t and 0<m<30 and 0<e<.36 and off(num(r,"ra"),num(r,"dec"))<=3:
                    event("CRTS","CV",t,m,e,"mag CRTS",source=col(r,"css_id","id","masterid"),note="Catalina unfiltered V calibration")
                    n+=1
            return record("CRTS mirrors","measured" if n else "no_usable_points",f"{n} matched photometry rows at {host}")
        except Exception as exc:
            LOG.append(("CRTS",host,str(exc)[:220]))
    record("CRTS mirrors","unavailable","Original and mirror cone-service endpoints did not supply parseable photometry.")


def run():
    jobs={"VLASS NRAO images":vlass,"Swift UVOT images":swift,"FIRST image":first,"GALEX matched products":ga_products,"WISE fixed aperture":aperture_wise,"CRTS mirrors":crts_retry}
    with ThreadPoolExecutor(max_workers=len(jobs)) as pool:
        futures={pool.submit(fn):a for a,fn in jobs.items()}
        for fut in as_completed(futures):
            try:
                fut.result()
            except Exception as e:
                record(futures[fut],"unavailable",e)
    with (OUT/"image-checks.csv").open("w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["archive","status","detail"])
        w.writerows(sorted(LOG))
