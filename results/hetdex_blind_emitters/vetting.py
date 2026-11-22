"""Post-selection tests for HETDEX features. All remote requests run inside Actions."""
import csv
import gzip
import hashlib
import heapq
import io
import json
import math
import os
from collections import Counter, defaultdict
from pathlib import Path
import threading

import numpy as np
import requests
from astropy.io import fits
from scipy.spatial import cKDTree
from vet_cubes import PDR, MPE, blank_cube, cube_test, mad
from photometry import imaging_ew

CAT_URL = "https://web.corral.tacc.utexas.edu/hetdex/HETDEX/catalogs/"
LAN = CAT_URL+"HETDEX_LAN/hetdex_lan_v0.3.fits"
DESI = CAT_URL+"HETDEX_DESI/DESI_HETDEX_SPEC_v1.6.fits"
DESI_BACKUP = "https://data.desi.lbl.gov/public/dr1/vac/dr1/hetdex/v1.6/DESI_HETDEX_SPEC_v1.6.fits"
# Table 1, Mirza Khanlari et al. 2026, doi:10.3847/1538-4357/aea6d9.
PRIOR = [(3005020033,165.19334,50.792854,2.213),(3007998098,221.16264,53.845905,1.994),
         (3009751099,19.451492,.471328,2.038),(3010547400,20.290672,-1.048341,2.053),
         (4019081078,233.91087,51.616085,1.994),(4028177445,228.20595,50.607033,2.184),
         (4028632963,211.51462,52.452236,2.106),(5001528748,170.90740,49.509910,2.204)]
REMOTE_PER_CHANNEL = 32


def atom(v):
    if isinstance(v,dict):
        return {str(k):atom(x) for k,x in v.items()}
    if isinstance(v,(list,tuple,np.ndarray)):
        return [atom(x) for x in v]
    if isinstance(v,(np.generic, bytes)):
        if isinstance(v, bytes):
            return v.decode("utf-8","replace").strip()
        return atom(v.item())
    if isinstance(v,float) and not math.isfinite(v):
        return None
    return v


def number(v, default=-math.inf):
    try:
        v = float(v)
        return v if np.isfinite(v) and v != -999 else default
    except (TypeError,ValueError):
        return default


def chord(ra,dec):
    a,d = np.deg2rad(ra),np.deg2rad(dec)
    return np.stack((np.cos(d)*np.cos(a),np.cos(d)*np.sin(a),np.sin(d)),axis=-1)


def separation(ra,dec,r2,d2):
    dist = np.linalg.norm(chord(ra,dec)-chord(r2,d2),axis=-1)
    return np.rad2deg(2*np.arcsin(np.minimum(dist/2,1)))*3600


class Cache:
    def __init__(self, base):
        self.root = Path(base)/"vet-cache"
        self.root.mkdir(parents=True,exist_ok=True)
        self.inputs = {}
        self.lock = threading.RLock()

    def fetch(self, kind, name, urls):
        with self.lock:
            return self._fetch(kind,name,urls)

    def _fetch(self, kind, name, urls):
        p = self.root/name
        meta = p.with_name(p.name+".json")
        if p.exists() and meta.exists():
            self.inputs[name] = json.loads(meta.read_text())
            p.touch()
            return p
        if kind == "cube":
            old = sorted(self.root.glob("dex_cube*.fits"),key=lambda q:q.stat().st_mtime)
            for retired in old[:-6]:
                retired.unlink(missing_ok=True)
        errors = []
        for url in urls:
            tmp = p.with_name(p.name+".part")
            try:
                digest = hashlib.sha256()
                with requests.get(url, timeout=(8,45), stream=True) as r, tmp.open("wb") as out:
                    r.raise_for_status()
                    for b in r.iter_content(1024*1024):
                        if b:
                            digest.update(b)
                            out.write(b)
                with tmp.open("rb") as inp:
                    if inp.read(8)!=b"SIMPLE  ":
                        raise ValueError("response is not a FITS file")
                tmp.replace(p)
                self.inputs[name] = dict(url=url,bytes=p.stat().st_size,sha256=digest.hexdigest())
                meta.write_text(json.dumps(self.inputs[name]))
                return p
            except (requests.RequestException,OSError,ValueError) as exc:
                errors.append(str(exc))
                tmp.unlink(missing_ok=True)
        self.inputs[name] = dict(urls=urls,error=errors)
        raise OSError("; ".join(errors))


def blank_apertures(catalog, cache_dir, n_cubes=8):
    cache = Cache(cache_dir)
    t = catalog.info
    det = catalog.detinfo
    rng = np.random.default_rng(31597)
    # Fixed random source-observation order; no filtering on class, wavelength or flux.
    indices = rng.permutation(len(t))
    seen = set()
    made = 0
    for i in indices:
        if made==n_cubes or len(seen)>=4*n_cubes:
            return
        m = catalog.metadata(int(i))
        try:
            shot,slot = int(m["shotid"]),str(int(m["ifuslot"])).zfill(3)
        except (KeyError,TypeError,ValueError):
            continue
        if (shot,slot) in seen:
            continue
        seen.add((shot,slot))
        name = f"dex_cube_{shot}_{slot}.fits"
        meta = dict(shotid=shot,ifuslot=slot,cube=name,coordinates=[])
        try:
            file = cache.fetch("cube",name,[f"{p}{shot}/{name}" for p in (PDR,MPE)])
            if det is None:
                raise ValueError("no detection table for source exclusion")
            sel = np.flatnonzero(det["shotid"]==shot)
            coords = np.stack((det["ra_det"][sel],det["dec_det"][sel]),axis=1)
            coords = coords[np.isfinite(coords).all(axis=1) & (np.abs(coords[:,1])<=90)]
            coords = coords[separation(coords[:,0],coords[:,1],m["ra"],m["dec"])<100]
            if not len(coords):
                raise ValueError("no catalogued detections in chosen pointing")
            centers,fs,es = blank_cube(file,coords,halo_meta=m,trial_output=meta)
            meta.update(cache.inputs[name],coordinates=centers,n_local_detections=len(coords))
            if len(fs):
                made += 1
            yield meta, fs, es
        except Exception as exc:
            meta.update(cache.inputs.get(name,{}),error=str(exc))
            yield meta,np.empty((0,len(catalog.wave))),np.empty((0,len(catalog.wave)))


def spatial(catalog):
    """Join independently located, common-wavelength knots; test diffuse emission later."""
    t,det = catalog.info,catalog.detinfo
    if det is None:
        return []
    ra,dec = det["ra_det"],det["dec_det"]
    wave,sn = det["wave"],det["sn"]
    valid = ((sn>=7) & (wave>=3500) & (wave<=5500) & np.isfinite(ra) &
             (np.abs(dec)<=90) & (ra>=0) & (ra<360))
    idx = np.flatnonzero(valid)
    shots = det["shotid"][idx]
    idx = idx[np.argsort(shots,kind="stable")]
    cuts = np.flatnonzero(np.diff(det["shotid"][idx]))+1
    sortedids = np.argsort(t["source_id"],kind="stable")
    ids = t["source_id"][sortedids]
    results=[]
    for group in np.split(idx,cuts):
        if len(group)<2:
            continue
        xyz = chord(ra[group],dec[group])
        pairs = cKDTree(xyz).query_pairs(np.deg2rad(14/3600),output_type="ndarray")
        if not len(pairs):
            continue
        a,b = group[pairs].T
        seps = separation(ra[a],dec[a],ra[b],dec[b])
        keep = (seps>=4.) & (np.abs(wave[a]-wave[b])<=5.)
        # One row per physical pair in a shot; neighbors are retained even if the
        # original spatial grouping assigned identical source IDs to both knots.
        for u,v,s in zip(a[keep],b[keep],seps[keep]):
            if sn[v]>sn[u]:
                u,v=v,u
            key=det["source_id"][u]
            l,r = np.searchsorted(ids,key,"left"), np.searchsorted(ids,key,"right")
            obs=sortedids[l:r]
            obs=obs[t["shotid"][obs]==det["shotid"][u]]
            if not len(obs):
                continue
            ref=int(obs[0])
            m=catalog.metadata(ref)
            if not m.get("ifuslot"):
                continue
            rec=dict(m,catalog=m.copy(),channel="blob",row_index=ref,
                     ra=float(ra[u]),dec=float(dec[u]),ra_det=float(ra[u]),dec_det=float(dec[u]),
                     knot_detectid=int(det["detectid"][u]),wave=float(wave[u]),snr=float(sn[u]),
                     z_lya=float(wave[u]*1.00028/1215.67-1),flux=float(det["flux"][u])*1e-17,
                     flux_err=float(det["flux_err"][u])*1e-17,
                     secondary_detectid=int(det["detectid"][v]),secondary_wave=float(wave[v]),
                     secondary_ra=float(ra[v]),secondary_dec=float(dec[v]),
                     secondary_snr=float(sn[v]),knot_separation_arcsec=float(s))
            results.append(rec)
    # Single-center arm is independent of optical size and LAN membership.
    # Draw up to 256 distinct source observations per field above catalog S/N 9.
    # Order is seeded independently of Ly-alpha halo size and published morphology.
    source_types=np.char.lower(np.char.strip(t["source_type"].astype("S10")))
    eligible=np.flatnonzero((source_types==b"lae") & (t["sn"]>=9) &
                            (t["z_hetdex"]>=1.88) & (t["z_hetdex"]<=3.52))
    rng=np.random.default_rng(91841)
    for field in np.unique(t["field"][eligible]):
        rows=eligible[t["field"][eligible]==field]
        for i in rng.permutation(rows)[:256]:
            m=catalog.metadata(int(i))
            wv=number(m.get("wave"))
            if not (3500<=wv<=5500 and m.get("ifuslot")):
                continue
            results.append(dict(m,catalog=m.copy(),channel="blob",row_index=int(i),
                      morphology="central",wave=wv,snr=number(m.get("sn")),
                      z_lya=float(wv*1.00028/1215.67-1),
                      flux=number(m.get("flux"),0.)*1e-17,
                      flux_err=number(m.get("flux_err"),0.)*1e-17))
    return results


class Ancillary:
    def __init__(self,cache):
        self.tables={}
        self.files=[]
        self.errors={}
        for key, urls, rac, decc in (("LAN",[LAN],"ra","dec"),
                                     ("DESI",[DESI,DESI_BACKUP],"ra_hetdex","dec_hetdex")):
            try:
                p=cache.fetch("catalog",urls[0].rsplit("/",1)[-1],urls)
                h=fits.open(p,memmap=True,character_as_bytes=True)
                self.files.append(h)
                tab=h["INFO"].data if "INFO" in h else h[1].data
                coords=chord(tab[rac],tab[decc])
                good=np.isfinite(coords).all(axis=1) & (np.abs(tab[decc])<=90)
                self.tables[key]=(h,tab,np.flatnonzero(good),cKDTree(coords[good]))
            except Exception as exc:
                self.errors[key]=str(exc)

    def lookup(self, key,ra,dec,radius):
        if key not in self.tables:
            return None,None
        _,tab,idx,tree = self.tables[key]
        dist, j = tree.query(chord(ra,dec))
        sep=float(np.rad2deg(2*np.arcsin(min(dist/2,1)))*3600)
        if sep>radius:
            return None,None
        i=int(idx[j])
        row={s.lower():atom(tab[s][i]) for s in tab.names if np.ndim(tab[s][i])==0}
        row["sep_arcsec"]=sep
        return i,row

    def close(self):
        for h in self.files:
            h.close()


def repeat_line(catalog,i,rec,rows=None):
    measurements=[]
    cen=rec["wave"]
    sigma = np.clip(number(rec.get("sigma_A"),2.5),1.9,8.)
    w=np.asarray(catalog.wave)
    core=np.abs(w-cen)<=3*sigma
    side=(np.abs(w-cen)>4*sigma)&(np.abs(w-cen)<max(30.,7*sigma))
    for j in (catalog.repeat_rows(i) if rows is None else rows):
        j=int(j)
        f,e=catalog.read(j,j+1)
        f,e=np.asarray(f).ravel(),np.asarray(e).ravel()
        ok=np.isfinite(f) & np.isfinite(e) & (e>0)
        use=(core|side)&ok
        if (core&ok).sum()<3 or (side&ok).sum()<6:
            continue
        g=np.exp(-0.5*((w[use]-cen)/sigma)**2)/(np.sqrt(2*np.pi)*sigma)
        X=np.stack((g,np.ones(len(g))),axis=1)/e[use,None]
        try:
            cov=np.linalg.inv(X.T@X)
        except np.linalg.LinAlgError:
            continue
        F=float((cov@X.T@(f[use]/e[use]))[0])
        err=float(np.sqrt(cov[0,0]))
        # Uncertainties may be underestimated in correlated blank pixels.
        scatter=mad((f[side&ok]-np.median(f[side&ok]))/e[side&ok])
        err *= max(1.,scatter) if np.isfinite(scatter) else 1.
        measurements.append(dict(row_index=j,shotid=int(catalog.info["shotid"][j]),
                                 flux=F*1e-17,flux_err=err*1e-17,snr=F/err))
    return measurements


def tap_request(service, items):
    endpoints = {
        "SIMBAD":("https://simbad.cds.unistra.fr/simbad/sim-tap/sync",
                  "main_id,otype,ra,dec,rvz_redshift,rvz_bibcode,nbref","basic","ICRS"),
        "NED":("https://ned.ipac.caltech.edu/tap/sync","*","objdir","J2000"),
        "Legacy":("https://datalab.noirlab.edu/tap/sync",
                  "ls_id,ra,dec,type,flux_g,flux_r,flux_z,flux_ivar_r,maskbits,nobs_r,fracmasked_r,fracflux_r,mw_transmission_r,"
                  "flux_w1,flux_w2,flux_ivar_w1,flux_ivar_w2,fracflux_w1,fracflux_w2,wise_coadd_id","ls_dr10.tractor","ICRS")}
    endpoints["SIMBAD mirror"]=("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",
                               *endpoints["SIMBAD"][1:])
    endpoints["NED mirror"]=("https://vo.ned.ipac.caltech.edu/tap/sync",
                            *endpoints["NED"][1:])
    endpoints["Legacy DR9"]=(endpoints["Legacy"][0],endpoints["Legacy"][1],"ls_dr9.tractor","ICRS")
    url,cols,table,frame=endpoints[service]
    clauses=[f"1=CONTAINS(POINT('{frame}',ra,dec),CIRCLE('{frame}',{r['ra']:.8f},{r['dec']:.8f},0.00083333334))"
             for r in items]
    query=f"SELECT {cols} FROM {table} WHERE "+" OR ".join(clauses)
    params=dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",QUERY=query,MAXREC="2000")
    try:
        res=requests.get(url,params=params,timeout=(8,32))
        res.raise_for_status()
        reader=csv.DictReader(io.StringIO(res.content.decode("utf-8-sig")))
        names={s.strip().lower() for s in (reader.fieldnames or [])}
        if not {"ra","dec"}<=names:
            raise ValueError("unexpected TAP result: "+res.text[:150])
        results=[{k.strip().lower():v for k,v in row.items() if k} for row in reader]
        if len(results)>=2000:
            raise ValueError("TAP result limit reached")
        answer=[]
        for rec in items:
            inside=[]
            for row in results:
                if not (math.isfinite(number(row["ra"])) and math.isfinite(number(row["dec"]))):
                    continue
                sep=float(separation(rec["ra"],rec["dec"],float(row["ra"]),float(row["dec"])))
                if sep<3.0001:
                    inside.append(dict(row,sep_arcsec=sep))
            answer.append(dict(status="queried",service=service,table=table,matches=sorted(inside,key=lambda r:r["sep_arcsec"]),
                               response_sha256=hashlib.sha256(res.content).hexdigest(),query_url=res.url))
        return answer
    except Exception as exc:
        return [dict(status="unavailable",error=str(exc),endpoint=url,query=query) for _ in items]


def tap_batch(service,items):
    primary=tap_request(service,items)
    alternate={"SIMBAD":"SIMBAD mirror","NED":"NED mirror","Legacy":"Legacy DR9"}[service]
    retry=[i for i,r in enumerate(primary) if r["status"]!="queried" or
           service=="Legacy" and not r.get("matches")]
    if retry:
        fallbacks=tap_request(alternate,[items[i] for i in retry])
        for i,other in zip(retry,fallbacks):
            if other["status"]=="queried" or primary[i]["status"]!="queried":
                other["first_attempt"]=primary[i]
                primary[i]=other
            else:
                primary[i]["second_attempt"]=other
    return primary


def hsc_photometry(rec,layer="hsc-dr2"):
    # HSC DR2 publicly redisplayed by the Legacy Surveys viewer. No credentials.
    params=dict(layer=layer,bands="r",ra=rec["ra"],dec=rec["dec"],
                size=96,pixscale=0.3)
    r=requests.get("https://www.legacysurvey.org/viewer/fits-cutout",params=params,timeout=(8,22))
    r.raise_for_status()
    digest=hashlib.sha256(r.content).hexdigest()
    with fits.open(io.BytesIO(r.content),memmap=False) as h:
        plane=next((x for x in h if x.data is not None and np.ndim(x.data) in (2,3)),None)
        if plane is None:
            raise ValueError("no HSC image plane")
        data=np.squeeze(plane.data)
        if data.ndim!=2:
            raise ValueError("unexpected HSC image shape")
        valid=np.isfinite(data)&(data!=0)
        cy,cx=(np.array(data.shape)-1)/2
        y,x=np.mgrid[:data.shape[0],:data.shape[1]]
        ap=np.hypot(x-cx,y-cy)*0.3<1.5
        if not np.all(valid[ap]):
            return dict(status="no central coverage",layer=layer,url=r.url,sha256=digest)
        bg=(np.hypot(x-cx,y-cy)*.3>6)&valid
        if bg.sum()<100:
            return dict(status="insufficient background",layer=layer,url=r.url,sha256=digest)
        sky=np.median(data[bg])
        blank=[]
        for ix in range(6,data.shape[1]-6,11):
            for iy in range(6,data.shape[0]-6,11):
                if np.hypot(ix-cx,iy-cy)*.3<=4.5:
                    continue
                b=np.hypot(x-ix,y-iy)*.3<1.5
                if np.all(valid[b]):
                    blank.append(float((data[b]-sky).sum())/b.sum()*ap.sum())
        noise=mad(blank)
        value=float((data[ap]-sky).sum())
        unit=str(plane.header.get("BUNIT","unreported"))
        return dict(status="measured",layer=layer,url=r.url,sha256=digest,unit=unit,flux=value,
                    ra=float(rec["ra"]),dec=float(rec["dec"]),fluxmag0=plane.header.get("FLUXMAG0"),
                    flux_err=float(noise) if np.isfinite(noise) else None,
                    zeropoint=plane.header.get("MAGZP",plane.header.get("MAGZERO")),
                    psfsize=plane.header.get("PSFSIZE"),
                    snr=float(value/noise) if noise>0 else None,aperture_radius_arcsec=1.5,
                    blank_apertures=len(blank))
    

def lowz_explanations(rec):
    # Only demand known correlated lines well away from the trigger wavelength.
    out=[]
    for name,hyp in (rec.get("low_z") or {}).items():
        lines=hyp.get("companions",{})
        hits=[k for k,d in lines.items() if number(d.get("snr"))>=4.5
              and abs(number(d.get("air_wave"))-rec["wave"])>12]
        if hits:
            out.append(dict(primary=name,z=hyp.get("z"),companions=hits))
    return out


def validate_local(rec, cat, aux):
    c=rec.get("catalog",rec)
    ch=rec["channel"]
    ra,dec=rec["ra"],rec["dec"]
    z=rec.get("z_heii") if ch=="heii" else rec.get("z_lya")
    result=dict(disposition="unresolved",reasons=[],passed=[],channel=ch)
    reject=result["reasons"]
    if rec.get("single_pixel") or rec.get("sky_near"):
        reject.append("spectral artifact or masked sky window")
    cz=number(c.get("z_hetdex"))
    st=str(c.get("source_type","")).lower()
    if st=="agn" or number(c.get("agn_flag"))==1 or str(c.get("cls_diagnose","")).upper()=="QSO":
        reject.append("HETDEX AGN classification")
    source=str(c.get("z_hetdex_src","")).lower()
    if st in ("oii","lzg","star") and ch!="unknown":
        independent_or_template=any(tag in source for tag in ("diagnos","sdss","desi","specz","2em"))
        if independent_or_template and number(c.get("z_hetdex_conf"))>=.8:
            reject.append("secure template/multiline HETDEX foreground or stellar classification")
        else:
            result["catalog_conflict"]="single-line foreground classification; retained as ambiguous"
    if ch=="heii" and number(c.get("agn_flag"))==0:
        reject.append("HETDEX broad-line AGN candidate")
    result["low_z_families"]=lowz_explanations(rec)
    if ch!="unknown" and any(len({k.rstrip("12") if k.startswith("oii") else k for k in h["companions"]})>=2 for h in result["low_z_families"]):
        reject.append("low-redshift pattern with at least two companion lines")
    if (not reject and ch=="heii" and rec.get("heii_family") and
        any(number(rec["heii_family"].get(k,{}).get("snr"))>=3 for k in
            ("civ1","civ2","ciii1","ciii2","oiii_uv1","oiii_uv2"))):
        reject.append("UV metal emission")
    if reject:
        result["disposition"]="excluded"
        return result
    il,lan=aux.lookup("LAN",ra,dec,3.0)
    if lan and abs(number(lan.get("z_hetdex"))-number(z))<.015:
        result["LAN"]={k:v for k,v in lan.items() if k in
            ("name","ra","dec","detectid","shotid","z_hetdex","area_r_iso_circ","flag_resolved","r_iso","r_s","flux_lya",
             "flux_lya_err","hsc-r_mag","hsc-r_mag_err","combined_eqw_rest_lya","sep_arcsec","source_type")}
        try:
            result["LAN"]["flux_lya_unit"]=aux.tables["LAN"][1].columns["flux_lya"].unit
        except (KeyError,AttributeError):
            pass
        if ch=="blob" and number(lan.get("flag_resolved"))==1:
            result.update(disposition="published",reasons=["HETDEX LAN already resolved"])
        if ch=="ew" and number(lan.get("combined_eqw_rest_lya"))>240:
            result.update(disposition="published",reasons=["HETDEX LAN already gives EW >240 A"])
        if str(lan.get("source_type","")).lower()=="agn":
            result.update(disposition="excluded",reasons=["HETDEX LAN AGN classification"])
    if ch=="heii":
        for ident,r,d,z0 in PRIOR:
            if int(c["detectid"])==ident or (separation(ra,dec,r,d)<=3. and abs(number(z)-z0)<.015):
                result.update(disposition="published",reasons=["published He II sample, ApJ 1010, 62"],
                              prior_heii_detectid=ident)
    if result["disposition"]!="unresolved":
        return result
    di,desi=aux.lookup("DESI",ra,dec,1.25)
    if desi:
        result["DESI"]={k:desi[k] for k in ("targetid","vi_z","vi_quality","source_type","detectid",
                          "sn_hetdex","wave_desi","flux_desi","flux_err_desi","sep_arcsec") if k in desi}
        fiber_sep=float(separation(ra,dec,number(desi.get("target_ra")),number(desi.get("target_dec"))))
        result["DESI"].update(target_sep_arcsec=fiber_sep,dex_flag=desi.get("dex_flag"))
        association=number(desi.get("dex_flag"))==1 and fiber_sep<=1.25
        if association and number(desi.get("vi_quality"))>=3:
            if ch!="unknown" and abs(number(desi.get("vi_z"))-number(z))>0.02 and number(desi.get("vi_z"))>=0:
                result.update(disposition="excluded",reasons=["independent DESI classification/redshift"])
                return result
            if abs(number(desi.get("vi_z"))-number(z))<0.008:
                result["passed"].append("independent DESI visual redshift")
        if association and ch in ("ew","heii","blob","unknown"):
            from desi_lines import desi_lines
            try:
                result["DESI_lines"]=desi_lines(aux.tables["DESI"][0],di,rec["wave"],ch)
            except Exception as exc:
                result["DESI_lines"]={"status":"unavailable","error":str(exc)}
    result["repeats"]=([r for r in rec["repeat_measurements"] if int(r["row_index"])!=int(rec["row_index"])]
                       if rec.get("repeat_measurements") else repeat_line(cat,int(rec["row_index"]),rec))
    if any(number(x.get("snr"))>=4 for x in result["repeats"]):
        result["passed"].append("independent HETDEX shot")
    if result["repeats"] and ch in ("heii","ew","blob"):
        pair=rec.get("lya_velocity_best") if ch=="heii" else max(
            (v for k,v in (rec.get("lya_family") or {}).items() if k in
             ("heii","oiii_uv1","oiii_uv2","ciii1","ciii2")),key=lambda d:number(d.get("snr")),default=None)
        if pair and number(pair.get("snr"))>=5 and number(pair.get("air_wave"))>0:
            probe={**rec,"wave":pair["air_wave"]}
            if rec.get("stream")=="repeat_stack" and rec.get("repeat_rows"):
                visits=repeat_line(cat,int(rec["row_index"]),probe,rows=rec["repeat_rows"])
                result["secondary_visit_measurements"]=visits
                result["secondary_repeats"]=[r for r in visits if r["shotid"]!=rec["shotid"]]
            else:
                result["secondary_repeats"]=repeat_line(cat,int(rec["row_index"]),probe)
    if st=="lae" and abs(cz-number(z)) < .012:
        result["passed"].append("HETDEX Ly-alpha redshift agrees")
    return result


def fixed_external_pairs(rec,result):
    """Both UV lines must recur at fixed centers in the same external fibre spectrum."""
    if rec["channel"] not in ("heii","ew","blob"):
        return []
    z=number(rec.get("z_heii" if rec["channel"]=="heii" else "z_lya"))
    candidates=[]
    vac=result.get("DESI") or {}
    if number(vac.get("dex_flag"))==1 and number(vac.get("target_sep_arcsec"),math.inf)<=1.25:
        candidates.append(("DESI_HETDEX_VAC",vac.get("targetid"),result.get("DESI_lines") or {}))
    for m in (result.get("SPARCL") or {}).get("matches",[]):
        if number(m.get("sep_arcsec"),math.inf)<=1.25 and m.get("diagnostics"):
            candidates.append((m.get("data_release"),m.get("spectrum_key") or m.get("sparcl_id") or m.get("specid"),
                               m["diagnostics"]))
    out=[]
    for source,ident,d in candidates:
        if abs(number(d.get("z_topocentric"),-99)-z)>.0003:
            continue
        lines=d.get("lines") or {}
        l,h=lines.get("lya",{}),lines.get("heii",{})
        if (l.get("coverage")=="ok" and h.get("coverage")=="ok" and
            number(l.get("snr"))>=5 and number(h.get("snr"))>=4):
            out.append(dict(source=source,identifier=ident,lya_snr=l["snr"],heii_snr=h["snr"],
                            lya_wave_vac=l.get("wave_vac"),heii_wave_vac=h.get("wave_vac")))
    return out


def stack_pair_shots(rec,result):
    """Distinct individual visits with primary ≥3σ and secondary ≥4σ in the same extraction."""
    primary=rec.get("repeat_measurements") or result.get("repeats") or []
    secondary=result.get("secondary_visit_measurements",result.get("secondary_repeats",[]))
    def detections(rows,threshold):
        return {(r["row_index"],r["shotid"]) for r in rows if r.get("row_index") is not None and
                r.get("shotid") is not None and number(r.get("snr"))>=threshold}
    return sorted({shot for _,shot in detections(primary,3) & detections(secondary,4)})


def ordinary_agn_disfavored(rec,result,z):
    ch=rec["channel"]
    ext=result.get("SPARCL") or {}
    galaxy=any(str(m.get("spectype","")).upper()=="GALAXY" and number(m.get("redshift_warning"))==0
               and abs(number(m.get("redshift"))-z)<.006 for m in ext.get("matches",[]))
    fam=rec.get("heii_family" if ch=="heii" else "lya_family") or {}
    narrow=2.355*number(rec.get("sigma_A"),math.inf)/rec["wave"]*299792.<900
    limits=[number(fam.get(k,{}).get("limit3"),math.inf)/rec["flux"] for k in ("nv1","civ1")]
    return narrow and (max(limits)<(.3 if ch=="heii" else .15) or (galaxy and min(limits)<.15))


PAIR_EXTERNAL="fixed Ly-alpha and He II recur in one independent instrument spectrum"
PAIR_REPEAT="both source lines recur in independent HETDEX shots"
AGN_PASSED="ordinary broad/UV-line AGN disfavored; obscured or weak-line AGN remains possible"
AGN_GAP="ordinary AGN insufficiently constrained by width and UV limits"


def external_check(rec,result,cache):
    from desi_lines import external_lines
    from concurrent.futures import ThreadPoolExecutor
    def safe(fn,*args):
        try:
            return fn(*args)
        except Exception as exc:
            return dict(status="unavailable",error=str(exc))
    with ThreadPoolExecutor(max_workers=2) as pool:
        supplied=rec.get("rescue_audit") or {}
        f=pool.submit(lambda:supplied.get("HSC") or safe(hsc_photometry,rec))
        ext=supplied.get("SPARCL") or safe(external_lines,rec["ra"],rec["dec"],rec["wave"],rec["channel"],cache.root)
        hsc=f.result()
    result.update(SPARCL=ext,HSC=hsc)
    if supplied.get("Legacy_forced"):
        result["Legacy_forced"]=supplied["Legacy_forced"]
    elif rec.get("channel")=="ew" and hsc.get("status")!="measured":
        result["Legacy_forced"]=safe(hsc_photometry,rec,"ls-dr9")
    m=rec.get("catalog",{})
    if "hsc" in str(m.get("counterpart_catalog_name","")).lower():
        dist=number(m.get("counterpart_dist"))
        mag=number(m.get("counterpart_mag"))
        err=number(m.get("counterpart_mag_err"))
        result["HETDEX_HSC"]=dict(survey=m.get("counterpart_catalog_name"),band=m.get("counterpart_filter_name"),
            ab_mag=mag if 0<mag<30 else None,ab_mag_err=err if err>0 else None,
            sep_arcsec=dist if dist>=0 else None,link="HETDEX representative-detection imaging extraction",
            measurement="published AB magnitude; aperture may differ from the spectral extraction")
        if 0<mag<30 and 0<=dist<=1.25:
            result["HETDEX_HSC"]["flux_nanomaggy"]=10**(.4*(22.5-mag))
    if hsc.get("status") in ("no central coverage","insufficient background"):
        hsc["footprint"]="no usable coverage at candidate position"
    elif hsc.get("status")=="unavailable":
        hsc["footprint"]="undetermined after request failure"
    refs=result["passed"]
    ch=rec["channel"]
    z=number(rec.get("z_heii" if ch=="heii" else "z_lya"))
    hard, gaps=[],[]
    mid=result["Legacy"].get("matches",[])
    if mid:
        nearest=mid[0]
        f1,f2=[number(nearest.get("flux_"+k)) for k in ("w1","w2")]
        i1,i2=[number(nearest.get("flux_ivar_"+k)) for k in ("w1","w2")]
        if min(f1,f2,i1,i2)>0:
            s1,s2=f1*math.sqrt(i1),f2*math.sqrt(i2)
            m1,m2=22.5-2.5*math.log10(f1)-2.699,22.5-2.5*math.log10(f2)-3.339
            dcolor=(2.5/math.log(10))*math.hypot(1/s1,1/s2)
            isolated=(len(mid)==1 and nearest["sep_arcsec"]<=1.25 and
                      number(nearest.get("maskbits"))==0 and
                      all(0<=number(nearest.get("fracflux_"+k))<.1 for k in ("w1","w2")))
            flag=isolated and min(s1,s2)>=5 and m2<=15 and m1-m2-2*dcolor>=.8
            result["WISE"]=dict(ls_id=nearest.get("ls_id"),wise_coadd_id=nearest.get("wise_coadd_id"),
                        w1_minus_w2_vega=m1-m2,color_err=dcolor,w2_vega=m2,
                        w1_snr=s1,w2_snr=s2,passes_isolation=isolated,possible_agn=flag,
                        source="Stern et al. 2012 doi:10.1088/0004-637X/753/1/30; Tractor forced WISE AB flux")
            if flag:
                gaps.append("isolated bright red WISE counterpart supports competing AGN interpretation")
    for service in ("NED","SIMBAD"):
        cone=result[service]
        if cone["status"]!="queried":
            gaps.append(service+" search unavailable")
        for m in cone.get("matches",[]):
            zz=number(m.get("z",m.get("rvz_redshift")))
            if (ch!="unknown" and m["sep_arcsec"]<=1.25 and zz>=0 and abs(zz-z)>.02 and
                    (m.get("rvz_bibcode") or m.get("refcode"))):
                gaps.append(service+" published counterpart redshift disagrees; association needs review")
            ty=str(m.get("pretype",m.get("otype",""))).strip().upper()
            if m["sep_arcsec"]<=1.25 and ty in ("QSO","AGN","SY1","SY2","SEYFERT","BLL","BLAZAR","BLA"):
                hard.append(service+" known active nucleus")
            # Published redshifts may be photometric or from HETDEX itself;
            # they inform the ledger but do not certify independent spectroscopy.
    external_redshift="independent DESI visual redshift" in refs
    external_primary=False
    metal_tests={}
    if ch=="heii":
        # HETDEX covers C III] only at the lowest redshifts. A missing component
        # cannot certify an ion limit. Summing component bounds is conservative
        # for emission; nominal individual fits to a blend are correlated.
        family=rec.get("heii_family") or {}
        for ion,keys in (("civ",("civ1","civ2")),("oiii",("oiii_uv1","oiii_uv2")),
                         ("ciii",("ciii1","ciii2"))):
            values=[family.get(k,{}) for k in keys]
            if all(d.get("status")=="measured" and d.get("limit3") is not None for d in values):
                ratio=sum(d["limit3"] for d in values)/rec["flux"]
                metal_tests[ion]={"source":"HETDEX","ratio_upper":ratio}
    metal_ok=False
    diagnostics=[]
    if result.get("DESI_lines"):
        diagnostics.append(result["DESI_lines"])
    for m in ext.get("matches",[]):
        if number(m.get("redshift_warning"))==0 and number(m.get("redshift"))>=0:
            if abs(number(m["redshift"])-z)<.006:
                external_redshift=True
            elif ch!="unknown":
                gaps.append("external spectroscopic redshift disagrees")
            if str(m.get("spectype","")).upper() in ("QSO","STAR"):
                hard.append("external spectroscopy classifies a QSO or star")
        if m.get("diagnostics"):
            diagnostics.append(m["diagnostics"])
    for d in diagnostics:
        lines=d.get("lines",{})
        primary=lines.get("primary",{})
        if number(primary.get("snr"))>=4.0:
            external_primary=True
        if ch=="heii" and number(primary.get("snr"))>=4.0:
            upper=[lines.get(k,{}).get("max_upper3sigma") for k in ("civ","oiii","ciii")]
            for ion,limit in zip(("civ","oiii","ciii"),upper):
                if limit is not None and limit>0 and number(primary.get("flux"))>0:
                    ratio=limit/primary["flux"]
                    if ratio<metal_tests.get(ion,{}).get("ratio_upper",math.inf):
                        metal_tests[ion]={"source":"independent spectrum","ratio_upper":ratio}
            if any(number(lines.get(k,{}).get("snr"))>=3 for k in ("civ","oiii","ciii")):
                hard.append("external metal emission")
        hits=[k for k in ("alt_hbeta","alt_halpha","alt_oiii5007","alt_oiii4959","alt_neiii")
              if number(lines.get(k,{}).get("snr"))>=5]
        if hits and number(lines.get("alt_oii",{}).get("snr"))>=4:
            hard.append("external [O II] and low-redshift companion")
    if ch=="heii":
        result["three_ion_limits_relative_to_heii"]=metal_tests
        metal_ok=all(metal_tests.get(k,{}).get("ratio_upper",math.inf)<=0.5 for k in ("civ","oiii","ciii"))
    if external_primary:
        refs.append("independent instrument confirms line")
    if external_redshift and "independent DESI visual redshift" not in refs:
        refs.append("independent survey redshift agrees")
    try:
        result["cube"]=supplied.get("cube") or cube_test(rec,cache.fetch)
    except Exception as exc:
        result["cube"]=dict(status="unavailable",error=str(exc))
    result["imaging_EW"]=imaging_ew(rec,result)
    cube=result["cube"]
    sky=cube.get("background_wavelength_test") or {}
    if cube.get("background_suspect") or abs(number(sky.get("sn"),0))>=5:
        gaps.append("same-wavelength cube sky annulus shows coherent residual")
    if cube.get("central_above_5sigma"):
        refs.append("masked cube central line above 5 empirical sigma")
    else:
        gaps.append("no independently noise-calibrated central cube line")
    if not external_primary and "independent HETDEX shot" not in refs:
        gaps.append("independent-visit line confirmation absent")
    independent_primary=external_primary or "independent HETDEX shot" in refs
    # Re-measure the other member of a common-redshift pair at fixed wavelength in
    # independent visits. Multiple reductions of the same shot never count.
    two_line=False
    other=None
    if ch=="heii":
        ly=rec.get("lya_velocity_best") or {}
        if number(ly.get("snr"))>=5:
            other=ly
    elif ch in ("ew","blob"):
        eligible=[d for k,d in (rec.get("lya_family") or {}).items()
                  if k in ("heii","oiii_uv1","oiii_uv2","ciii1","ciii2") and number(d.get("snr"))>=5]
        if eligible:
            other=max(eligible,key=lambda d:d["snr"])
    if other and result.get("secondary_repeats"):
        two_line=any(number(m.get("snr"))>=4 for m in result["secondary_repeats"])
    if rec.get("stream")=="repeat_stack":
        result["stack_pair_shots"]=stack_pair_shots(rec,result)
        two_line=other is not None and len(result["stack_pair_shots"])>=2
    pairs=fixed_external_pairs(rec,result)
    if pairs:
        result["external_fixed_uv_pairs"]=pairs
        refs.append(PAIR_EXTERNAL)
    id_supported=external_redshift or bool(pairs) or (independent_primary and two_line)
    if not id_supported and ch in ("ew","blob") and rec["wave"]<3714 and independent_primary:
        # Ly-alpha is below all optical OII/Balmer/OIII assignments. Require a
        # second high-z line and the survey classification to discriminate UV IDs.
        id_supported=(other is not None and (rec.get("stream")!="repeat_stack" or two_line)
                      and not result.get("low_z_families") and
                      "HETDEX Ly-alpha redshift agrees" in refs)
    if two_line:
        refs.append(PAIR_REPEAT)
    if ch=="unknown":
        from unknown_lines import reassess_unknown
        audit=reassess_unknown(rec,result,cache)
        result["unassigned_line_audit"]=audit
        stat=str(audit.get("status","")).lower()
        if stat=="common alternative has a supporting companion":
            hard.append("a known transition has independent-wavelength companion support")
        elif stat=="transition coincidence at independent redshift":
            gaps.append("published atomic transition close to this independently anchored wavelength")
        elif stat=="unmatched at independent redshift":
            refs.append("independently anchored host has no matching transition in tested library")
            gaps.append("foreground or projected Ly-alpha source still fits this wavelength")
        else:
            gaps.append("candidate lacks a secure identification across plausible redshifts")
    if ch=="heii":
        if number((rec.get("lya_velocity_best") or {}).get("snr"))<5:
            gaps.append("Ly-alpha counterpart below 5 sigma")
        if not metal_ok:
            gaps.append("joint CIV, OIII] and CIII] limits below half He II unavailable")
        if not id_supported:
            gaps.append("independent line identification unavailable")
    if ch=="ew":
        if rec.get("stream")=="imaging_rescue":
            c=number(rec.get("continuum"),0.)
            ce=number(rec.get("continuum_err"),math.inf)
            if c-3*ce>0 and (rec["flux"]+3*rec["flux_err"])/((c-3*ce)*(1+z))<=240:
                hard.append("spectral continuum excludes extreme EW")
            if not id_supported:
                gaps.append("independent Ly-alpha identification required for imaging-led EW")
        elif number(rec.get("ew0_lower"))<=240:
            gaps.append("spectral EW lower bound does not exceed 240 A")
        if number((result["imaging_EW"] or {}).get("ew0_lower_A"))<=240:
            gaps.append("redward imaging does not independently constrain EW above 240 A")
        if not id_supported:
            gaps.append("external redshift unavailable for low-z and AGN discrimination")
    if ch=="blob":
        mode = "annular" if rec.get("morphology")=="central" else "bridge"
        if not cube.get(mode+"_above_5sigma"):
            gaps.append("extended emission above the fitted seeing profile unconfirmed")
        else:
            refs.append(mode+" PSF excess above 5 empirical sigma")
        if not id_supported:
            gaps.append("independent redshift unavailable")
        # Second-line spatial blends and compact sources require repeat confirmation
        # of the bridge, in addition to an integrated spectral detection.
        repeats=result.get("repeats",[])
        corroborated=False
        for rep in sorted(repeats,key=lambda m:-number(m.get("snr")))[:2]:
            if number(rep.get("snr"))<3 or "repeat_cube_meta" not in rep:
                continue
            alt={**rec,**rep["repeat_cube_meta"]}
            rep["bridge_cube"]=safe(cube_test,alt,cache.fetch)
            corroborated |= bool(rep["bridge_cube"].get(mode+"_above_5sigma"))
        if not corroborated:
            gaps.append(mode+" PSF excess unconfirmed in separate shot")
    if result.get("Legacy",{}).get("status")!="queried" and hsc.get("status")!="measured":
        gaps.append("direct imaging checks unavailable")
    if result.get("LAN_catalog_unavailable") and ch in ("ew","blob"):
        gaps.append("published HETDEX LAN properties unavailable")
    if (result.get("catalog_conflict") or result.get("low_z_families")) and not id_supported:
        gaps.append("foreground interpretation unresolved")
    if not hard and not gaps:
        if ordinary_agn_disfavored(rec,result,z):
            refs.append(AGN_PASSED)
        else:
            gaps.append(AGN_GAP)
    result["pair_evidence_version"]=2
    return dict(result,disposition="excluded" if hard else "unresolved" if gaps else "candidate",
                reasons=sorted(set(hard+gaps)))


def enrich(rec,catalog):
    if rec.get("shuffled_detail_complete"):
        return {}
    from spectral import detail, search
    i=int(rec["row_index"])
    if rec.get("stream")=="repeat_stack" and rec.get("repeat_rows"):
        from repeats import stack_for_rows
        f,e=stack_for_rows(catalog,rec["repeat_rows"])
    else:
        if rec.get('stream')=='shuffled':
            from scan import BATCH
            from spectral import shuffle_fibres
            start=i//BATCH*BATCH
            f,e=catalog.read(start,min(start+BATCH,catalog.size))
            f,e=shuffle_fibres(f,e,410087+start)
            f,e=f[i-start],e[i-start]
        else:
            f,e=catalog.read(i,i+1)
    f,e=np.asarray(f).ravel(),np.asarray(e).ravel()
    if "noise_scale" not in rec or "wave_vac" not in rec:
        # Catalog knots are not necessarily at the center of the representative
        # 1D extraction. A missed peak is kept distinct from a measured line.
        peaks=search(catalog.wave,f,e,density_unit=1e-17)
        same=[p for p in peaks if abs(p["wave"]-rec["wave"])<=5.]
        if not same:
            return {"forced_detail_status":"no blind peak within 5 A in representative extraction"}
        peak=min(same,key=lambda p:abs(p["wave"]-rec["wave"]))
    else:
        peak=rec
    measurements=detail(catalog.wave,f,e,peak,density_unit=1e-17)
    if rec.get("channel")=="unknown":
        from unknown_lines import unknown_spectrum
        measurements["unknown_spectrum"]=unknown_spectrum(catalog.wave,f,e,peak,density_unit=1e-17)
    if "noise_scale" not in rec:
        measurements["sigma_A"]=peak["sigma_A"]
        measurements["spectral_EW0_A"]=peak.get("ew0")
        measurements["spectral_EW0_lower_A"]=peak.get("ew0_lower")
    return measurements


def vet(records,catalog,output_dir):
    """Audit all preliminary rows. Return the auditable status summary."""
    import itertools
    from concurrent.futures import ThreadPoolExecutor
    out=Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    cache=Cache(Path(os.environ.get("RUNNER_TEMP","/tmp"))/"hetdex")
    aux=Ancillary(cache)
    counts=defaultdict(Counter)
    selected=defaultdict(list)
    remote=Counter()
    cubes=Counter()
    total=Counter()
    def write(stream,rec,result):
        ch=rec["channel"]
        counts[ch][result["disposition"]]+=1
        # The full incoming forced-line ledger remains available in the scan products.
        projection={k:v for k,v in rec.items() if k not in rec["catalog"] or k in
            ("ra","dec","ra_det","dec_det","release","shotid","source_id","source_name","detectid","ifuslot",
             "wave","snr","flux","flux_err","source_type","z_hetdex","z_hetdex_src","z_hetdex_conf",
             "agn_flag","flag_best","fiber_id","fwhm","weight","counterpart_mag","counterpart_mag_err",
             "counterpart_catalog_name","counterpart_dist","counterpart_filter_name")}
        for k in ("catalog","low_z","heii_family","lya_family"):
            projection.pop(k,None)
        if rec.get("heii_family"):
            projection["heii_metal_components"]={k:{field:v for field,v in rec["heii_family"][k].items()
                 if field in ("status","flux","flux_err","snr","limit3")} for k in
                 ("civ1","civ2","oiii_uv1","oiii_uv2","ciii1","ciii2")}
        if "SPARCL" in result:
            for key in ("low_z","heii_family","lya_family"):
                if key in rec:
                    projection[key]=rec[key]
        stream.write(json.dumps(atom(dict(peak=projection,vet=result)),allow_nan=False,
                                separators=(",",":"))+"\n")
    try:
        with gzip.open(out/"vetted.jsonl.gz","wt") as stream:
            for idx,rec in enumerate(records):
                ch=rec["channel"]
                total[ch]+=1
                meta=catalog.metadata(int(rec["row_index"]))
                rec={**meta,**rec,"catalog":meta}
                if ch!="blob" and rec.get("ra_det") is not None and rec.get("dec_det") is not None:
                    rec.update(source_group_ra=rec["ra"],source_group_dec=rec["dec"],
                               ra=rec["ra_det"],dec=rec["dec_det"])
                if ch in ("ew","heii") and not rec.get("low_z"):
                    rec.update(enrich(rec,catalog))
                res=validate_local(rec,catalog,aux)
                if aux.errors.get("LAN"):
                    res["LAN_catalog_unavailable"]=True
                if res["disposition"]!="unresolved":
                    write(stream,rec,res)
                    continue
                # Spatial controls must use actual independent shots, not SC1/SC2 reductions of one shot.
                for rep in res.get("repeats",[]):
                    if ch=="blob" and number(rep.get("snr"))>=3:
                        m=catalog.metadata(rep["row_index"])
                        rep["repeat_cube_meta"]={k:m.get(k) for k in ("shotid","ifuslot","fwhm")}
                evidence=("independent HETDEX shot" in res["passed"]) + 2*("DESI" in res)
                if ch=="heii":
                    partner=number((rec.get("lya_velocity_best") or {}).get("snr"),0.0)
                    # The final He II gate requires a real hydrogen companion.
                    # Put objects that can pass it ahead of chance alignments.
                    evidence += 6*int(partner>=5) + 2*int("independent DESI visual redshift" in res["passed"])
                # Among equally corroborated rows, rank by the specific selection observable.
                rank=(number(rec.get("secondary_snr",rec.get("snr"))) if ch=="blob" else
                      number(rec.get("imaging_rank_A",rec.get("ew0_lower"))) if ch=="ew" else number(rec.get("snr")))
                # Avoid spending all external slots on duplicated features of the same object.
                coord=(round(rec["ra"],4),round(rec["dec"],4),round(rec["wave"]/4.))
                arm=(ch+"-"+("central" if rec.get("morphology")=="central" else "pair")) if ch=="blob" else ch
                if ch=="ew" and rec.get("stream")=="imaging_rescue":
                    arm="ew-imaging"
                pool=selected[arm]
                cap=REMOTE_PER_CHANNEL//2 if ch=="blob" else REMOTE_PER_CHANNEL
                duplicate=next((j for j,e in enumerate(pool) if e[2]==coord),None)
                entry=((evidence,rank,-idx),idx,coord,rec,res)
                if duplicate is not None:
                    old=pool[duplicate]
                    if old[0]>=entry[0]:
                        res["reasons"]=["duplicate position and wavelength in ranked follow-up"]
                        write(stream,rec,res)
                        continue
                    pool[duplicate]=entry
                    heapq.heapify(pool)
                elif len(pool)<cap:
                    heapq.heappush(pool,entry)
                    continue
                else:
                    old=heapq.heappushpop(pool,entry)
                old[4]["reasons"]=["outside ranked external follow-up; retained for review"]
                write(stream,old[3],old[4])
            heads=sorted(itertools.chain.from_iterable(selected.values()),key=lambda x:-x[0][0])
            services=("NED","SIMBAD","Legacy")
            for j in range(0,len(heads),12):
                batch=heads[j:j+12]
                rows=[e[3] for e in batch]
                with ThreadPoolExecutor(max_workers=3) as tp:
                    answers=list(tp.map(lambda s:tap_batch(s,rows),services))
                for k,e in enumerate(batch):
                    _,_,_,rec,res=e
                    if rec["channel"] in ("blob","unknown"):
                        extra=enrich(rec,catalog)
                        rec.update(extra)
                        res["low_z_families"]=lowz_explanations(rec)
                    for name,ans in zip(services,answers):
                        res[name]=ans[k]
                def assess(e):
                    rec,res=e[3],e[4]
                    return rec,external_check(rec,res,cache)
                with ThreadPoolExecutor(max_workers=3) as pool:
                    for rec,checked in pool.map(assess,batch):
                        remote[rec["channel"]]+=1
                        if checked.get("cube",{}).get("status")=="measured":
                            cubes[rec["channel"]]+=1
                        write(stream,rec,checked)
                        stream.flush()
        summary=dict(shortlist_counts=dict(total),status_counts={k:dict(v) for k,v in counts.items()},
                     remote_assessed=dict(remote),cube_assessed=dict(cubes),
                     max_remote_per_arm=REMOTE_PER_CHANNEL,
                     spatial_slots_per_arm=REMOTE_PER_CHANNEL//2, independent_imaging_EW_slots=REMOTE_PER_CHANNEL,catalog_errors=aux.errors)
        (out/"vet-summary.json").write_text(json.dumps(summary,indent=2)+"\n")
        (out/"vet-inputs.json").write_text(json.dumps(cache.inputs,indent=2)+"\n")
        return summary
    finally:
        aux.close()
