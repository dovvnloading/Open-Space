"""Position-specific public archive follow-up in the observed 0.2–2.0 keV band."""
import csv
import gzip
import hashlib
import importlib.util
import io
import json
import math
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse_single_table
from astropy.time import Time

from common import DEC, HOST_DEC, HOST_RA, OUT, RA, event, fetch, status

BASE = Path("results/xmmsl2-erass1-rass")
PREV = BASE / "outputs"
NASA = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
CSC = "https://cda.cfa.harvard.edu/csc21tap/sync"
ESA = "https://sky.esa.int/esasky-tap/tap/sync"
ENERGY = "0.2–2.0 keV"
UNIT = "erg cm^-2 s^-1"
ledger, hashes, notes = [], [], []
target = SkyCoord(HOST_RA, HOST_DEC, unit="deg")


def number(v):
    try:
        v = float(v)
        return v if math.isfinite(v) else None
    except (TypeError, ValueError):
        return None


def mjd(v):
    if v is None or str(v).strip() == "":
        return ""
    n = number(v)
    if n is not None:
        return n if 39000 < n < 90000 else ""
    try:
        return float(Time(str(v).strip().replace(" ", "T").rstrip("Z")).mjd)
    except ValueError:
        return ""


def clean(v):
    if np.ma.is_masked(v):
        return None
    if isinstance(v, bytes):
        return v.decode("utf-8", "replace")
    if isinstance(v, np.generic):
        return v.item()
    return v


def first(row, *keys):
    return next((row[k.lower()] for k in keys if row.get(k.lower()) not in (None, "")), None)


def scale(low, high):
    return model.band(.2, 2., 2., 3e20) / model.band(low, high, 2., 3e20)


def native(archive, ident, start, end, f, err, bounds, method, kind="detection", factor=None):
    f, err = number(f), number(err)
    if f is None or f <= 0:
        return
    mult = scale(*bounds) if factor is None else factor
    projected = f * mult
    e = abs(err) * mult if err is not None else ""
    ledger.append(dict(archive=archive, identifier=ident, start_mjd=start, end_mjd=end,
                       kind=kind, native_lo=bounds[0], native_hi=bounds[1], native_flux=f,
                       native_err=err, common_flux=projected, common_err=e, factor=mult, method=method))
    event(archive, ENERGY, start, projected, e, UNIT, kind, str(ident), method, end_mjd=end)


def input_file(path):
    raw = path.read_bytes()
    hashes.append((str(path), hashlib.sha256(raw).hexdigest()))
    return raw


def archived(key, data):
    """Keep small derived inputs from the UKSSDC client with deterministic hashes."""
    data = data if isinstance(data, bytes) else json.dumps(data, default=str).encode()
    out = OUT / "xray-inputs"
    out.mkdir(exist_ok=True, parents=True)
    (out / (key+".gz")).write_bytes(gzip.compress(data, mtime=0))
    hashes.append(("xray-inputs/"+key+".gz [uncompressed bytes]", hashlib.sha256(data).hexdigest()))


def tap(tag, url, query, limit=1000):
    raw = fetch("xr-"+tag, url, params=dict(REQUEST="doQuery", LANG="ADQL",
                FORMAT="votable", MAXREC=limit, QUERY=query), timeout=30)
    if re.search(rb'name=["\']QUERY_STATUS["\'][^>]*value=["\']ERROR', raw, re.I) or \
            re.search(rb'value=["\']ERROR["\'][^>]*name=["\']QUERY_STATUS', raw, re.I):
        raise ValueError("TAP query failed: "+re.sub(r"<[^>]+>", " ", raw.decode(errors="replace"))[-450:])
    if re.search(rb'value=["\']OVERFLOW', raw, re.I):
        raise ValueError("TAP row limit reached")
    t = parse_single_table(io.BytesIO(raw)).to_table(use_names_over_ids=True)
    return [{k.lower():clean(r[k]) for k in t.colnames} for r in t]


def cone(ra, dec, rad, x=HOST_RA, y=HOST_DEC):
    return (f"1=CONTAINS(POINT('ICRS',{ra},{dec}),"
            f"CIRCLE('ICRS',{x},{y},{rad/3600:.9f}))")


def sep(r, ra="ra", dec="dec"):
    return target.separation(SkyCoord(float(r[ra]), float(r[dec]), unit="deg")).arcsec


def guard(label, fn):
    try:
        fn()
    except Exception as exc:
        status(label, "unavailable", f"{type(exc).__name__}: {str(exc)[:380]}")


def baseline():
    global model
    spec = importlib.util.spec_from_file_location("slew_calibration", BASE/"catalogs.py")
    input_file(BASE/"catalogs.py")
    model = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(model)
    rows = csv.DictReader(io.StringIO(input_file(PREV/"candidates.csv").decode()))
    r = next(r for r in rows if r["name"] == "XMMSL2 J125851.2+235522")
    # The catalogue lacks a source-specific RASS date in the comparison. The dates delimit the main survey.
    native("ROSAT RASS",r["rass_name"],mjd("1990-07-30"),mjd("1991-01-25"),
           r["rass_flux"],r["rass_err"],(.2,2),"Main-survey date window only. PR 218 PSPC response; ROSAT beam includes neighbours.",factor=1)
    native("XMM slew",r["xmm_obs"],r["xmm_mjd_min"],r["xmm_mjd_max"],r["xmm_flux"],r["xmm_err"],
           (.2,2),"Released XMMSL2 band 6.",factor=1)
    native("eROSITA eRASS1",r["ero_name"],r["ero_mjd_min"],r["ero_mjd_max"],r["ero_flux"],r["ero_err"],
           (.2,2),"Source-specific catalogue time window; PR 218 model-projected from 0.2–2.3 keV.",factor=1)
    s = json.loads(input_file(PREV/"rass_hardness/swift2sxps_review.json"))
    f = {k:d["value"] for k,d in s["measurement_fields"].items()}
    native("Swift XRT 2SXPS",s["returned_name"],f["time"],f["end_time"],f["pow_flux_a"],
           f["pow_flux_a_pos_err"],(.3,10),"Observed catalogue mean over two visits; overlaps individual XRT light-curve rows.")
    notes.append(f"Swift 2SXPS gives HR1 {f['hardness_ratio_1']:+.2f} "
                 f"(+{f['hardness_ratio_1_pos_err']:.2f}/{f['hardness_ratio_1_neg_err']:.2f}) "
                 f"and HR2 {f['hardness_ratio_2']:+.2f} "
                 f"(+{f['hardness_ratio_2_pos_err']:.2f}/{f['hardness_ratio_2_neg_err']:.2f}). "
                 "The two visits are pooled. Its band edges differ from ROSAT's.")
    s = json.loads(input_file(PREV/"paper_context/dr2hc.json"))
    v = next(r["values"] for r in s["rows"] if r["values"]["UIdCSC"]=="2CXO J125851.4+235526")
    native("eROSITA 3eRASS",v["IAUName"],"","",v["MLFlux0"],v["e_MLFlux0"],(.2,5),
           "Later hard catalogue, unique detection. Source-specific date absent in this cross-match.")
    status("Survey baseline","measured","PR 218 detections, stacked Swift and later eROSITA linked by identifier.")


def footprint():
    """Point-in-polygon for actual footprints from the ESA multi-mission TAP."""
    tables = tap("esa-observations",ESA,"SELECT table_name FROM TAP_SCHEMA.tables WHERE schema_name='observations'")
    for rec in tables:
        table = rec["table_name"]
        tag = table.split(".")[-1]
        if not any(k in tag.lower() for k in ("chandra","xmm","swift_xrt","einstein_probe")):
            continue
        def one():
            cols = {r["column_name"].lower() for r in tap(tag+"-schema",ESA,
                    f"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='{table}'")}
            if "s_region" not in cols:
                raise ValueError("Archive does not expose an s_region footprint")
            wanted = [c for c in ("observation_id","obs_id","obsid","observation_oid","instrument_name",
                      "observation_date","start_time","end_time","t_min","t_max","duration","public") if c in cols]
            if not wanted:
                raise ValueError("No observation identifier or epoch columns recognized")
            rows = tap(tag,ESA,f"SELECT {','.join(wanted)} FROM {table} WHERE "
                       f"1=CONTAINS(POINT('ICRS',{HOST_RA},{HOST_DEC}),s_region)")
            for r in rows:
                start = mjd(first(r,"t_min","start_time","observation_date"))
                end = mjd(first(r,"t_max","end_time"))
                obs = first(r,"observation_id","obs_id","obsid","observation_oid")
                event("ESA "+tag,"X-ray",start,kind="footprint",source=str(obs),end_mjd=end,
                      note="Point in archived observation footprint; no photometry implied.")
            status("ESA "+tag,"coverage_only" if rows else "no_match",f"{len(rows)} point-in-footprint rows")
        guard("ESA "+tag,one)


def nasa_obs(table, arcmin):
    # Conservative centre search; these are possible visits until an on-detector
    # footprint, exposure map or a matched source row establishes coverage.
    rows = tap(table,NASA,f"SELECT * FROM heasarc.{table} WHERE "+cone("ra","dec",60*arcmin))
    for r in rows:
        ident = first(r,"obsid","obs_id","sequence_number","name")
        event(table,"X-ray coverage",mjd(first(r,"time","start_time")),kind="pointing candidate",
              source=str(ident),end_mjd=mjd(first(r,"end_time","stop_time")),
              note=f"Boresight {sep(r)/60:.2f} arcmin from nucleus; chip coverage/flux unverified.")
    status(table,"coverage_only" if rows else "no_match",
           f"{len(rows)} boresights within {arcmin} arcmin; pointing candidates only.")


def xmm():
    # The HEASARC rolling XMMSSC table contains individual source detections.
    found = False
    for table in ("xmmssc", "xmm4xmm"):
        try:
            rows = tap(table,NASA,f"SELECT * FROM heasarc.{table} WHERE "+cone("ra","dec",12))
        except Exception as exc:
            status("HEASARC "+table,"unavailable",str(exc)[:280])
            continue
        for r in rows:
            if sep(r)>6:
                continue
            bands = [number(first(r,f"ep_{i}_flux",f"ep{i}_flux",f"flux_band_{i}",f"flux_b{i}")) for i in (1,2,3)]
            errs = [number(first(r,f"ep_{i}_flux_err",f"ep{i}_flux_err",f"flux_band_{i}_error",f"flux_b{i}_error")) for i in (1,2,3)]
            if all(v is not None for v in bands):
                native("XMM pointed",first(r,"obsid","obs_id","name"),mjd(first(r,"time","mjd_start")),
                       mjd(first(r,"end_time","mjd_stop")),sum(bands),math.hypot(*errs) if all(v is not None for v in errs) else None,
                       (.2,2),"EPIC catalogue bands 1+2+3; errors in quadrature. Detection flags in raw table.",factor=1)
                found=True
        status("XMM pointed catalogue","measured" if found else ("no_match" if not rows else "coverage_only"),
               f"{len(rows)} detections within 12 arcsec. Only sources within 6 arcsec of the optical nucleus projected.")
        return


def csc():
    for idx, url in enumerate((CSC,"https://cda.cfa.harvard.edu/csc2tap/sync")):
        try:
            tables = tap("csc-tables"+str(idx),url,"SELECT table_name FROM TAP_SCHEMA.tables")
            keys = [r["table_name"] for r in tables if "source" in r["table_name"].lower()]
            if not keys:
                raise ValueError("CSC metadata has no source tables")
            break
        except Exception as exc:
            status("Chandra CSC service "+str(idx+1),"unavailable",str(exc)[:260])
    else:
        return
    linkage = []
    def order(table):
        return 0 if table.lower().endswith(".stack_source") else 1
    for t in sorted(keys, key=order):
        short=t.lower().split(".")[-1]
        if short not in ("master_source","stack_source","observation_source","obs_source","per_observation"):
            continue
        def one():
            columns=tap("csc-"+short+"-columns",url,
                        f"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='{t}'",5000)
            cols={r["column_name"].lower() for r in columns}
            wanted=[k for k in ("name","obsid","obi","region_id","detect_stack_id","ra","dec","gti_start","gti_end",
                                "gti_mjd_obs","gti_mjd_end","mjd_obs","mjd_end","date_obs","date_end","instrument",
                                "flux_aper_b","flux_aper_lolim_b","flux_aper_hilim_b",
                                "flux_aper_s","flux_aper_m","flux_aper_h","flux_aper_u",
                                "flux_aper_lolim_s","flux_aper_hilim_s","flux_aper_lolim_m","flux_aper_hilim_m",
                                "sat_src_flag","pileup_flag","conf_flag") if k in cols]
            if {"ra","dec"}<=cols:
                constraint=cone("ra","dec",8)
            elif {"detect_stack_id","region_id"}<=cols and linkage:
                conditions=[f"(detect_stack_id='{sid}' AND region_id={rid})" for sid,rid,dist in linkage]
                constraint=" OR ".join(conditions)
            else:
                raise ValueError("No independent coordinates or compatible stack/region key")
            rows=tap("csc-"+short,url,f"SELECT {','.join(wanted)} FROM {t} WHERE "+constraint)
            count=0
            for r in rows:
                if r.get("ra") is not None and r.get("dec") is not None:
                    dist=sep(r)
                else:
                    dist=next(d for s,i,d in linkage if s==r["detect_stack_id"] and i==int(r["region_id"]))
                if dist>2.5:
                    continue
                if short=="stack_source" and r.get("detect_stack_id") and r.get("region_id") is not None:
                    sid=str(r["detect_stack_id"])
                    if re.fullmatch(r"[\w.+:-]+",sid):
                        linkage.append((sid,int(r["region_id"]),dist))
                start=mjd(first(r,"gti_mjd_obs","mjd_obs","date_obs","gti_start"))
                end=mjd(first(r,"gti_mjd_end","mjd_end","date_end","gti_end"))
                ident=first(r,"obsid","name","detect_stack_id") or "CSC"
                label="Chandra CSC "+short
                # Use observed aperture energy fluxes. HRC broad/wide differs from ACIS broad.
                soft, medium = number(r.get("flux_aper_s")), number(r.get("flux_aper_m"))
                if all(v is not None and v>0 for v in (soft,medium)):
                    er=[number(r.get("flux_aper_"+edge+"_"+band)) for edge in ("hilim","lolim") for band in ("s","m")]
                    err=max(er[0]+er[1]-soft-medium,soft+medium-er[2]-er[3]) if None not in er else None
                    native(label,ident,start,end,soft+medium,err,(.5,2),
                           f"CSC ACIS soft+medium; conservative sum of interval widths; {dist:.2f} arcsec. "
                           "Master/stack rows can share photons.")
                    count+=1
                    continue
                for b,bounds in (("s",(.5,1.2)),("b",(.5,7))):
                    f=number(first(r,"flux_aper_"+b))
                    if f is None or f<=0:
                        continue
                    lo=number(first(r,"flux_aper_lolim_"+b))
                    hi=number(first(r,"flux_aper_hilim_"+b))
                    err=max(hi-f,f-lo) if lo is not None and hi is not None and lo<=f<=hi else None
                    native(label,ident,start,end,f,err,bounds,f"CSC aperture band {b}; {dist:.2f} arcsec from nucleus. "
                           f"Native interval [{lo},{hi}]; common error uses larger side. Master/stack rows can share photons.")
                    count+=1
                    break
            status("Chandra "+short,"measured" if count else ("coverage_only" if rows else "no_match"),
                   f"{len(rows)} nearby rows, {count} converted within 2.5 arcsec of the nucleus.")
        guard("Chandra "+short,one)


def chandra_footprint():
    raw=fetch("xr-chandra-footprint","https://cda.cfa.harvard.edu/footprint/get_vo_table.pl",
              params=dict(pos=f"{HOST_RA},{HOST_DEC}",size="0.0001"),timeout=25)
    tab=parse_single_table(io.BytesIO(raw)).to_table(use_names_over_ids=True)
    for entry in tab:
        row={k.lower():clean(entry[k]) for k in tab.colnames}
        ident=first(row,"obsid","obs_id","observation_id","obs_id_int")
        start=mjd(first(row,"obs_start","date_obs","start_date","t_min"))
        end=mjd(first(row,"obs_end","date_end","end_date","t_max"))
        event("Chandra footprint","X-ray",start,kind="footprint",source=str(ident or ""),
              note="Footprint service at optical nucleus. No exposure-map aperture limit implied.",end_mjd=end)
    status("Chandra footprint","coverage_only" if len(tab) else "no_match",f"{len(tab)} intersecting footprints.")


def swift():
    from swifttools.ukssdc.data import SXPS
    from swifttools.ukssdc.query import SXPSQuery
    from swifttools.ukssdc import main as uk
    import requests

    # The upstream SDK handles its catalogue-specific API. Record its exact HTTP
    # response bodies and impose a finite timeout without patching global requests.
    class ObservedClient:
        def __init__(self):
            self.counter = 0

        def __getattr__(self, name):
            return getattr(requests, name)

        def post(self, url, **kwargs):
            kwargs.setdefault("timeout", 25)
            response = requests.post(url, **kwargs)
            self.counter += 1
            archived("swift-api-"+str(self.counter),response.content)
            req=kwargs.get("json",kwargs.get("data",{}))
            archived("swift-query-"+str(self.counter),dict(url=url,parameters=req))
            return response

    if hasattr(uk,"requests"):
        uk.requests = ObservedClient()
    # The 2SXPS measurement identifies the earlier source and supplies the observed
    # total-band ECF independently of any SDK-specific field naming.
    prior = json.loads((PREV/"rass_hardness/swift2sxps_review.json").read_text())
    ecf = prior["measurement_fields"]["powcf_flux_a"]["value"]
    def lightcurve(cat, key):
        product = SXPS.getLightCurves(**key, cat=cat, bands="total", binning="observation",
                   timeFormat="MJD", returnData=True, saveData=False, groupRates=False, groupULs=False)
        if not isinstance(product,dict):
            raise ValueError("Unexpected UKSSDC product")
        datasets=[]
        def walk(v, path=""):
            if hasattr(v,"to_dict") and hasattr(v,"columns"):
                datasets.append((path,v))
            elif isinstance(v,dict):
                for k, item in v.items():
                    walk(item,path+"/"+str(k))
        walk(product)
        obs=set()
        for dataset, df in datasets:
            if "total" not in dataset.lower() and len(datasets)>1:
                continue
            archived("swift-"+cat+"-"+re.sub(r"[^a-zA-Z0-9]+","-",dataset),df.to_csv(index=False).encode())
            for original in df.to_dict("records"):
                r={re.sub(r"[^a-z0-9]","",str(k).lower()):v for k,v in original.items()}
                t=mjd(first(r,"mjd","time"))
                if t=="":
                    continue
                uplimit=number(first(r,"upperlimit","ratelimit","ul","upperlimitrate"))
                rate=number(first(r,"rate","countrate"))
                if rate is not None and "ul" in dataset.lower() and uplimit is None:
                    uplimit=rate
                is_ul=uplimit is not None and uplimit>0
                v=uplimit if is_ul else rate
                if v is None or v<=0:
                    continue
                # API bins span visits. Time uncertainties are offsets, not absolute dates.
                lower=number(first(r,"tneg","mjdneg","timeneg","timeerrneg"))
                upper=number(first(r,"tpos","mjdpos","timepos","timeerrpos"))
                start=t-abs(lower) if lower is not None else t
                end=t+abs(upper) if upper is not None else t
                ident=first(r,"obsid","obsids","observationid") or key.get("sourceName",str(key))
                pos=number(first(r,"ratepos","rateerrpos","error","rateerr"))
                neg=number(first(r,"rateneg","rateerrneg"))
                err=max((abs(x) for x in (pos,neg) if x is not None),default=None)
                native("Swift XRT "+cat,ident,start,end,v*ecf,None if is_ul or err is None else err*ecf,
                       (.3,10),f"{dataset} observation bin, total-band catalogue rate × 2SXPS observed ECF {ecf:.3g}; "
                       "then projected to the common band. Overlapping releases are repeated estimates.",
                       kind="upper limit" if is_ul else "detection")
                obs.add(str(ident)+str(t))
        if not obs:
            raise ValueError("No usable total-band bins; returned product keys "+str(list(product.keys()))[:200])
        status("Swift XRT "+cat,"measured",f"{len(obs)} dated total-band observation bins with rate or limit")
    guard("Swift XRT 2SXPS light curve",lambda:lightcurve("2SXPS",{"sourceName":prior["returned_name"]}))
    def live():
        q=SXPSQuery(cat="LSXPS")
        q.table="sources"
        q.addConeSearch(ra=HOST_RA,dec=HOST_DEC,radius=12,units="arcsec")
        q.submit()
        df=q.results
        archived("swift-LSXPS-cone",df.to_csv(index=False).encode())
        ids=[]
        for r in df.to_dict("records"):
            r={str(k).lower():v for k,v in r.items()}
            if not first(r,"lsxps_id","sourceid","source_id","iau_name"):
                continue
            if number(first(r,"ra")) is not None and number(first(r,"decl","dec")) is not None:
                d=target.separation(SkyCoord(float(r["ra"]),float(first(r,"decl","dec")),unit="deg")).arcsec
                if d>10:
                    continue
            sid=first(r,"lsxps_id","sourceid","source_id")
            ids.append({"sourceID":int(sid)} if sid is not None else {"sourceName":r["iau_name"]})
        if len(ids)==1:
            lightcurve("LSXPS",ids[0])
        else:
            status("Swift XRT LSXPS","no_match" if not ids else "coverage_only",
                   f"Found {len(ids)} identifications in 12 arcsec; require a unique source for photometry.")
    guard("Swift XRT LSXPS",live)


def ep():
    # Query national archive metadata first, then the IVOA registry for any
    # independently registered SCS archive. Save the service capabilities.
    refs=tap("ep-heasarc-schema",NASA,
             "SELECT table_name,description FROM TAP_SCHEMA.tables WHERE "
             "description LIKE '%Einstein Probe%' OR description LIKE '%Einstein-Probe%' "
             "OR table_name LIKE 'heasarc.epwxt%' OR table_name LIKE 'heasarc.epfxt%'")
    usable=[]
    for r in refs:
        t=r["table_name"]
        if not re.fullmatch(r"[A-Za-z0-9_.]+",t):
            continue
        col=tap("ep-col-"+t.split(".")[-1],NASA,
                f"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='{t}'")
        keys={a["column_name"].lower() for a in col}
        if {"ra","dec"}<=keys:
            selection=tap("ep-pos-"+t.split(".")[-1],NASA,
                          f"SELECT * FROM {t} WHERE "+cone("ra","dec",180))
            usable.append((t,len(selection)))
            for a in selection:
                event("Einstein Probe "+t,"X-ray",mjd(first(a,"time","mjd")),kind="archive position",
                      source=str(first(a,"name","obsid","trigger_id")),note=f"{sep(a):.1f} arcsec; association unverified.")
    status("Einstein Probe HEASARC","coverage_only" if usable else "unavailable",
           f"Mission tables with sky positions and counts in 3 arcmin: {usable}. Source/alert-table census only.")
    registry=tap("ep-vo-registry","https://dc.g-vo.org/tap/sync",
        "SELECT DISTINCT res_title,standard_id,access_url FROM rr.resource NATURAL JOIN rr.capability "
        "NATURAL JOIN rr.interface WHERE 1=ivo_nocasematch(res_title,'%Einstein Probe%')")
    endpoints=[x for x in registry if "conesearch" in str(x["standard_id"]).lower()]
    for i,r in enumerate(endpoints):
        raw=fetch("xr-ep-reg-"+str(i),r["access_url"],params=dict(RA=HOST_RA,DEC=HOST_DEC,SR=.05))
        n=len(parse_single_table(io.BytesIO(raw)).array)
        status("EP "+str(r["res_title"]),"coverage_only" if n else "no_match",f"{n} entries within 3 arcmin.")
    if not endpoints:
        status("Einstein Probe VO registry","unavailable",
               f"{len(registry)} registered capabilities, no matching cone service. No WXT/FXT exposure limit established.")


def ep_portal():
    from bs4 import BeautifulSoup
    url="https://ep.bao.ac.cn/ep/"
    soup=BeautifulSoup(fetch("xr-ep-public-portal",url,timeout=22),"html.parser")
    links=[(a.get_text(" ",strip=True),a.get("href","")) for a in soup.select("a[href]")]
    links=[(a,h) for a,h in links if re.search("observation|archive|public data|data center|upper limit", a+" "+h,re.I)]
    archived("ep-portal-links",links)
    status("Einstein Probe Science Center","coverage_only" if links else "unavailable",
           f"{len(links)} archive/navigation links. Portal availability alone supplies no position-specific exposure or upper limit.")


def run():
    baseline()
    tasks={"Chandra CSC":csc,"Chandra footprint":chandra_footprint,"XMM source catalogues":xmm,"ESA observing footprints":footprint,
           "Swift archive products":swift,"Einstein Probe observations":ep,"EP science center":ep_portal}
    for t,a in (("chanmaster",22),("xmmmaster",18),("swiftmastr",17),("rosmaster",60)):
        tasks[t]=lambda t=t,a=a:nasa_obs(t,a)
    with ThreadPoolExecutor(max_workers=7) as pool:
        list(pool.map(lambda pair:guard(*pair),tasks.items()))
    with (OUT/"xray-fluxes.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(ledger[0]))
        w.writeheader()
        w.writerows(sorted(ledger,key=lambda r:(mjd(r["start_mjd"]) or 1e7,r["archive"],str(r["identifier"]))))
    with (OUT/"xray-extra-inputs.tsv").open("w") as f:
        f.write("input\tsha256\n")
        for path,digest in sorted(hashes):
            f.write(f"{path}\t{digest}\n")
    text=["# X-ray epochs and sensitivity","",
          "ROSAT and eRASS1 establish the long-baseline soft X-ray change. "
          "The XMM slew is intermediate. The catalogue interval for this target's eRASS1 scans falls in December 2019. "
          "The Swift catalogue mean covers two sparse visits. Later eROSITA emission establishes continued activity.","",
          "The table reports observed 0.2–2.0 keV energy flux in erg cm⁻² s⁻¹. "
          "A band projection uses the previous survey's absorbed power law with photon index 2 and column 3 × 10²⁰ cm⁻². "
          "It multiplies the released observed energy flux by the ratio of absorbed energy integrals. "
          "Instrument response and spectral assumptions of the native catalogue persist. Statistical errors exclude model uncertainty. "
          "The X-ray-only measurements cannot separate luminosity and absorption changes.","",
          "| Archive | ID | First UTC | Last UTC | Common flux | Statistical error | Type |",
          "|---|---|---|---|---:|---:|---|"]
    def utc(v):
        return Time(float(v),format="mjd").utc.isot[:10] if mjd(v)!="" else "unreported"
    for r in sorted(ledger,key=lambda r:(mjd(r["start_mjd"]) or 1e7,r["archive"])):
        err=f"{r['common_err']:.3g}" if r['common_err']!="" else "—"
        text.append(f"| {r['archive']} | {r['identifier']} | {utc(r['start_mjd'])} | {utc(r['end_mjd'])} | "
                    f"{r['common_flux']:.3g} | {err} | {r['kind']} |")
    text+=["",*notes,"","A boresight match is a possible exposure. Detector gaps, field edges, observing modes and background affect "
           "sensitivity. An absent catalogue row or an unavailable service supplies no calibrated upper limit. "
           "Pointing candidates and footprints are in the timeline. Their flux cells are intentionally empty. "
           "The RASS dates delimit its main survey. CSC master/stack and repeated Swift products can share events.",
           "","[Measurement ledger](xray-fluxes.csv) · [Extra input fingerprints](xray-extra-inputs.tsv)"]
    (OUT/"xray.md").write_text("\n".join(text)+"\n")
