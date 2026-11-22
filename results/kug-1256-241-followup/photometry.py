"""Public archive photometry at the 2MASX nucleus. Run via the follow-up driver."""
import csv
import gzip
import io
import json
import math
import re
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from statistics import median
from urllib.parse import quote

import numpy as np
from astropy.io import fits, votable
from astropy.time import Time
from bs4 import BeautifulSoup

from common import DEC, HOST_DEC, HOST_RA, OUT, RA, event, fetch, status

POINTS = []
DETAIL = []
AUDIT = []
CENTRE = f"{HOST_RA:.7f} {HOST_DEC:+.7f}"
RADIUS = 3 / 3600
IRSA = "https://irsa.ipac.caltech.edu/TAP/sync"
VIZ = "https://vizier.cds.unistra.fr/viz-bin/votable"
HEASARC = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
MAST = "https://mast.stsci.edu/api/v0/invoke"


def col(row, *names, default=""):
    for n in names:
        x = row.get(n.lower())
        if x is not None and str(x).strip() not in ("", "--", "nan", "None", "null"):
            return str(x).strip()
    return default


def num(r, *keys):
    try:
        return float(col(r, *keys))
    except (ValueError, TypeError):
        return float("nan")


def rows_csv(raw):
    s = raw.decode("utf-8-sig", errors="replace").strip()
    if s.startswith("<") or not s or ("," not in s.splitlines()[0]):
        raise ValueError("expected a CSV table: " + re.sub(r"\s+", " ", s)[:180])
    return [{k.lower(): v for k, v in r.items() if k} for r in csv.DictReader(io.StringIO(s))]


def tap(key, server, q):
    p = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}
    return rows_csv(fetch(key, server, params=p, timeout=65))


def note(name, state, message):
    status(name, state, message)
    DETAIL.append((name, state, message))


def epoch(v):
    if not v or str(v) in ("--", "nan"):
        return ""
    try:
        a = float(v)
        if 30000 < a < 100000:
            return a
        if 2400000 < a < 2500000:
            return a - 2400000.5
    except (TypeError, ValueError):
        pass
    try:
        return float(Time(str(v).strip().replace("/", "-")).mjd)
    except Exception:
        return ""


def off(ra, dec):
    return 3600 * math.hypot((ra - HOST_RA) * math.cos(math.radians(HOST_DEC)), dec - HOST_DEC)


def add(archive, band, t, value, err, unit, src, text="", kind="detection", end=""):
    if not np.isfinite(value):
        return
    e = float(err) if np.isfinite(err) and err > 0 else ""
    event(archive, band, t, round(value, 7), e, unit, kind, source=str(src), note=text, end_mjd=end)
    if t != "" and kind == "detection":
        POINTS.append((archive, band, float(t), value, e, unit))


def wise(table, archive, multi):
    query = (f"SELECT * FROM {table} WHERE CONTAINS(POINT('ICRS',ra,dec),"
             f"CIRCLE('ICRS',{HOST_RA},{HOST_DEC},{RADIUS}))=1")
    rows = tap(table, IRSA, query)
    counts = [0, 0]
    if rows:
        AUDIT.append((archive, "response fields", ",".join(rows[0])))
    for r in rows:
        t = epoch(col(r, "mjd", "mjd_ep"))
        if not t or not 55197 <= t < 60676:
            continue
        if any(num(r, k) <= 0 for k in ("qual_frame", "qi_fact")) or num(r, "saa_sep") < 0:
            continue
        for i in (1, 2):
            w = f"w{i}"
            suf = "_ep" if multi else ""
            flags = col(r, "cc_flags", "cc_flags_ep").upper()
            qual = col(r, "ph_qual", "ph_qual_ep").upper()
            moon = col(r, "moon_masked")
            if (flags[i-1:i] not in ("", "0") or qual[i-1:i] == "U" or
                    moon[i-1:i] not in ("", "0") or num(r, w+"rchi2"+suf) > 5):
                continue
            m = num(r, w+"mpro"+suf)
            err = num(r, w+"sigmpro"+suf)
            if 0 < err <= 0.36 and 0 < m < 30:
                src = col(r, "cntr", "source_id_mf", "source_id", "frame_id")
                add(archive, w.upper(), t, m, err, "mag Vega", src, table+" profile-fit")
                counts[i-1] += 1
    note(archive, "measured" if sum(counts) else "no_usable_points",
         f"{len(rows)} rows within 3 arcsec, retained {counts[0]} W1 and {counts[1]} W2. "+table)


def ztf():
    url = "https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves"
    rows = rows_csv(fetch("ztf_dr", url, params={"POS": f"CIRCLE {HOST_RA} {HOST_DEC} {RADIUS}", "FORMAT":"csv"}, timeout=90))
    groups = defaultdict(list)
    for r in rows:
        if num(r, "catflags") != 0 or not (0 < num(r, "magerr") <= .36):
            continue
        groups[(col(r, "filtercode", "filter"), col(r, "oid", "objectid"))].append(r)
    accepted = []
    for (band, oid), rr in groups.items():
        dist = median([off(num(r, "ra"), num(r, "dec")) for r in rr])
        if not dist <= 2:
            continue
        for r in rr:
            accepted.append((band, epoch(col(r,"mjd","hjd")), num(r, "magerr"), oid, r))
    # Overlapping fields / duplicate object IDs can share exposure times.
    seen = set()
    for band, t, err, oid, r in sorted(accepted, key=lambda x: x[2]):
        if not t or (band, round(t, 4)) in seen:
            continue
        seen.add((band, round(t, 4)))
        add("ZTF", band, t, num(r, "mag"), err, "mag AB", oid, "catflags=0; median centroid within 2 arcsec")
    note("ZTF", "measured" if seen else "no_usable_points", f"{len(rows)} public DR rows, {len(seen)} unique quality-filtered epochs in band; coordinate query.")


def ps1():
    base = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/"
    mean = rows_csv(fetch("ps1_mean", base+"mean.csv", params={"ra":HOST_RA,"dec":HOST_DEC,"radius":RADIUS,
                 "columns":"[objID,raMean,decMean,nDetections]","pagesize":100}))
    hits = sorted((off(num(r,"ramean"),num(r,"decmean")),col(r,"objid")) for r in mean)
    if not hits or hits[0][0] > 2:
        return note("Pan-STARRS1 DR2", "no_match", f"{len(mean)} mean objects within 3 arcsec; no nucleus within 2 arcsec")
    d, oid = hits[0]
    rr = rows_csv(fetch("ps1_detections",base+"detection.csv",params={"objID":oid,"pagesize":50000},timeout=70))
    ct = 0
    for r in rr:
        f, e, q = num(r,"psfflux"), num(r,"psffluxerr"), num(r,"psfqfperfect")
        fl = int(num(r,"filterid")) if col(r,"filterid").isdigit() else 0
        if not (f > 0 and e > 0 and f/e >= 3 and q >= .85 and 1 <= fl <= 5):
            continue
        t = epoch(col(r,"obstime"))
        if t:
            add("Pan-STARRS1 DR2","grizy"[fl-1],t, -2.5*math.log10(f/3631), 1.085736*e/f,
                "mag AB",col(r,"detectid"),"PSF flux in Jy; psfQfPerfect>=0.85; S/N>=3")
            ct += 1
    note("Pan-STARRS1 DR2", "measured" if ct else "no_usable_points", f"Nearest mean object {oid} at {d:.2f} arcsec; {ct}/{len(rr)} PSF measurements retained.")


def crts():
    params = {"RA":HOST_RA, "Dec":HOST_DEC, "Rad":0.05, "DB":"photcat", "OUT":"csv"}
    failures = []
    for domain in ("nesssi.cacr.caltech.edu","nunuku.caltech.edu"):
        try:
            raw = fetch("crts_"+domain.split(".")[0], "https://"+domain+"/cgi-bin/getcssconedb_release_img.cgi",params=params,timeout=20)
            ss = BeautifulSoup(raw, "html.parser").get_text("\n") if b"<html" in raw.lower() else raw.decode(errors="replace")
            lines = [x.strip() for x in ss.splitlines() if "," in x]
            h = next((i for i,x in enumerate(lines) if "mjd" in x.lower() and "mag" in x.lower()), None)
            if h is None:
                raise ValueError("no recognized MJD/magnitude CSV header")
            rr = rows_csv(("\n".join(lines[h:])).encode())
            n = 0
            for r in rr:
                r = {k.strip():v for k,v in r.items()}
                m,e = num(r,"mag","v"),num(r,"magerr","mag_err","mag error")
                t = epoch(col(r,"mjd"))
                if t and 0 < m < 30 and 0 < e < .36 and off(num(r,"ra"),num(r,"dec")) <= 3:
                    add("CRTS","CV",t,m,e,"mag CRTS",col(r,"css_id","id","masterid"),"Unfiltered, calibrated to V; centroid within 3 arcsec")
                    n += 1
            return note("CRTS","measured" if n else "no_usable_points",f"{len(rr)} Catalina release rows; {n} retained.")
        except Exception as e:
            failures.append(str(e)[:140])
    raise RuntimeError(" / ".join(failures))


def gaia():
    server = "https://gea.esac.esa.int/tap-server/tap/sync"
    q = ("SELECT source_id,ra,dec,has_epoch_photometry FROM gaiadr3.gaia_source "
         f"WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{HOST_RA},{HOST_DEC},{RADIUS}))")
    hits = sorted((off(num(r,"ra"),num(r,"dec")), col(r,"source_id"), col(r,"has_epoch_photometry"))
                  for r in tap("gaia_source", server, q))
    if not hits or hits[0][0] > 1.5:
        return note("Gaia DR3", "no_match", f"{len(hits)} sources in 3 arcsec; no nucleus match within 1.5 arcsec.")
    d, sid, has = hits[0]
    if has.lower() in ("false","f","0"):
        return note("Gaia DR3","no_epoch_product",f"{sid} at {d:.2f} arcsec, has_epoch_photometry={has}")
    p = {"ID":sid,"RELEASE":"Gaia DR3","RETRIEVAL_TYPE":"EPOCH_PHOTOMETRY","DATA_STRUCTURE":"COMBINED","FORMAT":"csv"}
    raw = fetch("gaia_epoch","https://gea.esac.esa.int/data-server/data",params=p,timeout=65)
    if raw.startswith(b"PK"):
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            bodies = [z.read(n) for n in z.namelist() if n.lower().endswith((".csv",".csv.gz"))]
    else:
        bodies = [raw]
    count = 0
    for b in bodies:
        if b.startswith(b"\x1f\x8b"):
            b = gzip.decompress(b)
        for r in rows_csv(b):
            if any(col(r,k).lower() in ("true","1","t") for k in ("rejected_by_photometry","rejected_by_variability")):
                continue
            t, m, f, e = (num(r,k) for k in ("time","mag","flux","flux_error"))
            if 0 < t < 15000 and 0 < m < 30 and e > 0 and f/e > 3:
                add("Gaia DR3", col(r,"band"), 55197+t, m, 1.085736*e/f, "mag Vega",sid,"Barycentric MJD(TCB); archive passband")
                count += 1
    note("Gaia DR3","measured" if count else "no_usable_points",f"{sid} at {d:.2f} arcsec; {count} retained epochs.")


def viz(key, cat, radius):
    p = {"-source":cat, "-c":CENTRE, "-c.rs":radius,"-out.max":1000,"-out.all":"1","-out.add":"_r","-oc.form":"d"}
    for i, url in enumerate((VIZ, VIZ.replace("vizier.cds.unistra.fr","vizier.cfa.harvard.edu"))):
        try:
            raw = fetch(f"{key}_{i}",url,params=p,timeout=35)
            tree = votable.parse(io.BytesIO(raw))
            if any(str(x.value).upper() == "ERROR" or str(x.name).lower() == "error"
                   for container in [tree, *tree.resources] for x in container.infos):
                raise ValueError("VizieR query error")
            tables = []
            for tb in tree.iter_tables():
                fields = [f.name.lower() for f in tb.fields]
                rr = [{k: (v.decode() if isinstance(v,bytes) else str(v)) for k,v in zip(fields,row)}
                      for row in tb.array]
                tables.append((str(tb.name or tb.ID),rr))
            return tables
        except Exception:
            if i:
                raise


def mast(key, service, params):
    payload = {"service":service,"params":params,"format":"json","pagesize":50000,"page":1}
    r = json.loads(fetch(key,MAST,data={"request":json.dumps(payload)},timeout=65))
    if r.get("status") not in ("COMPLETE", None) or "data" not in r:
        raise ValueError(str(r.get("msg", r.get("status")))[:180])
    return [{k.lower():v for k,v in x.items()} for x in r["data"]]


def galex():
    rr = mast("galex_mcat","Mast.Catalogs.Galex.Cone",{"ra":HOST_RA,"dec":HOST_DEC,"radius":RADIUS})
    count = 0
    for r in rr:
        if off(num(r,"ra"),num(r,"dec")) > 3:
            continue
        for b in ("fuv","nuv"):
            m = num(r,b+"_mag",b+"mag")
            e = num(r,b+"_magerr",b+"_mag_err",b+"magerr")
            if not (0 < m < 35 and 0 < e < .36) or num(r,b+"_artifact",b+"artifact") > 0:
                continue
            t = epoch(col(r,b+"_mjd","mjd_obs","t_min","date_obs"))
            add("GALEX MCAT", b.upper(), t, m,e,"mag AB",col(r,"objid","id"),"Catalog exposure/coadd; date only if supplied")
            count += 1
    if rr:
        AUDIT.append(("GALEX MCAT","response fields",",".join(rr[0])))
    note("GALEX MCAT", "measured" if count else "no_usable_points", f"{len(rr)} catalog records, {count} valid FUV/NUV magnitudes.")

def galex_visits():
    rr = mast("galex_observations", "Mast.Caom.Filtered.Position",
              {"columns":"obsid,obs_id,t_min,t_max,filters,target_name,obs_collection",
               "position":f"{HOST_RA}, {HOST_DEC}, 0.0001","filters":[{"paramName":"obs_collection","values":["GALEX"]}]})
    n = 0
    for r in rr:
        event("GALEX footprint",col(r,"filters"),epoch(col(r,"t_min")),kind="coverage",source=col(r,"obs_id"),
              note="MAST footprint intersection; catalog/coadd observation",end_mjd=epoch(col(r,"t_max")))
    if rr:
        ids = ",".join(col(r,"obsid") for r in rr)
        products = mast("galex_products","Mast.Caom.Products",{"obsid":ids})
        chosen = {}
        for p in products:
            fn = col(p,"productfilename").lower()
            if fn.endswith(("-xd-mcat.fits.gz","-xd-mcat.fits")):
                chosen[col(p,"datauri")] = p
        bounds = {col(r,"obsid"):r for r in rr}
        for idx,(uri,p) in enumerate(chosen.items()):
            try:
                b = fetch(f"galex_mcat_fits_{idx}","https://mast.stsci.edu/api/v0.1/Download/file",params={"uri":uri},timeout=50)
                with fits.open(io.BytesIO(gzip.decompress(b) if b.startswith(b"\x1f\x8b") else b),memmap=False) as hd:
                    tab = hd[1].data
                    keys = {k.lower():k for k in tab.names}
                    sep = 3600*np.hypot((tab[keys["ra"]]-HOST_RA)*math.cos(math.radians(HOST_DEC)),tab[keys["dec"]]-HOST_DEC)
                    if not len(sep) or np.nanmin(sep)>3:
                        continue
                    row = tab[int(np.nanargmin(sep))]
                    r = {k:str(row[v]) for k,v in keys.items()}
                    par = bounds.get(col(p,"parent_obsid"), bounds.get(col(p,"obsid"), {}))
                    for band in ("fuv","nuv"):
                        m,e = num(r,band+"_mag"), num(r,band+"_magerr")
                        if 0 < m < 35 and 0 < e < .36 and num(r,band+"_artifact",band+"_artifact_flag") <= 0:
                            add("GALEX exposure catalog",band.upper(),epoch(col(par,"t_min")),m,e,"mag AB",uri,"Tile catalog; date interval may span visits",end=epoch(col(par,"t_max")))
                            n += 1
            except Exception as exc:
                AUDIT.append(("GALEX",uri,str(exc)[:200]))
        note("GALEX footprint","measured" if n else "coverage_only",f"{len(rr)} footprint records, {len(chosen)} distinct MCAT products, {n} matched band measurements.")
    else:
        note("GALEX footprint","no_match","No MAST CAOM GALEX footprints returned.")


def uvot():
    # Search mission serendipitous detections, whose catalog aperture and flags are calibrated.
    tb = viz("uvotssc","II/339",3)
    ct = 0
    rows = 0
    for name, rr in tb:
        rows += len(rr)
        if rr:
            AUDIT.append((name, "response fields", ",".join(rr[0])))
        for r in rr:
            t = epoch(col(r,"mjd","mjd_start","date","date_obs"))
            for b in ("uvw2", "uvm2", "uvw1", "u", "b", "v"):
                # Catalogue field names vary between HEASARC and VizieR.
                m = num(r,b+"mag",b+"_ab",b+"magab",b+"_ab_mag",b)
                e = num(r,"e_"+b+"mag","e_"+b,"e_"+b+"_ab", b+"_ab_mag_err")
                flags = col(r,"q_"+b+"mag","q_"+b,b+"_quality_flag")
                if 0 < m < 32 and 0 < e < .36 and (not flags or not flags.strip("0Ff ")):
                    add("Swift UVOTSSC",b.upper(),t,m,e,"mag (catalog system)",col(r,"obsid","iauname","srcnum"),name)
                    ct += 1
    note("Swift UVOTSSC","measured" if ct else ("no_usable_points" if rows else "no_match"),
         f"{rows} catalog rows within 3 arcsec, {ct} retained bands. Serendipitous catalog has limited mission dates.")


def uvot_log():
    # A cone on observation centres gives candidate fields; it is not per-filter footprint/exposure proof.
    names = tap("uvot_table_names",HEASARC,"SELECT table_name FROM TAP_SCHEMA.tables WHERE table_name LIKE '%uv%log%'")
    matches = [col(r,"table_name") for r in names if "swift" in col(r,"table_name") or "swuv" in col(r,"table_name")]
    if not matches:
        raise ValueError("HEASARC did not return a Swift UVOT observation-log table")
    name = next((k for k in matches if "uvlog" in k),matches[0])
    rows = tap("uvot_log",HEASARC, f"SELECT * FROM {name} WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{HOST_RA},{HOST_DEC},0.21))")
    for r in rows:
        t = epoch(col(r,"time","start_time","mjd","obs_date"))
        event("Swift UVOT fields",col(r,"filter","filters"),t,kind="possible_coverage",source=col(r,"obsid","obs_id","name"),
              note=f"{name}; pointing centre within 12.6 arcmin; exposure at source unverified",end_mjd=epoch(col(r,"end_time","stop_time")))
    note("Swift UVOT fields", "coverage_only" if rows else "no_match", f"{len(rows)} log fields within 12.6 arcmin; pointings alone do not give flux or limits.")


def radio():
    jobs = [("FIRST","VIII/92/first14",4),("VLASS","J/ApJS/255/30",4)]
    for archive, cat, radius in jobs:
        try:
            tables = viz(archive.lower(),cat,radius)
            ct = 0
            for name, rr in tables:
                if rr:
                    AUDIT.append((name,"response fields",",".join(rr[0])))
                for r in rr:
                    flux = num(r,"fpeak","speak","peakflux","peak_flux","fp")
                    rms = num(r,"rms","islrms","e_speak","e_fpeak","e_peakflux")
                    if flux > 0:
                        t = epoch(col(r,"mjdmean","mjd","date-obs","dateobs","date","mean"))
                        band = "1.4 GHz" if archive=="FIRST" else "3 GHz"
                        add(archive,band,t,flux,rms,"mJy/beam",col(r,"first","vlass1qlc","id","field","name","comp","component_name"),
                            f"Peak, {name}; catalog position within {radius} arcsec")
                        ct += 1
            note(archive,"measured" if ct else ("no_usable_points" if any(rr for _,rr in tables) else "no_match"),
                 f"{ct} peak fluxes; VizieR {cat} within {radius} arcsec. Empty cone gives no flux limit. VLASS component catalog is epoch 1." if archive=="VLASS"
                 else f"{ct} peak fluxes; VizieR {cat} within {radius} arcsec. Empty cone gives no flux limit.")
        except Exception as e:
            note(archive,"unavailable",str(e)[:230])


def vlass_fields():
    server = "https://ws.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/argus/sync"
    query = ("SELECT TOP 500 obs_id,t_min,t_max,access_url,obs_collection FROM ivoa.ObsCore "
             f"WHERE INTERSECTS(s_region,POINT('ICRS',{HOST_RA},{HOST_DEC}))=1 "
             "AND obs_collection='VLASS'")
    rr = tap("vlass_cadc_coverage",server,query)
    seen = set()
    for r in rr:
        key = (col(r,"obs_id"),col(r,"t_min"),col(r,"t_max"))
        if key in seen:
            continue
        seen.add(key)
        event("VLASS CADC fields","3 GHz",epoch(key[1]),kind="coverage",source=key[0],
              note="Public footprint; catalog point flux reported separately where available",end_mjd=epoch(key[2]))
    note("VLASS CADC fields","coverage_only" if seen else "no_match",
         f"{len(seen)} distinct VLASS observation/time records; footprint intersection with CADC ObsCore.")


def summary():
    series = defaultdict(list)
    for archive,b,t,m,e,unit in POINTS:
        if unit.startswith("mag"):
            series[(archive,b)].append((t,m,e))
    lines = ["# Optical, infrared, ultraviolet and radio measurements","",
             "Positions refer to the 2MASX nucleus. WISE uses a 3 arcsec cone. PS1 uses the closest mean object within 2 arcsec and ZTF requires a median centroid within 2 arcsec for each light-curve ID. The host contributes unresolved light. Different filters, zero points, apertures and angular resolutions preclude direct cross-survey magnitude differences.","",
             "## Same-survey temporal comparison","",
             "Medians before the 2014 slew (MJD < 56658, through 2013) and near eRASS1 (MJD 58484–59214, calendar 2019–2020). Positive differences mean fainter. Dates are MJD. Values pool individual detections, so uneven cadence, correlated calibration errors and non-detections limit interpretation. First and last refer to medians in the first and last occupied 180-day bins, including sparse bins. The seasonal table gives counts and scatter.","",
             "| Archive | Band | N | MJD range | First bin | Last bin | Before 2014 (N) | 2019–20 (N) | Difference |",
             "|---|---|---:|---|---:|---:|---:|---:|---:|"]
    for band in ("W1", "W2"):
        early = series.get(("AllWISE",band), [])
        late = series.get(("NEOWISE",band), [])
        if early and late:
            series[("WISE (both releases)",band)] = early + late
    binsout = []
    def fmt(vals):
        return f"{median(vals):.3f} ({len(vals)})" if vals else "—"
    for (a,b),rr in sorted(series.items()):
        rr.sort()
        bins = defaultdict(list)
        for t,m,e in rr:
            bins[math.floor((t-55000)/180)].append((t,m,e))
        for k,entries in sorted(bins.items()):
            tt = [t for t,m,e in entries]
            mm = [m for t,m,e in entries]
            med = median(mm)
            binsout.append([a,b,min(tt),max(tt),len(mm),med,median([abs(x-med) for x in mm]),median([float(e) for t,m,e in entries if e] or [float("nan")])])
        pre = [m for t,m,e in rr if t<56658]
        near = [m for t,m,e in rr if 58484<=t<59215]
        first = [m for t,m,e in bins[min(bins)]]
        last = [m for t,m,e in bins[max(bins)]]
        diff = f"{median(near)-median(pre):+.3f}" if pre and near else "—"
        lines.append(f"| {a} | {b} | {len(rr)} | {rr[0][0]:.2f}–{rr[-1][0]:.2f} | {fmt(first)} | {fmt(last)} | {fmt(pre)} | {fmt(near)} | {diff} |")
    with (OUT/"photometry-seasons.csv").open("w",newline="") as f:
        w = csv.writer(f)
        w.writerow(["archive","band","mjd_first","mjd_last","n","median_mag","mad_mag","median_formal_error_mag"])
        w.writerows(binsout)
    lines.extend(["","## Archive results",""])
    for a,s,d in sorted(DETAIL):
        lines.append(f"- **{a}** ({s}). {d}")
    lines.extend(["","The WISE release boundary changes processing. W1/W2 magnitudes use Vega profile fits. The 2010–2011 AllWISE multiepoch catalog and 2013–2024 NEOWISE single-exposure catalog are reported separately. Residual blending, scan angle, host light and calibration changes can mimic small offsets. A shared slow change across bands warrants inspection of seasonal sampling. Absence of a catalog row alone gives no upper limit.",""])
    (OUT/"photometry.md").write_text("\n".join(lines))
    with (OUT/"photometry-archive-fields.csv").open("w",newline="") as f:
        w = csv.writer(f)
        w.writerow(["archive","item","detail"])
        w.writerows(sorted(AUDIT))


def run():
    jobs = {"AllWISE":lambda: wise("allwise_p3as_mep","AllWISE",True),
            "NEOWISE":lambda: wise("neowiser_p1bs_psd","NEOWISE",False),"ZTF":ztf,
            "Pan-STARRS1 DR2":ps1,"CRTS":crts,"Gaia DR3":gaia,"GALEX MCAT":galex,
            "GALEX footprint":galex_visits,"Swift UVOTSSC":uvot,"Swift UVOT fields":uvot_log,"radio":radio,"VLASS CADC fields":vlass_fields}
    with ThreadPoolExecutor(max_workers=8) as pool:
        futures = {pool.submit(fn):a for a,fn in jobs.items()}
        for future in as_completed(futures):
            a = futures[future]
            try:
                future.result()
            except Exception as exc:
                note(a,"unavailable",str(exc)[:250])
    summary()
