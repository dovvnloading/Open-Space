"""Host and novelty checks for the slew search. Run inside the science job."""
import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
import re
import threading
import traceback
import subprocess
import tarfile
import xml.etree.ElementTree as ET
from html.parser import HTMLParser
from urllib.parse import urljoin, urlparse
import zipfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

XMATCH = ("https://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync",
          "https://xmatch.cds.unistra.fr/xmatch/api/v1/sync")
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
NED = "https://ned.ipac.caltech.edu/tap/sync"
DATALAB = "https://datalab.noirlab.edu/query/query"
VIZIER_TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
SEARCH_RADIUS = 35.0
# Search terms apply only to papers associated with positional SIMBAD matches.
RELEVANT = re.compile(r"tidal.{0,12}disrupt|transient|changing.{1,4}(look|state)|"
                      r"(x.ray|nuclear).{0,50}(variab|flar|fading|outburst)|"
                      r"(variab|flar|fading|outburst).{0,50}(agn|nuclei|x.ray)|"
                      r"switch(ing|ed).{1,4}off", re.I)
STELLAR = re.compile(r"^(\*|Star|V\*|Fl\*|Er\*|RS\*|BY\*|Ro\*|TT\*|PM\*|"
                     r"YSO|CV\*|EB\*|Mira|PulsV\*|Flare\*|V\*\?|BYDra|RSCVn)$", re.I)
TRANSIENT = re.compile(r"^(TDE|TDe|SN|Nova|XB|XRB|HXB|LXB|CV\*|DN\*)\??$", re.I)
PUBLISHED_SEARCHES = {
    "TDE": "tidal disruption events",
    "changing-look": "changing-look AGN",
    "slew variability": "XMM slew variability",
    "slew transients": "XMM slew transient",
    "Strotjohann": "Strotjohann",
    "Sazonov": "Sazonov tidal disruption",
    "Saxton": "Saxton tidal disruption",
    "ROSAT flares": "ROSAT X-ray flares",
    "eROSITA transients": "eROSITA X-ray transients",
}


PRIOR_PAPERS = {
    "10.1093/mnras/stac714": ("Li2022", "J/MNRAS"),
    "10.1051/0004-6361/201628241": ("Strotjohann", "J/A+A"),
    "10.1111/j.1365-2966.2010.18024.x": ("Starling", "J/MNRAS"),
}


def chunks(values, n):
    for i in range(0, len(values), n):
        yield values[i:i + n]


def scalar(row, *names):
    keys = {str(k).strip().casefold(): v for k, v in row.items() if k is not None}
    for name in names:
        val = keys.get(name.casefold())
        if val is not None and str(val).strip() not in ("", "--", "nan"):
            return str(val).strip()
    return ""


def number(row, *names):
    try:
        value = float(scalar(row, *names))
        return value if math.isfinite(value) else None
    except (ValueError, TypeError):
        return None


def sky(points):
    return SkyCoord([p["ra"] for p in points] * u.deg, [p["dec"] for p in points] * u.deg)


def clipped_failure(exc):
    return re.sub(r"\s+", " ", str(exc))[:420]


class Audit:
    def __init__(self):
        self.events = []
        self.mu = threading.Lock()

    def record(self, payload):
        with self.mu:
            self.events.append(payload)

    def get(self, name, url, *, method="GET", **kwargs):
        params = kwargs.get("data") or kwargs.get("params") or {}
        # Include the submitted coordinates/query or their digest; the original
        # input and code reproduce multipart uploads deterministically.
        upload_digests = {n:hashlib.sha256((v[1].encode() if isinstance(v[1],str) else v[1])).hexdigest()
                          for n,v in kwargs.get("files",{}).items()}
        signature = hashlib.sha256(json.dumps([params,upload_digests], sort_keys=True).encode()).hexdigest()
        try:
            response = requests.request(method, url, timeout=(10, 65), **kwargs)
            raw = response.content
            self.record({"service": name, "utc":datetime.now(timezone.utc).isoformat(), "url": url, "request_sha256": signature,
                         "status": response.status_code, "bytes": len(raw),
                         "sha256": hashlib.sha256(raw).hexdigest()})
            response.raise_for_status()
            return raw
        except requests.RequestException as exc:
            self.record({"service": name, "url": url, "error": clipped_failure(exc)})
            raise


def csv_rows(raw, expected=()):
    src = raw.decode("utf-8-sig", "replace").lstrip()
    reader = csv.DictReader(io.StringIO(src))
    header = reader.fieldnames or []
    if src.startswith("<") or len(header) < 2 or any(k.lower() not in
            {h.lower() for h in header} for k in expected):
        raise ValueError("Unrecognized response: " + src[:220])
    return list(reader)


def tap(audit, name, url, query):
    return csv_rows(audit.get(name, url, method="POST", data={
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": query,
        "MAXREC": "100001"}))


def select_points(rows, controls):
    points = []
    for i, row in enumerate(rows):
        ra, dec = float(row["ra"]), float(row["dec"])
        sig = number(row, "poserr_arcsec")
        gate = min(SEARCH_RADIUS, max(5, 2.45 * (sig if sig and sig > 0 else 8.0)))
        points.append(dict(k=i * 5, i=i, shift=0, ra=ra, dec=dec, gate=gate,
                           ero_name=scalar(row,"ero_name"), ero_ra=number(row,"ero_ra"),ero_dec=number(row,"ero_dec")))
        if controls:
            c = SkyCoord(ra * u.deg, dec * u.deg)
            for j, pa in enumerate((0, 90, 180, 270), 1):
                s = c.directional_offset_by(pa * u.deg, 240 * u.arcsec)
                points.append(dict(k=i * 5 + j, i=i, shift=j, ra=s.ra.deg,
                                   dec=s.dec.deg, gate=gate))
    return points


def match_fields(input_points, rows, names, *, radius=SEARCH_RADIUS):
    """Return all in-cone rows with a great-circle separation in arcseconds."""
    results = defaultdict(list)
    fields, coords = [], []
    for row in rows:
        ra, dec = number(row, *names[0]), number(row, *names[1])
        if ra is None or dec is None or abs(dec) > 90:
            continue
        coords.append({"ra": ra, "dec": dec})
        fields.append(row)
    if rows and not coords:
        raise ValueError("Catalogue rows have no recognized decimal coordinates")
    if coords and input_points:
        ip, ic, d, _ = search_around_sky(sky(input_points), sky(coords), radius * u.arcsec)
        for i, j, sep in zip(ip, ic, d.arcsec):
            results[input_points[i]["k"]].append(dict(fields[j], sep_arcsec=round(sep, 3)))
    return results


def xmatch(audit, cat, points, fields):
    result, status = defaultdict(list), {}
    failures,reason = 0,""
    for part in chunks(points, 900):
        if failures >= 2:
            status.update({p["k"]:"service unavailable after repeated failures: "+reason for p in part})
            continue
        buf = io.StringIO()
        writer = csv.writer(buf)
        writer.writerow(("qkey", "input_ra", "input_dec"))
        for p in part:
            writer.writerow((p["k"], p["ra"], p["dec"]))
        try:
            data = {"request": "xmatch", "distMaxArcsec": str(SEARCH_RADIUS), "RESPONSEFORMAT": "csv",
                    "cat2": "vizier:" + cat, "colRA1": "input_ra", "colDec1": "input_dec"}
            data["RESPONSEFORMAT"] = "csv"
            for mirror in XMATCH:
                try:
                    raw = audit.get(cat, mirror, method="POST", data=data,
                                    files={"cat1": ("positions.csv", buf.getvalue(), "text/csv")})
                    rows = csv_rows(raw, ("qkey","angDist"))
                    break
                except Exception:
                    if mirror == XMATCH[-1]:
                        raise
            failures = 0
            for row in rows:
                k = int(row["qkey"])
                item = {field: scalar(row, *aliases) for field, aliases in fields.items()}
                item["sep_arcsec"] = float(row["angDist"])
                result[k].append(item)
            status.update({p["k"]: "ok" for p in part})
        except Exception as exc:
            reason = clipped_failure(exc)
            failures += 1
            status.update({p["k"]: reason for p in part})
    return result, status


def region_query(audit, name, points, make_query, names):
    result, status = defaultdict(list), {}
    failures,reason = 0,""
    for part in chunks(points, 35):
        if failures >= 2:
            status.update({p["k"]:"service unavailable after repeated failures: "+reason for p in part})
            continue
        try:
            rows = make_query(audit, name, part)
            if len(rows) >= 100000:
                raise ValueError("service row limit reached")
            matched = match_fields(part, rows, names)
            failures = 0
            for k, found in matched.items():
                result[k].extend(found)
            status.update({p["k"]: "ok" for p in part})
        except Exception as exc:
            reason = clipped_failure(exc)
            failures += 1
            status.update({p["k"]: reason for p in part})
    return result, status


def cones(part):
    return " OR ".join(f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{p['ra']:.9f},"
                       f"{p['dec']:.9f},{SEARCH_RADIUS / 3600:.9f}))" for p in part)


def simbad_query(audit, name, part):
    return tap(audit, name, SIMBAD, "SELECT oid,main_id,otype,ra,dec FROM basic WHERE " + cones(part))


def ned_query(audit, name, part):
    rows = tap(audit, name, NED, "SELECT * FROM NEDTAP.objdir WHERE " + cones(part))
    return [{"ra": scalar(r, "ra"), "dec": scalar(r, "dec"),
             "preferred_name": scalar(r, "prefname", "objname", "name"),
             "type": scalar(r, "prefphytype", "pretype", "objtype", "type"), "z": scalar(r, "z", "redshift")}
            for r in rows]


def legacy_query(audit, name, part):
    # DR10 includes earlier northern imaging. Deduplication uses stable IDs.
    where = " OR ".join(f"q3c_radial_query(ra,dec,{p['ra']:.9f},{p['dec']:.9f},"
                         f"{SEARCH_RADIUS / 3600:.9f})" for p in part)
    query = "SELECT ls_id,type,ra,dec,flux_r,flux_z FROM ls_dr10.tractor WHERE " + where
    raw = audit.get(name, DATALAB, method="GET", params={"sql": query, "ofmt": "csv",
                    "async": "False"}, headers={"X-DL-AuthToken": "anonymous.0.0.anon_access"})
    return csv_rows(raw, ("ls_id", "ra", "dec"))


GAIA = {
    "id": ["Source", "source_id"], "ra": ["RAdeg", "RA_ICRS", "RAJ2000", "_RAJ2000"],
    "dec": ["DEdeg", "DE_ICRS", "DEJ2000", "_DEJ2000"], "plx": ["Plx"], "e_plx": ["e_Plx"],
    "pmra": ["pmRA"], "e_pmra": ["e_pmRA"], "pmdec": ["pmDE"], "e_pmdec": ["e_pmDE"],
    "ruwe": ["RUWE"], "gmag": ["Gmag"], "variable": ["VarFlag", "VARIABLE"]}

WISE = {"id": ["AllWISE", "WISE", "designation"], "ra": ["RAJ2000", "RAdeg"],
        "dec": ["DEJ2000", "DEdeg"], "w1": ["W1mag"], "w2": ["W2mag"], "ccf": ["ccf"]}

VSX = {"id": ["Name", "OID"], "type": ["Type"], "ra": ["RAJ2000"], "dec": ["DEJ2000"]}


# Independent public catalogue front ends. Each route must recover a named
# object before empty cone results on that route are accepted.
SIMBAD_HOSTS = ("simbad.cfa.harvard.edu","simbad.u-strasbg.fr","simbad.cds.unistra.fr")
SIMBAD_TAPS = ["https://"+h+"/simbad/sim-tap/sync" for h in SIMBAD_HOSTS]
PROBES = {"Gaia":(269.4485,4.73942,"4472832130942575872"),
          "WISE":(187.277915,2.052388,"122906"),
          "VSX":(291.3663,42.78436,"rrlyr"),
          "SIMBAD":(187.277915,2.052388,"3c273")}
GAIA_NATIVE = {"id":["source_id"],"ra":["ra"],"dec":["dec"],"plx":["parallax"],
               "e_plx":["parallax_error"],"pmra":["pmra"],"e_pmra":["pmra_error"],
               "pmdec":["pmdec"],"e_pmdec":["pmdec_error"],"ruwe":["ruwe"],
               "gmag":["phot_g_mean_mag"],"variable":["phot_variable_flag"]}
WISE_NATIVE = {"id":["designation","source_id"],"ra":["ra"],"dec":["dec"],
               "w1":["w1mpro"],"w2":["w2mpro"],"ccf":["cc_flags"]}


def native_route(schema,fields,service,url,columns=None):
    selected = columns or ",".join(v[0] for v in fields.values())
    def run(audit,name,part):
        if service=="datalab":
            where=" OR ".join(f"q3c_radial_query(ra,dec,{p['ra']:.8f},{p['dec']:.8f},{SEARCH_RADIUS/3600:.8f})"
                             for p in part)
            query=f"SELECT {selected} FROM {schema} WHERE "+where
            raw=audit.get(name+" "+schema,url,params={"sql":query,"ofmt":"csv","async":"False"},
                          headers={"X-DL-AuthToken":"anonymous.0.0.anon_access"})
            rows=csv_rows(raw)
        else:
            rows=tap(audit,name+" "+schema,url,f"SELECT {selected} FROM {schema} WHERE "+cones(part))
        cleaned=[{k:scalar(row,*names) for k,names in fields.items()} for row in rows]
        if cleaned and not all(r.get("id") and r.get("ra") and r.get("dec") for r in cleaned):
            raise ValueError("unrecognized native source or coordinate fields")
        return cleaned
    return run


def simbad_tap_route(url):
    def run(audit,name,part):
        rs=tap(audit,name,url,"SELECT oid,main_id,otype,ra,dec FROM basic WHERE "+cones(part))
        return [dict(r,bib_tap=url) for r in rs]
    return run


def votable_rows(raw):
    from astropy.io.votable import parse
    xml = ET.fromstring(raw)
    if any(str(info.get("value","")).upper() in ("ERROR","OVERFLOW")
           for info in xml.iter() if info.tag.split("}")[-1]=="INFO"):
        raise ValueError("error or truncated VO response")
    if xml.tag.split("}")[-1].upper()!="VOTABLE":
        raise ValueError("non-VOTable response")
    tables=list(parse(io.BytesIO(raw),verify="ignore").iter_tables())
    if not tables:
        if re.search(r"(no (?:astronomical )?object|no matching|no result|no record)",str(raw[:16000]),re.I):
            return []
        raise ValueError("VOTable has no catalogue table")
    out=[]
    for extension in tables:
        tab=extension.to_table(use_names_over_ids=True)
        for row in tab:
            out.append({k:str(row[k]).strip() for k in tab.colnames if not np.ma.is_masked(row[k])})
    return out


def harvard_vsx_tsv(audit,name,part):
    rows=[]
    for p in part:
        raw=audit.get(name,"https://vizier.cfa.harvard.edu/viz-bin/asu-tsv",params={
            "-source":"B/vsx/vsx","-c":f"{p['ra']:.8f} {p['dec']:+.8f}","-c.rs":SEARCH_RADIUS,
            "-out":"Name,RAJ2000,DEJ2000,Type","-out.max":1000})
        txt=raw.decode("utf-8-sig","replace")
        scoped="#Name: B/vsx" in txt and "#RESOURCE=yCat_" in txt
        if not scoped or re.search(r"^#(?:ERROR|INFO.*(?:QUERY_STATUS=ERROR|truncated|overflow))",txt,re.I|re.M):
            raise ValueError("unrecognized or truncated VizieR VSX response")
        lines=[line for line in txt.splitlines() if line and not line.startswith("#")]
        if not lines:
            continue
        header=next(csv.reader([lines[0]],delimiter="\t"))
        if not {"Name","RAJ2000","DEJ2000","Type"} <= set(header):
            raise ValueError("VizieR VSX changed column names")
        for r in csv.DictReader(lines[1:],fieldnames=header,delimiter="\t"):
            ra,dec=number(r,"RAJ2000"),number(r,"DEJ2000")
            if ra is not None and dec is not None:
                rows.append(dict(id=scalar(r,"Name"),ra=ra,dec=dec,type=scalar(r,"Type")))
    return rows


def vizier_vsx_route(server):
    def run(audit,name,part):
        rs=[]
        for p in part:
            raw=audit.get(name,server+"/viz-bin/votable",params={"-source":"B/vsx/vsx",
                "-c":f"{p['ra']:.8f} {p['dec']:+.8f}","-c.rs":SEARCH_RADIUS,
                "-out":"Name,Type,_RAJ2000,_DEJ2000,RAJ2000,DEJ2000","-out.max":1000})
            for row in votable_rows(raw):
                ra,dec=parse_coord(scalar(row,"_RAJ2000","RAJ2000"),scalar(row,"_DEJ2000","DEJ2000"))
                rs.append(dict(id=scalar(row,"Name","OID"),type=scalar(row,"Type"),ra=ra,dec=dec))
        return rs
    return run


def aavso_vsx_route(audit,name,part):
    rs=[]
    for p in part:
        raw=audit.get(name,"https://www.aavso.org/vsx/index.php",params={"view":"api.list","format":"json",
                      "ra":p["ra"],"dec":p["dec"],"radius":SEARCH_RADIUS/3600})
        payload=json.loads(raw)
        if "VSXObjects" not in payload:
            raise ValueError("unknown AAVSO response schema")
        content=payload["VSXObjects"] or {}
        records=content.get("VSXObject",[]) if isinstance(content,dict) else content
        if isinstance(records,dict):
            records=[records]
        for r in records or []:
            ra,dec=parse_coord(scalar(r,"RA2000","RA","RAJ2000"),scalar(r,"Declination2000","Declination","Dec","DEJ2000"))
            rs.append(dict(id=scalar(r,"Name"),type=scalar(r,"VariabilityType","Type"),ra=ra,dec=dec))
    return rs


def sim_coo_route(host):
    def run(audit,name,part):
        rs=[]
        for p in part:
            raw=audit.get(name,"https://"+host+"/simbad/sim-coo",params={
                    "Coord":f"{p['ra']:.8f} {p['dec']:+.8f}","CooFrame":"ICRS","CooEpoch":2000,
                    "Radius":SEARCH_RADIUS,"Radius.unit":"arcsec","output.format":"VOTable",
                    "output.params":"main_id,otype,ra(d),dec(d)"})
            for r in votable_rows(raw):
                ra,dec=parse_coord(scalar(r,"RA_d","RA(d)","RA","ra_d_ICRS_J2000.0"),
                                   scalar(r,"DEC_d","DEC(d)","DEC","dec_d_ICRS_J2000.0"))
                rs.append(dict(main_id=scalar(r,"MAIN_ID","identifier"),otype=scalar(r,"OTYPE","OTYPE_S"),
                               oid=scalar(r,"oid","oidref"),ra=ra,dec=dec))
        return rs
    return run


def public_routes(name):
    if name=="Gaia":
        return [("Data Lab Gaia DR3",native_route("gaia_dr3.gaia_source",GAIA_NATIVE,"datalab",DATALAB)),
                ("ESA Gaia DR3",native_route("gaiadr3.gaia_source",GAIA_NATIVE,"tap","https://gea.esac.esa.int/tap-server/tap/sync")),
                ("AIP Gaia DR3",native_route("gaiadr3.gaia_source",GAIA_NATIVE,"tap","https://gaia.aip.de/tap/sync"))]
    if name=="WISE":
        return [("Data Lab AllWISE",native_route("allwise.source",WISE_NATIVE,"datalab",DATALAB)),
                ("IRSA AllWISE",native_route("allwise_p3as_psd",WISE_NATIVE,"tap","https://irsa.ipac.caltech.edu/TAP/sync"))]
    if name=="VSX":
        return [("Harvard VizieR ASU TSV",harvard_vsx_tsv)]+[(h,vizier_vsx_route("https://"+h)) for h in
                 ("vizier.cfa.harvard.edu","vizier.iucaa.in","vizier.cds.unistra.fr")]+[("AAVSO direct",aavso_vsx_route)]
    return [(SIMBAD_TAPS[0],simbad_tap_route(SIMBAD_TAPS[0])),
            (SIMBAD_HOSTS[0],sim_coo_route(SIMBAD_HOSTS[0]))]+[(url,simbad_tap_route(url)) for url in SIMBAD_TAPS[1:]]+[
                (h,sim_coo_route(h)) for h in SIMBAD_HOSTS[1:]]


def verified_service(audit,name,points):
    ra,dec,expected=PROBES[name]
    center=[dict(ra=ra,dec=dec,k=-1)]
    result,status={},{}
    reasons=[]
    for tag,query in public_routes(name):
        missing=[p for p in points if status.get(p["k"])!="ok"]
        if not missing:
            break
        try:
            sample=query(audit,name,center)
            def is_control(r):
                # SIMBAD often returns the quasar's alternate display ID
                # QSO B1226+023, rather than its 3C identifier.
                label=re.sub(r"[^a-z0-9]","",scalar(r,"id","main_id").lower())
                good=expected in label or (name=="SIMBAD" and ("1226023" in label or "122906" in label))
                return good and distance(r,center[0])<5
            if not any(is_control(r) for r in sample):
                raise ValueError("named positive control absent; fields="+str(sample[:1])[:400])
            audit.record({"service":name,"route":tag,"control":expected,"status":"passed"})
            if name=="VSX" or "sim-coo" in query.__qualname__ or (name=="SIMBAD" and tag in SIMBAD_HOSTS):
                # Single-cone services can fail on one coordinate. Keep other
                # successful cones and avoid dropping an entire 35-source block.
                def one(p):
                    try:
                        rr=query(audit,name,[p])
                        m=match_fields([p],rr,(("ra",),("dec",)))
                        return m.get(p["k"],[]),"ok"
                    except Exception as exc:
                        return [],clipped_failure(exc)
                with ThreadPoolExecutor(max_workers=3) as pool:
                    for p,(hits,state) in zip(missing,pool.map(one,missing)):
                        result[p["k"]],status[p["k"]]=hits,state
            else:
                hits,checks=region_query(audit,name,missing,query,(("ra",),("dec",)))
                for k,state in checks.items():
                    if state=="ok":
                        result[k]=hits.get(k,[])
                    status[k]=state
        except Exception as exc:
            reasons.append(tag+": "+clipped_failure(exc))
            audit.record({"service":name,"route":tag,"status":"failed","error":clipped_failure(exc)})
    if not status or any(status.get(p["k"])!="ok" for p in points):
        missing=[p for p in points if status.get(p["k"])!="ok"]
        if name in ("Gaia","WISE","VSX"):
            cat,fields={"Gaia":("I/355/gaiadr3",GAIA),"WISE":("II/328/allwise",WISE),"VSX":("B/vsx/vsx",VSX)}[name]
            probe,ps=xmatch(audit,cat,center,fields)
            if ps.get(-1)=="ok" and any(expected in re.sub(r"[^a-z0-9]","",q["id"].lower()) and
                                         q["sep_arcsec"]<5 for q in probe.get(-1,[])):
                hits,st=xmatch(audit,cat,missing,fields)
                result.update(hits)
                status.update(st)
            else:
                reasons.append("CDS XMatch positive control failed")
        for p in missing:
            status.setdefault(p["k"]," | ".join(reasons) or "validated service unavailable")
    return result,status


def stellar_astrometry(row):
    # Gaia proper motions are the two projected components in mas/yr. Requiring
    # one significant component avoids treating a noisy vector length as Gaussian.
    plx, ep = number(row, "plx"), number(row, "e_plx")
    ruwe = number(row, "ruwe")
    if ruwe and ruwe > 1.4:
        return False
    if plx and ep and ep > 0 and plx / ep >= 5:
        return True
    for field in ("pmra", "pmdec"):
        v, e = number(row, field), number(row, "e_" + field)
        if v and e and e > 0 and abs(v) >= 5 and abs(v) / e >= 5:
            return True
    return False


def associate_refs(audit, simbad, points):
    ids = sorted({int(r["oid"]) for p in points for r in simbad.get(p["k"], [])
                  if r.get("oid") and float(r["sep_arcsec"]) <= p["gate"]})
    refs, checks = defaultdict(list), {}
    failures = 0
    for part in chunks(ids, 60):
        if failures >= 2:
            checks.update({k:"service unavailable after repeated bibliography failures" for k in part})
            continue
        keys = ",".join(map(str, part))
        try:
            query = ("SELECT h.oidref,r.bibcode,r.title FROM has_ref AS h JOIN ref AS r "
                     f"ON h.oidbibref=r.oidbib WHERE h.oidref IN ({keys})")
            mirrors = list(dict.fromkeys([r["bib_tap"] for hits in simbad.values() for r in hits if r.get("bib_tap")] + SIMBAD_TAPS))
            for url in mirrors:
                try:
                    rows = tap(audit,"SIMBAD-bibliography",url,query)
                    break
                except Exception:
                    if url==mirrors[-1]:
                        raise
            if len(rows) >= 100000:
                raise ValueError("bibliography truncated")
            for r in rows:
                if RELEVANT.search(r.get("title") or ""):
                    refs[int(r["oidref"])].append({"bibcode": r["bibcode"], "title": r["title"]})
            failures = 0
            checks.update({k: "ok" for k in part})
        except Exception as exc:
            failures += 1
            checks.update({k: clipped_failure(exc) for k in part})
    return refs, checks


def parse_coord(ra, dec):
    try:
        rd, dd = float(ra), float(dec)
        if 0 <= rd < 360 and abs(dd) <= 90:
            return rd, dd
    except (TypeError, ValueError):
        pass
    pos = SkyCoord(str(ra).replace(" ", ":"), str(dec).replace(" ", ":"),
                   unit=(u.hourangle, u.deg))
    return pos.ra.deg, pos.dec.deg


def tns_snapshot(audit,points):
    """Search public Fink TNS snapshots for positive evidence if the official cone is down.

    A broker can preserve only a subset or a dated snapshot. Its zero-hit
    positions remain unverified against the complete TNS, even if our positive
    control is present in the snapshot.
    """
    import pandas as pd
    results=defaultdict(list)
    seen=set()
    for repo in ("astrolabsoftware/fink-science","astrolabsoftware/fink-broker"):
        try:
            api="https://api.github.com/repos/"+repo
            rev=json.loads(audit.get("TNS-mirror-revision",api+"/commits/HEAD"))["sha"]
            tree=json.loads(audit.get("TNS-mirror-inventory",api+"/git/trees/"+rev,
                                      params={"recursive":"1"})).get("tree",[])
            assets=[t for t in tree if t["type"]=="blob" and "tns" in t["path"].lower() and
                    re.search(r"\.(parquet|csv|csv\.gz|csv\.zip)$",t["path"],re.I)]
            for t in assets:
                if t["sha"] in seen:
                    continue
                seen.add(t["sha"])
                path=t["path"]
                url="https://raw.githubusercontent.com/"+repo+"/"+rev+"/"+path
                data=audit.get("TNS-mirror",url)
                if data.startswith(b"version https://git-lfs.github.com/"):
                    url="https://media.githubusercontent.com/media/"+repo+"/"+rev+"/"+path
                    data=audit.get("TNS-mirror-LFS",url)
                if path.lower().endswith(".parquet"):
                    table=pd.read_parquet(io.BytesIO(data))
                else:
                    compression="gzip" if data.startswith(b"\x1f\x8b") else ("zip" if data.startswith(b"PK") else None)
                    table=pd.read_csv(io.BytesIO(data),compression=compression,low_memory=False)
                cols={str(c).lower().strip():c for c in table.columns}
                ra=next((cols[c] for c in ("ra","ra_deg","radeg") if c in cols),None)
                dec=next((cols[c] for c in ("dec","dec_deg","decdeg","decl") if c in cols),None)
                name=next((cols[c] for c in ("name","fullname","tns_name","tnsname") if c in cols),None)
                if not ra or not dec or not name:
                    continue
                names=table[name].astype(str).str.lower().str.replace(" ","",regex=False)
                cow=table[names.str.endswith("2018cow")]
                if not len(cow):
                    continue
                a,b=parse_coord(cow.iloc[0][ra],cow.iloc[0][dec])
                if distance({"ra":a,"dec":b},{"ra":244.000925,"dec":22.268008})>5:
                    continue
                objs=[]
                typ=next((cols[c] for c in ("type","obj. type","tns_type") if c in cols),None)
                date=next((cols[c] for c in ("discovery date (ut)","discoverydate","discovery_date") if c in cols),None)
                for values in table[[name,ra,dec]+([typ] if typ else [])].itertuples(index=False,name=None):
                    try:
                        a,b=parse_coord(values[1],values[2])
                    except (ValueError,TypeError):
                        continue
                    objs.append({"name":str(values[0]),"ra":a,"dec":b,"type":str(values[3]) if typ else "",
                                 "source_snapshot":url})
                hits=match_fields(points,objs,(("ra",),("dec",)))
                for k,rows in hits.items():
                    results[k].extend(rows)
                audit.record({"service":"TNS-broker-coverage","exact_snapshot":url,"rows":len(objs),
                              "total":len(table),"positive_control":"AT 2018cow",
                              "latest_event":str(table[date].dropna().astype(str).max()) if date else "",
                              "columns":list(table.columns),"full_TNS_coverage":False})
        except Exception as exc:
            audit.record({"service":"TNS-broker-snapshot","repo":repo,"error":clipped_failure(exc)})
    return results


def backup_tns(audit,points,official_error):
    # Independent routes can each contribute a real positive. An authoritative
    # no-hit on one source does not erase a match from a historic snapshot.
    from tns_html import tns_html_lookup
    from tns_mirror import tns_lookup
    hits=defaultdict(list)
    checks={}
    def partial():
        return tns_snapshot(audit,points),{}
    routes={"official_html":lambda:tns_html_lookup(points,audit),
            "recorded_census":lambda:tns_lookup(points,audit),
            "broker_snapshot":partial}
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures={tag:pool.submit(fn) for tag,fn in routes.items()}
        for tag,future in futures.items():
            try:
                found,states=future.result()
                for k,rs in found.items():
                    hits[k].extend(dict(r,tns_lookup_route=tag) for r in rs)
                for k,state in states.items():
                    checks.setdefault(k,{})[tag]=state
            except Exception as exc:
                audit.record({"service":"TNS-secondary","route":tag,"error":clipped_failure(exc),
                              "traceback":traceback.format_exc()})
    final={}
    for p in points:
        k=p["k"]
        states=checks.get(k,{})
        if "ok" in states.values():
            final[k]="ok"
        else:
            dated=sorted(s for s in states.values() if isinstance(s,str) and s.startswith("complete_through_"))
            final[k]=dated[-1] if dated else " | ".join(str(v) for v in states.values()) or official_error
        audit.record({"service":"TNS-field-coverage","field":k,"routes":states,"final":final[k],
                      "positional_matches":len(hits.get(k,[]))})
    return hits,final


def tns(audit, points, previous):
    try:
        raw = audit.get("TNS", "https://www.wis-tns.org/system/files/tns_public_objects/tns_public_objects.csv.zip")
        archive = zipfile.ZipFile(io.BytesIO(raw))
        file = next(s for s in archive.namelist() if s.lower().endswith(".csv"))
        # First line of the official export may contain a standalone timestamp.
        stream = io.TextIOWrapper(archive.open(file), encoding="utf-8-sig")
        header = next(stream)
        if "," not in header:
            header = next(stream)
        columns = next(csv.reader([header]))
        if "ra" not in {s.lower() for s in columns}:
            raise ValueError("TNS export has unrecognized columns")
        objs = []
        missed = 0
        for r in csv.DictReader(stream, fieldnames=columns):
            try:
                ra, dec = parse_coord(scalar(r, "RA"), scalar(r, "DEC", "declination"))
                objs.append({"ra": ra, "dec": dec, "name": scalar(r, "Name Prefix") + scalar(r, "Name"),
                             "type": scalar(r, "Obj. Type", "Type"), "date": scalar(r,"Discovery Date (UT)")})
            except (TypeError, ValueError):
                missed += 1
        if not any("2018cow" in obj["name"].replace(" ","") and
                   distance(obj,{"ra":244.000925,"dec":22.268008})<5 for obj in objs):
            raise ValueError("TNS export missing positive control AT 2018cow")
        matched = match_fields(points, objs, (("ra",), ("dec",)))
        audit.record({"service":"TNS", "parsed":len(objs), "unparsed":missed})
        return matched, {p["k"]: "ok" if not missed else f"partial ({missed} unreadable)" for p in points}
    except Exception as exc:
        archive_error = clipped_failure(exc)
    # AT 2018cow validates the TNS coordinate convention, filters and returned
    # coordinate fields before any empty cone is counted as checked.
    def cone(ra, dec):
        c = SkyCoord(ra * u.deg, dec * u.deg)
        raw = audit.get("TNS-cone", "https://www.wis-tns.org/search", params={
            "ra": c.ra.to_string(unit=u.hourangle, sep=":", precision=3),
            "decl": c.dec.to_string(sep=":", alwayssign=True, precision=2),
            "radius": SEARCH_RADIUS, "coords_unit": "arcsec", "format": "csv", "num_page": 100})
        rows = csv_rows(raw)
        header = {s.strip().lower() for s in next(csv.reader([raw.decode("utf-8-sig").splitlines()[0]]))}
        if not {"ra", "dec"} <= header or "name" not in header or len(rows) >= 100:
            raise ValueError("TNS cone schema/truncation error")
        return rows

    result, status = defaultdict(list), {}
    try:
        control = cone(244.000925, 22.268008)
        if not any("2018cow" in scalar(r, "Name").replace(" ", "") for r in control):
            raise ValueError("TNS official search missed positive control AT 2018cow")
    except Exception as exc:
        return backup_tns(audit,points,archive_error+" | "+clipped_failure(exc))

    def lookup(p):
        hits = []
        try:
            rows = cone(p["ra"], p["dec"])
            for r in rows:
                ra, dec = parse_coord(scalar(r, "RA"), scalar(r, "DEC", "declination"))
                center = SkyCoord(p["ra"] * u.deg, p["dec"] * u.deg)
                sep = center.separation(SkyCoord(ra*u.deg, dec*u.deg)).arcsec
                if sep > SEARCH_RADIUS + 0.2:
                    raise ValueError("TNS search returned an object outside requested cone")
                hits.append({"name": scalar(r,"Name"), "type":scalar(r,"Type","Obj. Type"),
                             "ra":ra,"dec":dec,"sep_arcsec":round(sep,3)})
            return hits, "ok", False
        except Exception as exc:
            blocked = isinstance(exc, requests.HTTPError) and exc.response.status_code in (401,403,429)
            return hits, archive_error + " | " + clipped_failure(exc), blocked

    coords = defaultdict(list)
    for p in points:
        # A successful bulk download above covers even the disqualified inputs.
        # On the per-target fallback, avoid hammering TNS about known X-ray binaries.
        if previous.get(p["k"],False):
            status[p["k"]] = "previously classified; TNS bulk unavailable"
        else:
            coords[(p["ra"],p["dec"])].append(p)
    with ThreadPoolExecutor(max_workers=2) as pool:
        blocked = ""
        for group in chunks(list(coords.values()), 2):
            if blocked:
                for members in group:
                    for q in members:
                        status[q["k"]] = "not attempted after service block: "+blocked
                continue
            for members, (hits, state, stop) in zip(group,pool.map(lookup,[arr[0] for arr in group])):
                for q in members:
                    result[q["k"]].extend(hits)
                    status[q["k"]] = state
                if stop:
                    blocked = state
    missing=[p for p in points if status.get(p["k"])!="ok" and not previous.get(p["k"])]
    if missing:
        supplement,checks=backup_tns(audit,missing,"official cone unavailable")
        for p in missing:
            k=p["k"]
            result[k].extend(supplement.get(k,[]))
            status[k]=checks[k]
    return result, status


def exact_papers(audit):
    records = {}
    for doi,(label,journal) in PRIOR_PAPERS.items():
        try:
            msg = json.loads(audit.get("Crossref-DOI", "https://api.crossref.org/works/"+doi))["message"]
            page = str(msg.get("page") or msg.get("article-number") or "").split("-")[0].strip()
            volume = str(msg.get("volume","")).strip()
            resource = journal+"/"+volume+"/"+page if page and volume else ""
            records[doi] = (msg,resource,label)
            audit.record({"service":"DOI-to-CDS", "doi":doi, "title":msg.get("title"),
                          "volume":volume,"first_page":page,"cds_lookup":resource})
        except Exception as exc:
            audit.record({"service":"DOI-to-CDS","doi":doi,"error":clipped_failure(exc)})
    return records


def publication_catalogs(audit,papers):
    from astroquery.vizier import VizierClass
    viz = VizierClass(columns=["*", "_RAJ2000", "_DEJ2000"], row_limit=20001)
    viz.TIMEOUT = 28
    catalogues = {}
    coverage = {}
    mirror = VizierClass(columns=["*", "_RAJ2000", "_DEJ2000"],row_limit=20001,
                         vizier_server="vizier.cfa.harvard.edu")
    mirror.TIMEOUT = 28
    for label, terms in PUBLISHED_SEARCHES.items():
        try:
            try:
                found = viz.find_catalogs(terms)
                if not found:
                    raise ValueError("no primary metadata records")
            except Exception:
                found = mirror.find_catalogs(terms)
            for name, resource in (found or {}).items():
                # Some VizieR metadata replies include a nameless RESOURCE
                # for status/description. It has no catalogue to retrieve.
                if not isinstance(name,str) or "/" not in name:
                    continue
                desc = str(resource.description or "")
                if label == "Strotjohann" and not re.search(r"slew|tidal|AGN|X-ray", desc, re.I):
                    continue
                catalogues.setdefault(name, set()).add(label)
            coverage[label] = {"catalogue_count":len(found), "status":"ok"}
        except Exception as exc:
            coverage[label] = {"status":clipped_failure(exc)}
    # Independent published sample with explicit archive identifier. The returned
    # table shape and coordinates still have to pass the same parsing checks.
    # These release identifiers point to the van Velzen and Hammerstein
    # population studies. Only a successfully parsed positional table counts.
    catalogues.setdefault("J/ApJ/908/4",set()).add("TDE")
    catalogues.setdefault("J/ApJ/942/9",set()).add("TDE")
    catalogues.setdefault("J/MNRAS/508/3820",set()).add("Sazonov")
    for _,resource,label in papers.values():
        if resource:
            catalogues.setdefault(resource,set()).add(label)
    audit.record({"service":"VizieR-publication-discovery", "search_terms":PUBLISHED_SEARCHES,
                  "coverage":coverage, "catalogues":sorted(catalogues)})
    return (viz,mirror), catalogues


def positional_columns(tab):
    keys = set(tab.colnames)
    for pair in (("_RAJ2000", "_DEJ2000"), ("RAJ2000", "DEJ2000"), ("RA_ICRS", "DE_ICRS"),
                 ("RAdeg", "DEdeg"), ("RA", "DE"), ("RA", "Dec")):
        if set(pair) <= keys:
            return pair
    coords = {}
    for name in keys:
        ucd = str(tab[name].meta.get("ucd", "")).lower()
        if ucd in ("pos.eq.ra;meta.main", "pos.eq.dec;meta.main"):
            coords[ucd.split(";")[0]] = name
    return coords.get("pos.eq.ra"), coords.get("pos.eq.dec")


def hash_table(tab):
    digest = hashlib.sha256()
    digest.update(json.dumps(tab.colnames).encode() + b"\n")
    for row in tab:
        digest.update(json.dumps([str(x) for x in row], ensure_ascii=True).encode() + b"\n")
    return digest.hexdigest()


class ArticleHTML(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.parts = []

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if tag == "a" and a.get("href"):
            self.links.append(a["href"])
        if tag == "meta" and a.get("name") in ("citation_fulltext_html_url","citation_pdf_url"):
            self.links.append(a.get("content",""))

    def handle_data(self, data):
        self.parts.append(data)


IAU_NAME = re.compile(r"\b(XMMSL[12]?|[12]RXS|eRASS[0-9]?|[234]XMM|SRGe|SDSS|RX)"
                      r"\s*J?\s*(\d{6}(?:\.\d{1,3})?)\s*([+-])\s*(\d{6}(?:\.\d{1,2})?)",re.I)


def article_objects(audit, points, papers):
    """Find coordinate-coded objects in open primary papers; a hit flags review."""
    searches = [
        ("slew","The evolution of highly variable AGN in the XMM-Newton slew survey Strotjohann"),
        ("slew","XMM-Newton slew survey ROSAT X-ray transient sources tidal disruption Saxton"),
        ("slew","highly variable X-ray sources XMM-Newton Slew Survey Li"),
        ("TDE","First tidal disruption events discovered by SRG eROSITA Sazonov"),
    ]
    matches, searched, seendois = defaultdict(list),set(),set()
    allsearch = [("exact",doi) for doi in papers] + searches
    for kind, phrase in allsearch:
        try:
            if kind=="exact":
                records = [papers[phrase][0]]
            else:
                raw = audit.get("Crossref", "https://api.crossref.org/works",
                                params={"query.bibliographic":phrase,"rows":5,"select":"DOI,title,link,resource"})
                records = json.loads(raw)["message"]["items"]
            for rec in records:
                doi, title = rec["DOI"], rec.get("title",[""])[0]
                if doi in seendois:
                    continue
                if kind == "slew" and not ("slew" in title.lower() and re.search(
                        r"varia[bnt]|transient|disrupt|flare|outburst",title,re.I)):
                    continue
                if kind == "TDE" and not (re.search(r"tidal.{0,10}disrupt",title,re.I) and
                                           re.search(r"eROSITA|SRG",title,re.I)):
                    continue
                seendois.add(doi)
                urls = [rec.get("resource",{}).get("primary",{}).get("URL","")]
                urls.extend(x.get("URL","") for x in rec.get("link",[]))
                urls = [re.sub("^http://","https://",x) for x in dict.fromkeys(urls) if x.startswith(("http://","https://"))]
                if doi.lower()=="10.1051/0004-6361/201628241":
                    urls.append("https://www.aanda.org/articles/aa/full_html/2016/08/aa28241-16/aa28241-16.html")
                visited, seenids = set(), {}
                for link in urls:
                    if link in visited:
                        continue
                    visited.add(link)
                    try:
                        page = audit.get(doi,link)
                        if page.startswith(b"%PDF"):
                            run = subprocess.run(["pdftotext","-layout","-","-"],input=page,
                                                 capture_output=True,timeout=45,check=True)
                            body = run.stdout.decode("utf-8","replace")
                        else:
                            doc = ArticleHTML()
                            doc.feed(page.decode("utf-8","replace"))
                            body = " ".join(doc.parts)
                            # Full-text and stand-alone table URLs from the
                            # publisher's own article, excluding citation links.
                            for href in doc.links:
                                stem = re.search(r"(aa[0-9]+[-_][0-9]+)",link,re.I)
                                if (stem and stem.group().lower() in href.lower() and
                                    re.search(r"(_T\d+\.html|/full_html/|/pdf/)",href,re.I)) :
                                    dest = urljoin(link,href).split("#")[0]
                                    if urlparse(dest).hostname==urlparse(link).hostname and dest not in urls:
                                        urls.append(dest)
                        body = body.replace("−","-").replace("–","-")
                        for m in IAU_NAME.finditer(body):
                            prefix, a, sign, d = m.groups()
                            try:
                                ra,dc = parse_coord(a[:2]+":"+a[2:4]+":"+a[4:],
                                                    sign+d[:2]+":"+d[2:4]+":"+d[4:])
                            except ValueError:
                                continue
                            alias = f"{prefix} J{a}{sign}{d}"
                            seenids[alias] = {"name":alias,"ra":ra,"dec":dc,"doi":doi,"title":title}
                    except Exception as exc:
                        audit.record({"service":"paper-full-text","doi":doi,"url":link,"error":clipped_failure(exc)})
                if seenids:
                    searched.add("slew" if doi in PRIOR_PAPERS else kind)
                    mapped = match_fields(points,list(seenids.values()),(("ra",),("dec",)))
                    for k,hits in mapped.items():
                        matches[k].extend(hits)
                audit.record({"service":"paper-full-text","doi":doi,"title":title,
                              "distinct_coordinate_names":len(seenids),"urls_checked":sorted(visited)})
        except Exception as exc:
            audit.record({"service":"Crossref","query":phrase,"error":clipped_failure(exc)})
    return matches, searched


def byte_documents(raw, name, depth=0):
    """Read journal attachments and arXiv source in memory without extraction."""
    if depth > 3:
        return
    if raw.startswith(b"%PDF"):
        r = subprocess.run(["pdftotext","-layout","-","-"],input=raw,capture_output=True,
                           check=True,timeout=45)
        yield name,r.stdout.decode("utf-8","replace")
        return
    if zipfile.is_zipfile(io.BytesIO(raw)):
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            for f in z.infolist():
                if not f.is_dir() and f.file_size<=40_000_000:
                    yield from byte_documents(z.read(f),name+"/"+f.filename,depth+1)
        return
    if raw.startswith(b"\x1f\x8b") or raw[257:262]==b"ustar":
        try:
            with tarfile.open(fileobj=io.BytesIO(raw),mode="r:*") as archive:
                for m in archive:
                    if m.isfile() and m.size<=40_000_000:
                        yield from byte_documents(archive.extractfile(m).read(),name+"/"+m.name,depth+1)
            return
        except tarfile.TarError:
            if raw.startswith(b"\x1f\x8b"):
                yield from byte_documents(gzip.decompress(raw),name,depth+1)
                return
    if raw.startswith(b"SIMPLE  ="):
        from astropy.table import Table
        tab = Table.read(io.BytesIO(raw),format="fits")
        yield name,"\n".join(" ".join(str(v) for v in row) for row in tab)
        return
    if b"\x00" not in raw[:4096] and not name.lower().endswith((".eps",".ps",".jpg",".png",".bib")):
        yield name,raw.decode("utf-8","replace")


def coordinate_names(body):
    rows = {}
    body = body.replace("−","-").replace("–","-").replace("$","").replace("\,"," ").replace("~"," ")
    body = re.sub(r"(XMMSL[12]?)\s*\\[a-z]+\s*",r"\1 ",body,flags=re.I)
    for m in IAU_NAME.finditer(body):
        prefix, a, sign, d = m.groups()
        if not prefix.upper().startswith("XMMSL"):
            continue
        try:
            ra,dec = parse_coord(a[:2]+":"+a[2:4]+":"+a[4:],sign+d[:2]+":"+d[2:4]+":"+d[4:])
        except ValueError:
            continue
        sid = f"{a}{sign}{d}"
        rows[sid]={"name":prefix+" J"+sid,"ra":ra,"dec":dec}
    return rows


def li_supplement(audit,points,rec):
    """Search exact-DOI supplements. Require the published sample size for a clean non-hit."""
    doi="10.1093/mnras/stac714"
    matches, opened, expected, samples = defaultdict(list),set(),set(),{}
    candidates = [(r.get("URL",""),False) for r in rec.get("link",[])]
    candidates.append((rec.get("resource",{}).get("primary",{}).get("URL",""),False))
    abstract = re.sub("<[^>]+>"," ",rec.get("abstract",""))
    arxiv_ids=set()
    try:
        payload=json.loads(audit.get("Li-SemanticScholar","https://api.semanticscholar.org/graph/v1/paper/DOI:"+doi,
                           params={"fields":"externalIds,openAccessPdf,title,abstract"}))
        if payload.get("externalIds",{}).get("ArXiv"):
            arxiv_ids.add(payload["externalIds"]["ArXiv"])
        candidates.append(((payload.get("openAccessPdf") or {}).get("url",""),False))
        abstract+=" "+(payload.get("abstract") or "")
    except Exception as exc:
        audit.record({"service":"Li-SemanticScholar","error":clipped_failure(exc)})
    try:
        meta=json.loads(audit.get("Li-OpenAlex","https://api.openalex.org/works/https://doi.org/"+doi))
        for loc in meta.get("locations",[]):
            for url in (loc.get("landing_page_url"),loc.get("pdf_url")):
                if url:
                    if "arxiv.org/" in url:
                        m=re.search(r"/(?:pdf|abs)/([0-9]{4}\.[0-9]{4,5})",url)
                        if m:
                            arxiv_ids.add(m.group(1))
                    candidates.append((url,False))
    except Exception as exc:
        audit.record({"service":"Li-OpenAlex","error":clipped_failure(exc)})
    for query in ('doi:"'+doi+'"','ti:"Populations of highly variable X-ray sources"'):
        try:
            cat=ET.fromstring(audit.get("Li-arXiv-metadata","https://export.arxiv.org/api/query",
                          params={"search_query":query,"max_results":5}))
            for entry in cat.findall("{*}entry"):
                title=re.sub(r"\s+"," ",entry.findtext("{*}title","")).lower()
                if "slew" not in title or "variable" not in title:
                    continue
                match=re.search(r"([0-9]{4}\.[0-9]{4,5})",entry.findtext("{*}id",""))
                if match:
                    arxiv_ids.add(match.group(1))
                    abstract+=" "+entry.findtext("{*}summary","")
        except Exception as exc:
            audit.record({"service":"Li-arXiv-metadata","query":query,"error":clipped_failure(exc)})
    for ident in arxiv_ids:
        candidates.extend([("https://arxiv.org/src/"+ident,True),("https://arxiv.org/pdf/"+ident,False),
                           ("https://ar5iv.labs.arxiv.org/html/"+ident,False)])
    candidates.extend([("https://articles.adsabs.harvard.edu/pdf/2022MNRAS.512.3858L",False),
                       ("https://ui.adsabs.harvard.edu/link_gateway/2022MNRAS.512.3858L/EPRINT_PDF",False)])
    try:
        found=json.loads(audit.get("Li-figshare-search","https://api.figshare.com/v2/articles",
                                   params={"resource_doi":doi,"page_size":100}))
        for item in found:
            metadata=json.loads(audit.get("Li-figshare-record","https://api.figshare.com/v2/articles/"+str(item["id"])))
            if doi in str(metadata.get("resource_doi","")).lower() or doi in " ".join(metadata.get("references",[])).lower():
                candidates.extend((f.get("downloadUrl",""),True) for f in metadata.get("files",[])
                                  if 0<f.get("size",0)<=40_000_000)
    except Exception as exc:
        audit.record({"service":"Li-figshare","error":clipped_failure(exc)})
    try:
        meta=json.loads(audit.get("Li-DataCite","https://api.datacite.org/dois",params={
            "query":'relatedIdentifiers.relatedIdentifier:"'+doi+'"',"page[size]":80}))
        for item in meta.get("data",[]):
            fields=item["attributes"]
            related=[str(x.get("relatedIdentifier","")).lower() for x in fields.get("relatedIdentifiers",[])]
            if doi not in related:
                continue
            content=fields.get("contentUrl") or []
            candidates.extend((u,True) for u in content if isinstance(u,str))
            if fields.get("url"):
                candidates.append((fields["url"],True))
            audit.record({"service":"Li-related-deposit","doi":fields.get("doi"),"titles":fields.get("titles"),
                          "content_url":content,"landing":fields.get("url")})
    except Exception as exc:
        audit.record({"service":"Li-DataCite","error":clipped_failure(exc)})
    # Exact related DOI and record metadata bind any Zenodo files to the paper.
    try:
        meta=json.loads(audit.get("Li-Zenodo","https://zenodo.org/api/records",params={"q":'"'+doi+'"',"size":25}))
        for rec0 in meta["hits"]["hits"]:
            if doi.lower() in json.dumps(rec0.get("metadata",{})).lower():
                for f in rec0.get("files",[]):
                    if 0 < f.get("size",0) <= 40_000_000:
                        candidates.append((f.get("links",{}).get("self",""),True))
    except Exception as exc:
        audit.record({"service":"Li-Zenodo","error":clipped_failure(exc)})
    # Only explicit total-sample statements are used to corroborate completeness.
    total_pattern=re.compile(r"(?:identified|selected|obtained|found|sample\s+(?:comprises|includes|contains|of))"
            r"\s+(?:a\s+(?:total|sample)\s+of\s+)?([0-9][0-9,]+)\s+(?:highly\s+)?"
            r"(?:variable\s+)?(?:X.ray\s+)?sources",re.I)
    def collect_total(txt):
        for m in total_pattern.finditer(re.sub(r"\s+"," ",txt)):
            expected.add(int(m.group(1).replace(",","")))
    collect_total(abstract)
    for url,is_attachment in candidates:
        if not url:
            continue
        url=re.sub(r"^http://","https://",url)
        if not url.startswith("https://") or url in opened:
            continue
        opened.add(url)
        try:
            raw=audit.get("Li-source-material",url)
            for doc,body in byte_documents(raw,url):
                if body.lstrip().lower().startswith(("<html","<!doctype")):
                    html=ArticleHTML()
                    html.feed(body)
                    for href in html.links:
                        if re.search(r"suppl|stac714[_-]s|supporting.information",href,re.I) or (
                               is_attachment and re.search(r"\.(?:fits|fit|csv|tsv|dat|txt|zip)(?:\?|$)",href,re.I)):
                            dest=urljoin(url,href).split("#")[0]
                            if dest!=url and (dest,True) not in candidates:
                                candidates.append((dest,True))
                    body=" ".join(html.parts)
                # Published data availability statements sometimes link the
                # author's deposit without an electronic-journal attachment.
                before_bib=re.split(r"\\begin\{thebibliography\}|\nREFERENCES\s*\n",body,flags=re.I)[0]
                for match in re.finditer(r"https?://[^\s<>{}\\]+",before_bib):
                    link=match.group().rstrip(".,);]$'"+chr(34))
                    host=urlparse(link).hostname or ""
                    note=before_bib[max(0,match.start()-180):match.end()+120]
                    if (any(h in host for h in ("zenodo.org","figshare.com","china-vo.org","nadc.china",
                           "cdsarc","vizier","github.com","heasarc.gsfc")) and re.search(
                           r"supplement|machine.readable|data availab|electronic|full table|full catalog",note,re.I)
                           and (link,True) not in candidates):
                        candidates.append((link,True))
                names=coordinate_names(body)
                if re.search(r"XMMSL[12]",body,re.I) and (is_attachment or re.search(r"tab|suppl",doc,re.I)):
                    # Rows often omit the survey prefix printed in the
                    # column heading. Keep a leading J and both coordinates.
                    for m in re.finditer(r"(?m)(?:^|[|&>\s])J([0-9]{6}(?:\.[0-9]+)?[+−–-][0-9]{6}(?:\.[0-9]+)?)",body):
                        names.update(coordinate_names("XMMSL J"+m.group(1)))
                # A primary article may count its own sample in its abstract
                # and provide a long main/appendix table. Retain sample counts
                # by document and require independent agreement with total.
                if doi in body or "Populations of highly variable" in body:
                    collect_total(body[:18000])
                if names:
                    samples[doc]=names
            audit.record({"service":"Li-source-material","url":url,"is_attachment":is_attachment})
        except Exception as exc:
            audit.record({"service":"Li-source-material","url":url,"error":clipped_failure(exc)})
    allrows={}
    for name,data in samples.items():
        for key,row in data.items():
            allrows[key]=dict(row,table="doi:"+doi,document=name)
    # Illustrative names from unrelated sections cannot establish negative
    # coverage. Require one data/TeX table with the reported number of sources.
    eligible = [(doc,len(rows)) for doc,rows in samples.items() if
                re.search(r"tab|suppl|\.csv|\.dat|\.fits|\.tsv",doc,re.I)]
    complete = any(n>=20 and n in expected for _,n in eligible)
    audit.record({"service":"Li-exact-coverage","doi":doi,"arxiv_ids":sorted(arxiv_ids),
                  "expected_counts_from_primary":sorted(expected),"unique_slew_source_names":len(allrows),
                  "table_counts":eligible,"complete":complete,"abstract":abstract[:4800]})
    matched=match_fields(points,list(allrows.values()),(("ra",),("dec",)))
    return matched,complete


def published_catalogue_matches(audit, points, exact):
    result, inventory = defaultdict(list), []
    (viz,mirror), catalogues = publication_catalogs(audit,exact)
    for cat, terms in sorted(catalogues.items()):
        try:
            try:
                table_list = viz.get_catalogs(cat)
                if not len(table_list):
                    raise ValueError("empty primary resource")
            except Exception:
                table_list = mirror.get_catalogs(cat)
            if not len(table_list):
                raise ValueError("empty resource or service failure")
            for table_name in table_list.keys():
                tab = table_list[table_name]
                item = {"catalog":cat, "table":table_name, "terms":sorted(terms), "rows":len(tab)}
                inventory.append(item)
                if len(tab) >= 20001:
                    item["status"] = "truncated; excluded"
                    continue
                ra, dec = positional_columns(tab)
                if not ra or not dec:
                    item["columns"] = list(tab.colnames)
                    if terms & {"Li2022","Strotjohann","Starling","slew transients","slew variability"}:
                        identified=[]
                        for idx,row in enumerate(tab):
                            named={}
                            for field in tab.colnames:
                                value=str(row[field])
                                if "xmmsl" in field.lower() and "xmmsl" not in value.lower():
                                    value="XMMSL "+value
                                named.update(coordinate_names(value))
                            identified.extend(dict(obj,table=table_name,catalog=cat,record=idx+1)
                                              for obj in named.values())
                        if identified:
                            hits=match_fields(points,identified,(("ra",),("dec",)))
                            for k,rs in hits.items():
                                result[k].extend(rs)
                            item.update(status="coordinate names",usable=len(identified),rows_sha256=hash_table(tab))
                            continue
                    item["status"] = "table lacks paired sky coordinates"
                    continue
                # Hash the exact table projection used for the crossmatch.
                item["rows_sha256"] = hash_table(tab)
                objects, failures = [], 0
                idfields = [x for x in tab.colnames if x.lower() in {
                    "name","objid","id","source","sourceid","tns","tnsname","simbad","iau",
                    "xmmsl1","xmmsl2","rosat","erosita","object","atname","galaxy","event"}]
                for idx, r in enumerate(tab):
                    try:
                        x,y = parse_coord(r[ra],r[dec])
                    except (ValueError,TypeError):
                        failures += 1
                        continue
                    ids = " | ".join(str(r[k]) for k in idfields if str(r[k]) != "--")
                    objects.append({"ra":x, "dec":y, "catalog":cat, "table":table_name,
                                    "record":idx+1, "name":ids, "search_terms":sorted(terms)})
                matched = match_fields(points, objects, (("ra",),("dec",)))
                for k,v in matched.items():
                    result[k].extend(v)
                item.update(status="ok", usable=len(objects), invalid_coordinates=failures)
        except Exception as exc:
            inventory.append({"catalog":cat, "terms":sorted(terms), "status":clipped_failure(exc)})
    audit.record({"service":"VizieR-publication-tables", "tables":inventory})
    # A TDE sample and a slew variability sample must both have searchable
    # positions. Other omissions are reported in the full coverage inventory.
    complete = {t for i in inventory if i.get("status")=="ok" and i.get("usable",0)
                and not i.get("invalid_coordinates") for t in i["terms"]}
    return result,complete


def publications(audit,points):
    exact=exact_papers(audit)
    doi="10.1093/mnras/stac714"
    pieces={
        "published_catalogues":(published_catalogue_matches,(audit,points,exact)),
        "open_primary_papers":(article_objects,(audit,points,exact)),
        "Li_exact_electronic_sources":(li_supplement,(audit,points,exact.get(doi,({},))[0])),
    }
    result=defaultdict(list)
    coverage={}
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures={key:pool.submit(fn,*args) for key,(fn,args) in pieces.items()}
        for key,fut in futures.items():
            try:
                hits,complete=fut.result()
                coverage[key]=complete
                for k,rows in hits.items():
                    for h in rows:
                        if key=="open_primary_papers":
                            h=dict(h,table="doi:"+h["doi"])
                        result[k].append(h)
            except Exception as exc:
                audit.record({"service":key,"task_error":clipped_failure(exc),"traceback":traceback.format_exc()})
    catalogue=coverage.get("published_catalogues",set())
    paper=coverage.get("open_primary_papers",set())
    missing=[]
    if "TDE" not in catalogue and "TDE" not in paper:
        missing.append("TDE sample")
    if not (catalogue & {"Strotjohann","slew variability"}) and "slew" not in paper:
        missing.append("slew variability sample")
    if "Li2022" not in catalogue and not coverage.get("Li_exact_electronic_sources",False):
        missing.append("Li2022 electronic source list")
    state="ok" if not missing else "incomplete "+", ".join(missing)
    audit.record({"service":"published-coverage-summary","status":state,"catalogue_samples":sorted(catalogue),
                  "paper_samples":sorted(paper),"Li_supplement_complete":coverage.get("Li_exact_electronic_sources",False)})
    return result,{p["k"]:state for p in points}



def compressed_json(path, records):
    with path.open("wb") as stream, gzip.GzipFile(fileobj=stream, mode="wb", mtime=0) as gz:
        for record in records:
            gz.write((json.dumps(record, sort_keys=True, default=str) + "\n").encode())


def write_table(path, rows, columns):
    with path.open("w",newline="",encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def is_galaxy(name):
    return bool(re.match(r"^(G|Galaxy|AGN|QSO|Sy1|Sy2|Seyfert|LINER|Bla|BLLac|EmG|GiP|GiG|rG|RadioG|SBG)\??$",name,re.I))


def extended(source):
    return str(source.get("type","")).strip().upper() in ("EXP","DEV","SER","REX","COMP")


def distance(first, second):
    try:
        a=SkyCoord(float(first["ra"])*u.deg,float(first["dec"])*u.deg)
        b=SkyCoord(float(second["ra"])*u.deg,float(second["dec"])*u.deg)
        return a.separation(b).arcsec
    except (ValueError,TypeError,KeyError):
        return math.inf


def host_summary(found):
    candidates = []
    for key, ident, typ in (("SIMBAD","main_id","otype"),("NED","preferred_name","type")):
        for r in found[key]:
            if is_galaxy(scalar(r,typ)):
                candidates.append((scalar(r,ident), scalar(r,typ), key, r))
    for r in found["Legacy"]:
        if extended(r):
            candidates.append((scalar(r,"ls_id"),scalar(r,"type"),"Legacy",r))
    if not candidates:
        for key in ("SIMBAD", "NED", "Gaia", "Legacy", "WISE"):
            fields = {"SIMBAD":("main_id","otype"),"NED":("preferred_name","type"),
                      "Gaia":("id",""),"Legacy":("ls_id","type"),"WISE":("id","")}
            if found[key]:
                r = found[key][0]
                f,t = fields[key]
                return (scalar(r,f), scalar(r,t), key, r["sep_arcsec"], False, "unresolved host class")
        return ("", "", "", "", False, "no optical/IR counterpart identified")
    name,typ,key,first = candidates[0]
    rivals = [r for _,_,_,r in candidates[1:] if distance(first,r)>2.0]
    # Unclassified/pointlike optical alternatives also matter.
    rivals.extend(r for r in found["Legacy"] + found["Gaia"] if distance(first,r)>2.0)
    # A secure host needs an independently identified nucleus/generic galaxy,
    # or an extended Legacy source independently coincident with AllWISE.
    independent = key in ("NED","SIMBAD") or any(distance(first,r)<=2 for r in found["WISE"])
    reason = "competing optical counterparts" if rivals else (
        "one extended source without independent identification" if not independent else "")
    return (name,typ,key,first["sep_arcsec"],not reason,reason)

MPE_BASE = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/"
MPE_LISTS = {
    "HamStar-likely":"FreundS_DR1/HamStar_eRASS1_Main_Likely_Identifications_v1.1.fits.tgz",
    "HamStar-possible":"FreundS_DR1/HamStar_eRASS1_Main_Possible_Identifications_v1.1.fits.tgz",
    "BlazEr1":"HaemmerichS_DR1/BlazEr1_catalog.fits",
    "BLAZE":"HaemmerichS_DR1/BLAZE_catalog.fits",
    "HMXB":"KaltenbrunnerD_DR1/eRASS1_HMXB_LMC_v1.0.fits.tgz",
    "ULX":"WeberP_DR1/ulx_erass1_main.fits",
    "CV":"SchwopeA_DR1/allnewdr1cvs_fin_2cds.fits",
}


def ancillary(audit, points):
    import tarfile
    from astropy.io import fits
    result, states = defaultdict(list), {}
    if not points:
        return result, states
    centers = sky(points)
    for cat, path in MPE_LISTS.items():
        try:
            raw = audit.get(cat, MPE_BASE + path)
            if path.endswith(".tgz"):
                try:
                    with tarfile.open(fileobj=io.BytesIO(raw),mode="r:*") as tf:
                        options=[m for m in tf if m.isfile() and ".fits" in m.name.lower()]
                        if not options:
                            raise ValueError("archive contains no FITS members")
                        raw=tf.extractfile(max(options,key=lambda m:m.size)).read()
                except tarfile.TarError:
                    raw=gzip.decompress(raw)
            if raw.startswith(b"\x1f\x8b"):
                raw=gzip.decompress(raw)
            with fits.open(io.BytesIO(raw),memmap=False) as hdus:
                hdu=next(h for h in hdus if getattr(h,"columns",None) and h.data is not None)
                tab=hdu.data
                cols={re.sub(r"[^a-z0-9]","",n.lower()):n for n in tab.names}
                audit.record({"service":cat+"-schema","columns":list(tab.names),"units":{
                    c.name:c.unit for c in hdu.columns if c.unit}})
                pairs=[("erora","erodec"),("raero","decero"),("raerass1","decerass1"),
                       ("ra","dec"),("radeg","decdeg"),("rax","decx"),("raxray","decxray"),
                       ("racorr","deccorr"),("raerosita","decerosita"),
                       ("raj2000","dej2000"),("raopt","decopt"),("raicrs","deicrs"),
                       ("ragdr3","decgdr3"),("gaiara","gaiadec"),("raj2000","decj2000")]
                if cat.startswith("HamStar"):
                    pairs.insert(0,("ctpra","ctpdec"))
                for key in cols:
                    if key.startswith("ra"):
                        pairs.append((key,"dec"+key[2:]))
                        pairs.append((key,"de"+key[2:]))
                    if key.endswith("ra"):
                        pairs.append((key,key[:-2]+"dec"))
                selected=next(((cols[a],cols[b]) for a,b in pairs if a in cols and b in cols),None)
                if not selected:
                    # Some counterpart releases contain only the DR1 source
                    # designation and the Gaia donor ID. Use the published join.
                    def normalized_id(v):
                        hit=re.search(r"J[0-9]{6}(?:\.[0-9]+)?[+-][0-9]{6}(?:\.[0-9]+)?",str(v),re.I)
                        return hit.group().upper() if hit else ""
                    wanted=defaultdict(list)
                    for p in points:
                        if normalized_id(p.get("ero_name","")) and p["ero_ra"] is not None and p["ero_dec"] is not None:
                            wanted[normalized_id(p["ero_name"])].append(p)
                    named=[cols[n] for n in ("eroname","erass1","erass1name","iauname","eroiauname",
                             "erass1iauname","name","srcname","xname") if n in cols]
                    key=next((n for n in named if any(normalized_id(v) for v in np.asarray(tab[n]).ravel()[:20])),None)
                    if not key:
                        raise ValueError("no paired sky coordinates or eRASS1 designation in "+str(tab.names))
                    ids=[cols[n] for n in ("gdr3","sourceid","gaiasourceid","gaiadr3","gaiadr3sourceid",
                                           "name","id","xid") if n in cols]
                    for j,v in enumerate(np.asarray(tab[key]).ravel()):
                        for p in wanted.get(normalized_id(v),[]):
                            er={"ra":p["ero_ra"],"dec":p["ero_dec"]}
                            result[p["k"]].append(dict(catalog=cat,row=j+1,id=" | ".join(str(
                                np.asarray(tab[n]).ravel()[j]) for n in ids),ra=er["ra"],dec=er["dec"],
                                sep_arcsec=round(distance(er,p),3),method="exact DR1 designation"))
                    states[cat]="ok"
                    audit.record({"service":cat,"source_name_column":key,"rows":len(np.asarray(tab[key]).ravel())})
                    continue
                ra,dec=selected
                xx, yy = np.asarray(tab[ra], dtype=float).ravel(), np.asarray(tab[dec],dtype=float).ravel()
                mask = np.isfinite(xx) & np.isfinite(yy) & (np.abs(yy)<=90) & (xx>=0) & (xx<360)
                if not np.any(mask):
                    raise ValueError("zero usable positions")
                flat = np.flatnonzero(mask)
                catpos = SkyCoord(xx[mask]*u.deg,yy[mask]*u.deg)
                a,b,sep,_ = search_around_sky(centers,catpos,SEARCH_RADIUS*u.arcsec)
                names = [cols[n] for n in ("eroiauname","ctpid","simbadname","iauname","name","detuid","erass1","designation","source_id",
                                             "srcid","object","gdr3","gaiadr3","gaiadr3sourceid","gaiasourceid","gaiaid","sourceid","xid") if n in cols]
                for i,j,d in zip(a,b,sep.arcsec):
                    idx = flat[j]
                    ids = [str(np.asarray(tab[n]).ravel()[idx]).strip() for n in names]
                    hit={"catalog":cat,"id":" | ".join(ids),"row":int(idx+1),
                         "ra":xx[idx],"dec":yy[idx],"sep_arcsec":round(d,3)}
                    if cat.startswith("HamStar"):
                        for field in ("ERO_IAUNAME","CTP_ID","CTP_RANK","p_coronal","p_ij",
                                      "PLX","SIMBAD_NAME","SIMBAD_OTYPE","CORONAL","FLAG_OPT"):
                            col=cols.get(re.sub("[^a-z0-9]","",field.lower()))
                            if col:
                                val=np.asarray(tab[col]).ravel()[idx]
                                hit[field]=(val.decode("utf-8","replace").strip() if isinstance(val,bytes) else str(val).strip())
                        base=lambda value: str(value).lower().split("j")[-1].replace(" ","")
                        ename=hit.get("ERO_IAUNAME")
                        match=bool(ename and points[i].get("ero_name") and base(ename)==base(points[i]["ero_name"]))
                        hit["same_ero_name"]=match
                        hit["method"]="proper-motion corrected optical position at eRASS1 epoch"
                    result[points[i]["k"]].append(hit)
                states[cat] = "ok"
                audit.record({"service":cat,"ra_col":ra,"dec_col":dec,
                              "valid_rows":int(mask.sum()), "total_rows":len(xx)})
        except Exception as exc:
            states[cat] = clipped_failure(exc) or type(exc).__name__
    audit.record({"service":"MPE-ancillary-coverage", "coverage":states})
    return result, {p["k"]:"ok" if all(s=="ok" for s in states.values()) else
                             "unsearched: "+",".join(c for c,s in states.items() if s!="ok")
                    for p in points}


def confident_coronal(row):
    """Require one published coronal identification for this precise eROSITA source."""
    return (row["catalog"]=="HamStar-likely" and row.get("same_ero_name") and
            number(row,"CTP_RANK")==1 and (number(row,"p_coronal") or 0)>=0.8 and
            (number(row,"p_ij") or 0)>=0.8)


def vet_rows(rows, output_dir):
    """Return one verdict per input row. Write all positional alternatives and service outcomes."""
    rows = list(rows)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    columns = ["id","ra","dec","match_radius_arcsec","vet_status","host_id","host_type","host_catalog",
               "host_offset_arcsec","rationale","known_associations","prior_bibcodes",
               "prior_identity","stellar_evidence","event_review","bibliography_context",
               "shifted_host_controls","service_outcomes"]
    if not rows:
        write_table(out/"vetting.csv",[],columns)
        return []
    points = select_points(rows, controls=True)
    real = [p for p in points if not p["shift"]]
    audit = Audit()
    tasks = {
        "Gaia":(verified_service,(audit,"Gaia",points)),
        "WISE":(verified_service,(audit,"WISE",points)),
        "VSX":(verified_service,(audit,"VSX",real)),
        "Legacy":(region_query,(audit,"Legacy",points,legacy_query,(("ra",),("dec",)))),
        "NED":(region_query,(audit,"NED",real,ned_query,(("ra",),("dec",)))),
        "SIMBAD":(verified_service,(audit,"SIMBAD",real)),
        "papers":(publications,(audit,real)),
        "MPE":(ancillary,(audit,real)),
    }
    found, status = {}, {}
    with ThreadPoolExecutor(max_workers=len(tasks)) as executor:
        futures = {k:executor.submit(fn,*args) for k,(fn,args) in tasks.items()}
        for k,fut in futures.items():
            try:
                found[k],status[k] = fut.result()
            except Exception as exc:
                audit.record({"service":k,"task_error":clipped_failure(exc),"traceback":traceback.format_exc()})
                found[k] = {}
                status[k] = {p["k"]:clipped_failure(exc) or type(exc).__name__ for p in points}
    previously_classified = {}
    for p in real:
        k = p["k"]
        near = [h for h in found["MPE"].get(k,[]) if h["sep_arcsec"]<=p["gate"]]
        binary = any(r["catalog"] in ("CV","HMXB","ULX") for r in near) or any(
            TRANSIENT.fullmatch(r.get("otype","")) for r in found["SIMBAD"].get(k,[]) if r["sep_arcsec"]<=p["gate"])
        previously_classified[k] = not binary and any(confident_coronal(r) for r in near)
    found["TNS"],status["TNS"] = tns(audit,real,previously_classified)
    refs, ref_status = associate_refs(audit, found["SIMBAD"], real)
    output, evidence, controls = [], [], []
    for p in real:
        k = p["k"]
        source = rows[p["i"]]
        sid = str(source["id"])
        in_gate = {name:sorted((h for h in hits.get(k,[]) if h["sep_arcsec"]<=p["gate"]),
                               key=lambda h:h["sep_arcsec"]) for name,hits in found.items()}
        checks = {name:states.get(k,"not searched") for name,states in status.items()}
        sim_ids = [int(h["oid"]) for h in in_gate["SIMBAD"] if h.get("oid")]
        checks["SIMBAD_bibliography"] = ("ok" if checks["SIMBAD"]=="ok" and
                len(sim_ids)==len(in_gate["SIMBAD"]) and all(ref_status.get(n)=="ok" for n in sim_ids) else "incomplete")
        citations = [ref for oid in sim_ids for ref in refs.get(oid,[])]
        host,typ,cat,sep,secure,host_issue = host_summary(in_gate)
        known, uncertain, identity, event, context = [], [], [], [], []
        xtype = scalar(source,"xmm_xtyp")
        xident = scalar(source,"xmm_counterpart","xmm_altident")
        xidsep = number(source,"xmm_iddist")
        binary_re = re.compile(r"^(HXB|LXB|XB|XRB|HMXB|LMXB|ULX|CV\*)$|x.ray binary|cataclysmic",re.I)
        variablestar_re = re.compile(r"^(V\*|Fl\*|Er\*|RS\*|BY\*|RR\*|Mira|PulsV\*)$|flare star|variable star",re.I)
        binary = any(binary_re.search(h["otype"]) for h in in_gate["SIMBAD"]) or any(
                    h["catalog"] in ("CV","HMXB","ULX") for h in in_gate["MPE"])
        xmm_shift = distance({"ra":scalar(source,"xmm_ra") or p["ra"],"dec":scalar(source,"xmm_dec") or p["dec"]},p)
        if xident and xtype and xidsep is not None and xidsep >= 0:
            identity.append("slew "+xident+" "+xtype+f" ({xidsep:g} arcsec from XMM)")
            if xidsep+xmm_shift <= p["gate"]:
                binary = binary or bool(binary_re.search(xtype))
                if not binary and "?" not in xtype and variablestar_re.search(xtype):
                    known.append("catalogued variable/flare star "+xident+" "+xtype)
                elif not binary and (re.search(r"\bstar\b|dwarf|stellar",xtype,re.I) or STELLAR.fullmatch(xtype)):
                    uncertain.append("possible stellar X-ray emitter "+xident+" "+xtype)
        for h in in_gate["SIMBAD"]:
            identity.append("SIMBAD "+h["main_id"]+" "+h["otype"])
            if not binary and variablestar_re.search(h["otype"]):
                known.append("catalogued SIMBAD variable/flare star "+h["main_id"]+" "+h["otype"])
            elif not binary and (STELLAR.fullmatch(h["otype"]) or "*" in h["otype"]):
                uncertain.append("stellar SIMBAD counterpart "+h["main_id"])
            if re.fullmatch(r"(TDE|SN|Nova)\??",h["otype"],re.I):
                event.append("catalogued "+h["otype"]+" "+h["main_id"]+"; compare event epochs")
        for h in in_gate["MPE"]:
            tag = h["catalog"]+" "+(h["id"] or f"row {h['row']}")
            identity.append(tag)
            if not binary and confident_coronal(h):
                known.append("rank-one, exact eRASS1 coronal-star identification "+tag)
            elif not binary and h["catalog"].startswith("HamStar-"):
                uncertain.append("possible coronal-star identification "+tag)
        for h in in_gate["VSX"]:
            identity.append("VSX "+h["id"]+" "+h["type"])
            if not binary:
                uncertain.append("potential known variable "+h["id"])
        for h in in_gate["TNS"]:
            event.append("TNS "+h["name"]+" "+h["type"]+f" ({h['sep_arcsec']} arcsec); compare event epochs")
        for h in in_gate["papers"]:
            event.append(h["table"]+" "+(h["name"] or f"row {h['record']}"))
        for h in in_gate["Gaia"]:
            if stellar_astrometry(h):
                tag = "significant Gaia stellar astrometry "+h["id"]
                identity.append(tag)
                if not binary:
                    uncertain.append(tag)
                    if any(distance(h,s)<=2 for s in in_gate["VSX"]):
                        known.append("astrometric star coincident with VSX "+h["id"])
            if scalar(h, "variable").upper() == "VARIABLE" and not binary:
                uncertain.append("Gaia variable "+h["id"])
        # For a classified X-ray binary an isolated optical stellar counterpart
        # can be a donor; stellar astrometry alone cannot reject this hypothesis.
        if binary and not secure:
            opt = in_gate["Gaia"] or in_gate["Legacy"]
            if len(opt)==1 and all(distance(opt[0],h)<=2 for h in in_gate["Gaia"]+in_gate["Legacy"]):
                cat = "Gaia" if in_gate["Gaia"] else "Legacy"
                host,typ,sep = scalar(opt[0],"id","ls_id"),"possible donor in classified X-ray binary",opt[0]["sep_arcsec"]
                secure,host_issue = bool(host),""
        identifiers = [host,xident,scalar(source,"name")]+[h["main_id"] for h in in_gate["SIMBAD"]]
        norm = lambda x: re.sub(r"[^a-z0-9]","",x.lower())
        aliases = [norm(x) for x in identifiers if x and len(norm(x))>=6]
        for ref in citations:
            if any(x in norm(ref["title"]) for x in aliases):
                event.append("object-named variability paper "+ref["bibcode"])
            else:
                context.append("host bibliography "+ref["bibcode"]+" (title has no matched alias)")
        uncertain.extend(event)

        ctl = {}
        for label,catalog,qualifies in (("Gaia","Gaia",bool),("Legacy","Legacy",bool),
                    ("WISE","WISE",bool),("Gaia_significant_astrometry","Gaia",stellar_astrometry),
                    ("Legacy_extended","Legacy",extended)):
            flags = [any(x["sep_arcsec"]<=p["gate"] and qualifies(x) for x in found[catalog].get(k+j,[]))
                     for j in range(1,5) if status[catalog].get(k+j)=="ok"]
            ctl[label] = [sum(flags),len(flags)]
            controls.append({"id":sid,"catalog":label,"shifted_with_counterpart":sum(flags),
                             "shifted_checked":len(flags),"real_with_counterpart":int(any(qualifies(x) for x in in_gate[catalog])),
                             "real_checked":int(checks[catalog]=="ok"),"match_radius_arcsec":round(p["gate"],2)})
        incomplete = [n for n,s in checks.items() if s!="ok" and n!="MPE"]
        if known:
            outcome,why = "known_or_stellar","catalogued stellar variability or published coronal identification"
        elif uncertain or incomplete or not secure:
            outcome = "inconclusive"
            why = " | ".join(filter(None,[host_issue, "prior association or stellar/variable alternative" if uncertain else "",
                                          "incomplete "+",".join(incomplete) if incomplete else ""]))
        else:
            outcome,why = "candidate","optical/IR counterpart identified; no prior event in completed positional and object-named bibliography checks"
        output.append(dict(id=sid,ra=p["ra"],dec=p["dec"],match_radius_arcsec=round(p["gate"],2),
                           vet_status=outcome,host_id=host,host_type=typ,host_catalog=cat,host_offset_arcsec=sep,
                           rationale=why,known_associations=" | ".join(dict.fromkeys(known+uncertain)),
                           prior_bibcodes=" | ".join(sorted({r["bibcode"] for r in citations})),
                           prior_identity=" | ".join(dict.fromkeys(identity)),stellar_evidence=" | ".join(known),
                           event_review=" | ".join(dict.fromkeys(event)),bibliography_context=" | ".join(context),
                           shifted_host_controls=json.dumps(ctl,sort_keys=True),
                           service_outcomes=json.dumps(checks,sort_keys=True)))
        evidence.append({"id":sid,"ra":p["ra"],"dec":p["dec"],"gate_arcsec":p["gate"],
                         "alternatives_within_35arcsec":{n:h.get(k,[]) for n,h in found.items()},
                         "slew_identification":{"id":xident,"class":xtype,"separation_arcsec":xidsep},
                         "bibliography_hits":citations,"checks":checks})
    write_table(out/"vetting.csv",output,columns)
    write_table(out/"host_shifted_controls.csv",controls,list(controls[0]))
    compressed_json(out/"counterpart_evidence.jsonl.gz",evidence)
    compressed_json(out/"vetting_provenance.jsonl.gz",audit.events)
    (out/"publication_coverage.json").write_text(json.dumps([e for e in audit.events if
        e.get("service") in ("Li-exact-coverage","Li-related-deposit","published-coverage-summary",
        "VizieR-publication-tables","VizieR-publication-discovery","DOI-to-CDS") or e.get("task_error")],
        indent=2,default=str)+"\n")
    (out/"tns_replica_coverage.json").write_text(json.dumps([
        e for e in audit.events if e.get("service")=="TNS-replica-audit"],indent=2,default=str)+"\n")
    return output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input")
    parser.add_argument("output_dir")
    args = parser.parse_args()
    with open(args.input, newline="", encoding="utf-8") as inp:
        vet_rows(csv.DictReader(inp), args.output_dir)
