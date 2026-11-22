"""Anonymous catalog and independent photometry checks for the CHEOPS search.

Run only inside the data workflow. All network responses used in a scientific
decision have a URL (or archive query) and a SHA-256 digest. Time is BJD_TDB.
"""
import argparse
import csv
import copy
import gzip
import os
import hashlib
import io
import json
import math
import re
from urllib.parse import urlsplit
from collections import Counter
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import requests

NEA = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
VIZIER = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
SIMBAD = "https://simbad.cds.unistra.fr/simbad/sim-tap/sync"
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"
ARXIV = "https://export.arxiv.org/api/query"


# The run-2 US runner timed out connecting to the Strasbourg hosts. HTTP
# absence of evidence is always recorded separately from a successful empty cone.
ASU_MIRRORS = ("vizier.cfa.harvard.edu","vizier.nao.ac.jp","vizier.iucaa.in","vizier.cds.unistra.fr")
SIMBAD_MIRRORS = ("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync",SIMBAD)
OUTAGES = Counter()


def _route(key, operation):
    if OUTAGES[key] >= 2:
        raise ConnectionError(f"route unavailable after two previous failures: {key}")
    try:
        result = operation()
    except (requests.ConnectionError, requests.Timeout, requests.HTTPError) as exc:
        if not isinstance(exc,requests.HTTPError) or exc.response.status_code in (403,429,500,502,503,504):
            OUTAGES[key] += 1
        raise
    OUTAGES[key] = 0
    return result


def _hash(data):
    return hashlib.sha256(data).hexdigest()


def _filehash(path):
    with open(path, "rb") as src:
        return hashlib.file_digest(src, "sha256").hexdigest()


def _get(cache, url, params):
    query = requests.Request("GET", url, params=params).prepare().url
    folder = Path(cache) / "responses"
    folder.mkdir(parents=True, exist_ok=True)
    name = folder / _hash(query.encode())
    if not name.exists():
        def request():
            reply = requests.get(url, params=params, timeout=(4, 32),
                                 headers={"User-Agent": "Open-Space-CHEOPS/1.0"})
            reply.raise_for_status()
            return reply
        reply = _route(urlsplit(url).netloc, request)
        with tempfile.NamedTemporaryFile(dir=folder, delete=False) as temp:
            temp.write(reply.content)
        Path(temp.name).replace(name)
    data = name.read_bytes()
    return data, {"url": query, "sha256": _hash(data), "bytes": len(data)}


def _tap(cache, endpoint, sql):
    raw, source = _get(cache, endpoint, {
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": 200000, "QUERY": sql})
    decoded = raw.decode("utf-8-sig")
    if not decoded.strip() or decoded.lstrip().startswith("<"):
        raise ValueError(f"Expected TAP CSV; received {decoded[:160]}")
    reader = csv.DictReader(io.StringIO(decoded))
    source["columns"] = reader.fieldnames
    if not reader.fieldnames or any("exception" in col.lower()
                                     for col in reader.fieldnames):
        raise ValueError("Invalid TAP response")
    return list(reader), source


def _number(v):
    try:
        f = float(v)
        return f if math.isfinite(f) else None
    except (ValueError, TypeError):
        return None


def _uncertainty(row, field):
    return max(abs(_number(row.get(field + suffix)) or 0.0)
               for suffix in ("err1", "err2"))


def _norm(name):
    return re.sub(r"[^a-z0-9]", "", str(name).lower())


def _sep(ra1, dec1, ra2, dec2):
    if any(_number(v) is None for v in (ra1, dec1, ra2, dec2)):
        return float("inf")
    a, b, c, d = np.radians([float(v) for v in (ra1, dec1, ra2, dec2)])
    h = np.sin((b - d) / 2) ** 2 + np.cos(b) * np.cos(d) * np.sin((a - c) / 2) ** 2
    return float(np.degrees(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))) * 3600)


def _cone(ra, dec, radius_arcsec, racol, deccol):
    return ("1=CONTAINS(POINT('ICRS'," + racol + "," + deccol + "),"
            f"CIRCLE('ICRS',{float(ra):.9f},{float(dec):.9f},{radius_arcsec / 3600:.9f}))")


def _orbit_phase(ecc, omega):
    def mean_anomaly(f):
        eanom = 2 * math.atan2(math.sqrt(1 - ecc) * math.sin(f / 2),
                               math.sqrt(1 + ecc) * math.cos(f / 2))
        return eanom - ecc * math.sin(eanom)
    f1 = math.pi / 2 - math.radians(omega)
    return ((mean_anomaly(f1 + math.pi) - mean_anomaly(f1)) % (2 * math.pi)) / (2 * math.pi)


def _planet_phases(row, is_toi, unresolved):
    label = "TOI " + str(row.get("toi", "?")) if is_toi else row["pl_name"]
    t0 = _number(row.get("pl_tranmid"))
    p = _number(row.get("pl_orbper"))
    field = "pl_trandurh" if is_toi else "pl_trandur"
    dur = _number(row.get(field))
    if t0 is None or not 2400000 < t0 < 2600000:
        unresolved.append(label + " has no usable published transit epoch")
        return []
    if not dur or dur <= 0:
        unresolved.append(label + " has no published transit duration; 4 h masking window assumed")
        dur = 4
    if p is None or p <= 0:
        p = None
        unresolved.append(label + " has no published period")
    ep = dict(label=label, epoch_bjd=t0, period_days=p, duration_days=dur/24,
              epoch_error_days=_uncertainty(row, "pl_tranmid"),
              period_error_days=_uncertainty(row, "pl_orbper"),
              kind="published TOI" if is_toi else "planet transit")
    phases = [ep]
    if not p:
        return phases
    e = _number(row.get("pl_orbeccen"))
    w = _number(row.get("pl_orblper"))
    if e is None or not 0 <= e < 1 or (e and w is None):
        phases.append(dict(ep, kind="nominal half-period eclipse", epoch_bjd=t0+p/2))
        unresolved.append(label + " secondary phase is unconstrained; half-phase also masked")
        return phases
    phase = _orbit_phase(e, w or 0)
    # Project the published orbit errors into the occultation time. The
    # edge-on conjunction approximation still requires review near grazing orbits.
    de = _uncertainty(row, "pl_orbeccen")
    dw = _uncertainty(row, "pl_orblper")
    phase_err = max(abs(((_orbit_phase(ee, ww) - phase + .5) % 1) - .5)
                    for ee in (max(0, e-de), min(.999999, e+de), e)
                    for ww in ((w or 0) - dw, (w or 0) + dw, w or 0))
    factor = (1 + e * math.sin(math.radians(w or 0))) / (1 - e * math.sin(math.radians(w or 0)))
    phases.append(dict(ep, epoch_bjd=t0+p*phase, duration_days=dur*factor/24,
                       epoch_error_days=ep["epoch_error_days"] + p*phase_err,
                       kind="predicted planetary occultation"))
    return phases


def _binary_phases(row, unresolved):
    label = "TESS-EB TIC "+str(row["TIC"])
    p, t0 = _number(row.get("Per")), _number(row.get("BJD0"))
    if not p or p <= 0 or t0 is None:
        unresolved.append(label+" lacks an eclipse epoch or period")
        return []
    t0 += 2457000.0
    primary = dict(label=label, kind="published binary primary", epoch_bjd=t0, period_days=p,
                   epoch_error_days=abs(_number(row.get("e_BJD0")) or 0),
                   period_error_days=abs(_number(row.get("e_Per")) or 0))
    result = []
    for comp, default in (("p",0.0),("s",.5)):
        width = None
        phase = 0.0 if comp=="p" else None
        for method in ("pf","2g"):
            width = width or _number(row.get("W"+comp+"-"+method))
            if phase is None:
                phip = _number(row.get("Phip-"+method))
                phis = _number(row.get("Phis-"+method))
                if phip is not None and phis is not None:
                    phase = (phis-phip)%1
        if phase is None:
            phase = default
            unresolved.append(label+" has no fitted secondary phase; half-phase also masked")
        if not width or width <= 0:
            width = min(.2,.35/p)
            unresolved.append(label+" has no fitted "+comp+" eclipse width")
        result.append(dict(primary, epoch_bjd=t0+p*phase,duration_days=p*width,
                           kind="published binary primary" if comp=="p" else "binary secondary"))
    return result


class Catalogs:
    def __init__(self, cache, snapshot=None):
        self.cache = Path(cache)
        self.sources, self.errors, self.rows = {}, {}, []
        self.memo = {}
        self.vsx = snapshot.get("vsx") if snapshot else None
        if snapshot is None:
            self._load_tables()
        else:
            self.sources = {name:{"sha256":source["sha256"],"source":"catalog-receipts.json"}
                            for name,source in snapshot["sources"].items()}
            self.errors = snapshot["errors"]
            self.rows = snapshot["rows"]
        self.ras = np.array([_number(row.get("ra")) for _,row,_ in self.rows], dtype=float)
        self.decs = np.array([_number(row.get("dec")) for _,row,_ in self.rows], dtype=float)
        self.names = {}
        for i, (_, _, name) in enumerate(self.rows):
            self.names.setdefault(name, []).append(i)
        if self.vsx is not None:
            self.centers = np.asarray(self.vsx["centers"],dtype=float).reshape(-1,3)
            self.source_coords = np.array([[float(r["RAJ2000"]),float(r["DEJ2000"])] for r in self.vsx["rows"]]).reshape(-1,2)

    def _load_tables(self):
        for key, sql in (
                ("pscomppars", "SELECT pl_name,hostname,ra,dec,tran_flag,pl_orbper,pl_orbpererr1,"
                 "pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandur,"
                 "pl_orbeccen,pl_orbeccenerr1,pl_orbeccenerr2,pl_orblper,"
                 "pl_orblpererr1,pl_orblpererr2 FROM pscomppars"),
                ("toi", "SELECT toi,tid,ra,dec,tfopwg_disp,pl_orbper,pl_orbpererr1,"
                 "pl_orbpererr2,pl_tranmid,pl_tranmiderr1,pl_tranmiderr2,pl_trandurh FROM toi")):
            try:
                # Validate the explicit science column list against the live service;
                # ancillary-column differences may be recovered from the full table.
                schema, schema_receipt = _tap(self.cache, NEA,
                    f"SELECT column_name FROM TAP_SCHEMA.columns WHERE table_name='{key}'")
                available = {str(r["column_name"]).lower() for r in schema}
                requested = [c.strip() for c in sql.split(" FROM ")[0][7:].split(",")]
                required = {"ra","dec","pl_orbper","pl_tranmid"}
                required.update({"pl_name","hostname","pl_trandur"} if key=="pscomppars"
                                else {"toi","tid","pl_trandurh"})
                if not required <= available:
                    raise ValueError(f"NASA TAP schema missing {sorted(required-available)}")
                missing = set(requested)-available
                if missing:
                    self.errors[key+"_columns"] = "missing: "+", ".join(sorted(missing))
                cols = [c for c in requested if c in available]
                table, source = _tap(self.cache, NEA, "SELECT "+",".join(cols)+" FROM "+key)
                source["schema_check"] = schema_receipt
                if not table or not {"ra","dec","pl_tranmid"}.issubset(table[0]):
                    raise ValueError("empty or invalid planet table")
                count, check = _tap(self.cache, NEA, f"SELECT COUNT(*) AS total FROM {key}")
                if int(next(iter(count[0].values()))) != len(table):
                    raise ValueError("TAP returned a partial catalog")
                self.sources[key] = dict(source, rows=len(table), count_check=check)
                for row in table:
                    self.rows.append((key, row, _norm(row.get("hostname") or "TIC"+str(row.get("tid")))))
            except Exception as exc:
                self.errors[key] = str(exc)[:500]
        try:
            tabname = '"J/ApJS/258/16/tess-ebs"'
            table, source = _tap(self.cache, VIZIER, "SELECT * FROM "+tabname)
            n, check = _tap(self.cache, VIZIER, "SELECT COUNT(*) AS total FROM "+tabname)
            if int(next(iter(n[0].values()))) != len(table) or not table:
                raise ValueError("partial TESS eclipsing binary catalog")
            self.sources["tess_eb"] = dict(source,rows=len(table),count_check=check)
            for row in table:
                row.update(ra=row.get("RAJ2000") or row.get("RAdeg"),
                           dec=row.get("DEJ2000") or row.get("DEdeg"))
                self.rows.append(("tess_eb",row,_norm("TIC"+str(row["TIC"]))))
        except Exception as exc:
            self.errors["tess_eb"] = str(exc)[:500]

    def _vsx_rows(self, ra, dec):
        if self.vsx is None:
            rows, source = _tap(self.cache, VIZIER,
                'SELECT TOP 501 * FROM "B/vsx/vsx" WHERE '+_cone(ra,dec,125,'"RAJ2000"','"DEJ2000"'))
            if len(rows)>=501:
                raise ValueError("VSX cone exceeded row limit")
            return rows,source
        # A 125 arcsec cone is complete only if it lies entirely inside a
        # successfully retrieved footprint. No per-star mirror requests occur.
        covers = np.flatnonzero(np.abs(self.centers[:,1]-dec)<=(VSX_RADIUS-125)/3600)
        batch = next((int(self.centers[i,2]) for i in covers
                      if _sep(ra,dec,*self.centers[i,:2])<=VSX_RADIUS-125),None)
        if batch is None:
            raise ValueError("no complete 125 arcsec VSX coverage in the shared catalog seed")
        candidates = np.flatnonzero(np.abs(self.source_coords[:,1]-dec)<=125/3600)
        rows = [self.vsx["rows"][i] for i in candidates if _sep(ra,dec,*self.source_coords[i])<=125]
        receipt = self.vsx["batches"][batch]["receipt"]
        return rows,{"batch":batch,"sha256":receipt["sha256"],"source":"catalog-receipts.json"}

    def match(self, star, ra_deg, dec_deg):
        key = (str(star), _number(ra_deg), _number(dec_deg))
        if key in self.memo:
            return self.memo[key]
        name, ra, dec = key
        match = dict(star=name, ra_deg=ra, dec_deg=dec, sources=dict(self.sources),
                     errors=dict(self.errors), ephemerides=[], known_objects=[],
                     unresolved=[], nearby_variables=[])
        close = set()
        if ra is not None and dec is not None:
            # Latitude prefilter is independent of the RA wrap and poles.
            for i in np.flatnonzero(np.abs(self.decs-dec) <= 10/3600):
                if _sep(ra, dec, self.ras[i], self.decs[i]) <= 10:
                    close.add(i)
        normalized_name = _norm(name)
        close.update(self.names.get(normalized_name, []))
        for i in sorted(close):
            table, row, _ = self.rows[i]
            if table=="tess_eb":
                match["known_objects"].append("TESS-EB TIC "+str(row["TIC"]))
                match["ephemerides"].extend(_binary_phases(row,match["unresolved"]))
                continue
            is_toi = table == "toi"
            label = "TOI " + str(row.get("toi")) if is_toi else row["pl_name"]
            match["known_objects"].append(label)
            # False positives remain relevant as potential eclipsing blends.
            if not is_toi and str(row.get("tran_flag")) == "0" and not row.get("pl_tranmid"):
                continue
            match["ephemerides"].extend(_planet_phases(row, is_toi, match["unresolved"]))
        if ra is None or dec is None:
            match["errors"]["vsx"] = "Missing source coordinates"
        else:
            try:
                # VSX reports eclipse type, variability epoch (typically HJD)
                # and period. Preserve original values and pad for HJD/UTC.
                rows, source = self._vsx_rows(ra, dec)
                match["sources"]["vsx"] = source
                for row in rows:
                    sep = _sep(ra, dec, row.get("RAJ2000"), row.get("DEJ2000"))
                    vtype = str(row.get("Type",""))
                    label = str(row.get("Name","VSX source"))
                    if not math.isfinite(sep):
                        raise ValueError("VSX source missing astrometry")
                    match["nearby_variables"].append(dict(name=label, type=vtype, separation_arcsec=sep))
                    if sep > 10 or not re.search(r"(^|[^A-Z])(EA|EB|EW|EC|E|ELL)([^A-Z]|$)", vtype.upper()):
                        continue
                    t0, per = _number(row.get("Epoch")), _number(row.get("Period"))
                    if t0 is not None and 30000 < t0 < 100000:
                        t0 += 2400000.0
                    if not (per and per>0 and t0 and 2400000<t0<2600000):
                        match["unresolved"].append(label+" is a cataloged eclipsing variable with incomplete ephemeris")
                        continue
                    maskdur = min(.2*per,.35)
                    match["unresolved"].append(label+" secondary phase and widths require binary ephemeris review")
                    for phase in (0,.5):
                        match["ephemerides"].append(dict(label=label,kind="VSX binary eclipse",
                            epoch_bjd=t0+phase*per,period_days=per,duration_days=maskdur,
                            epoch_error_days=.01, period_error_days=0.0))
            except Exception as exc:
                match["errors"]["vsx"] = str(exc)[:500]
        if match["unresolved"]:
            match["errors"]["ephemerides"] = "Incomplete published mask; see unresolved entries"
        self.memo[key] = match
        return match


VSX_RADIUS = 300
VSX_COLS = ["RAJ2000","DEJ2000","Name","Type","Period","Epoch"]


def _votable_rows(raw, columns, limit):
    from astropy.table import Table
    xml=ET.fromstring(raw)
    for e in xml.iter():
        if e.tag.rsplit("}",1)[-1]=="INFO" and (e.attrib.get("value","").upper() in ("ERROR","OVERFLOW")
                or e.attrib.get("name","").lower()=="error"):
            raise ValueError("VizieR reported an error or overflow: "+str(e.text)[:150])
    tabs=[e for e in xml.iter() if e.tag.rsplit("}",1)[-1]=="TABLE"]
    if not tabs:
        raise ValueError("VizieR returned no data-table schema")
    table=Table.read(io.BytesIO(raw),format="votable",use_names_over_ids=True)
    if len(table)>=limit or not set(columns)<=set(table.colnames):
        raise ValueError("incomplete VizieR table response: "+str(table.colnames))
    return [{k:None if np.ma.is_masked(row[k]) else str(row[k]) for k in table.colnames}
            for row in table]


def _vsx_fallback(batch):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astroquery.vizier import Vizier
    coord=SkyCoord([p[0] for p in batch],[p[1] for p in batch],unit="deg")
    failures=[]
    for host in ASU_MIRRORS:
        try:
            cat=Vizier(columns=VSX_COLS,row_limit=25000,vizier_server=host)
            cat.TIMEOUT=16
            def query():
                res=cat.query_region_async(coord,radius=VSX_RADIUS*u.arcsec,catalog="B/vsx/vsx",cache=False)
                res.raise_for_status()
                return res
            res=_route(host,query)
            rows=_votable_rows(res.content,VSX_COLS,25000)
            return rows,{"url":res.url,"method":"VizieR.query_region_async(list of ICRS centers)",
                         "sha256":_hash(res.content),"bytes":len(res.content),"mirror_errors":failures}
        except Exception as exc:
            failures.append(f"{host}: {str(exc)[:220]}")
    raise RuntimeError(" | ".join(failures))


def _eb_bulk(cache):
    errors=[]
    for host in ASU_MIRRORS:
        try:
            raw,source=_get(cache,f"https://{host}/viz-bin/votable",{
                "-source":"J/ApJS/258/16/tess-ebs","-out.all":"1","-out.add":"_RAJ2000,_DEJ2000","-out.max":"20000"})
            rows=_votable_rows(raw,["TIC","BJD0","Per"],20000)
            if len(rows)!=4584:
                raise ValueError(f"TESS-EB published sample has 4584 entries; got {len(rows)}")
            for r in rows:
                r["ra"]=r.get("RAdeg") or r.get("_RAJ2000") or r.get("RAJ2000")
                r["dec"]=r.get("DEdeg") or r.get("_DEJ2000") or r.get("DEJ2000")
                if _number(r["ra"]) is None or _number(r["dec"]) is None:
                    raise ValueError("TESS-EB coordinates missing")
            return rows,dict(source,rows=len(rows),mirror_errors=errors)
        except Exception as exc:
            errors.append(f"{host}: {str(exc)[:240]}")
    # Published CDS fixed-width tables can also be retrieved independently
    # of ASU/TAP. The ReadMe supplies the units, column layout and row count.
    from astropy.io import ascii
    for base in ("https://cdsarc.cds.unistra.fr/ftp/J/ApJS/258/16",
                 "https://vizier.cfa.harvard.edu/ftp/cats/J/ApJS/258/16"):
        try:
            meta,mreceipt=_get(cache,base+"/ReadMe",{})
            for suffix in (".gz",""):
                try:
                    byt,source=_get(cache,base+"/tess-ebs.dat"+suffix,{})
                    payload=gzip.decompress(byt) if byt.startswith(b"\x1f\x8b") else byt
                    break
                except Exception as exc:
                    errors.append(str(exc)[:220])
            else:
                continue
            with tempfile.TemporaryDirectory() as d:
                folder=Path(d)
                (folder/"ReadMe").write_bytes(meta)
                (folder/"tess-ebs.dat").write_bytes(payload)
                table=ascii.read(folder/"tess-ebs.dat",format="cds",readme=str(folder/"ReadMe"))
            rows=[{k:None if np.ma.is_masked(r[k]) else str(r[k]) for k in table.colnames} for r in table]
            if len(rows)!=4584:
                raise ValueError("incomplete published EB table")
            for r in rows:
                r.update(ra=r["RAdeg"],dec=r["DEdeg"])
            return rows,dict(source,rows=len(rows),readme=mreceipt,mirror_errors=errors)
        except Exception as exc:
            errors.append(f"{base}: {str(exc)[:240]}")
    raise RuntimeError(" | ".join(errors))


def _vsx_seed(inventory,cache):
    coords = set()
    missing = 0
    for line in Path(inventory).read_text().splitlines():
        row = json.loads(line)
        ra, dec = _number(row.get("ra")),_number(row.get("dec"))
        if (ra is None or dec is None) and row.get("coords_hms_dms"):
            from astropy.coordinates import SkyCoord
            import astropy.units as u
            try:
                raw = str(row["coords_hms_dms"]).strip().replace(","," ")
                parts=raw.split("/")
                if len(parts)==2 and all(_number(v) is not None for v in parts):
                    pos=SkyCoord(float(parts[0]),float(parts[1]),unit="deg")
                else:
                    pos = SkyCoord(raw.replace("/"," "), unit=(u.hourangle, u.deg))
                ra, dec = float(pos.ra.deg), float(pos.dec.deg)
            except (ValueError, TypeError):
                pass
        if ra is None or dec is None:
            missing += 1
        else:
            coords.add((round(ra,5),round(dec,5)))
    # Pointings within 5 arcsec share the same 300 arcsec footprint. Retained
    # centers remain exact. Final coverage is determined geometrically per source.
    centers = []
    for ra,dec in sorted(coords,key=lambda p:p[1]):
        nearby = []
        for prior in reversed(centers):
            if dec-prior[1]>5/3600:
                break
            nearby.append(prior)
        if not any(_sep(ra,dec,*p)<=5 for p in nearby):
            centers.append((ra,dec))
    unique_rows, covered, batches = {}, [], []
    consecutive_failures = 0
    for start in range(0,len(centers),24):
        batch = centers[start:start+24]
        entry = {"centers":batch}
        try:
            if consecutive_failures>=2:
                raise ValueError("continued service outage; batch skipped after two failed batches on both services")
            sql = 'SELECT TOP 25001 '+",".join('"'+c+'"' for c in VSX_COLS)+' FROM "B/vsx/vsx" WHERE '
            sql += " OR ".join(_cone(ra,dec,VSX_RADIUS,'"RAJ2000"','"DEJ2000"') for ra,dec in batch)
            try:
                # Try the US/Japanese ASU services first; Strasbourg TAP
                # timed out on the run-2 runner. A failed multi-route probe
                # falls through to the independent Strasbourg TAP interface.
                rows,receipt=_vsx_fallback(batch)
            except Exception as exc:
                entry["asu_error"]=str(exc)[:1200]
                # Use VOTable to preserve the service's QUERY_STATUS=OVERFLOW
                # signal, which would be lost in plain CSV output.
                from astropy.table import Table
                raw, receipt = _get(cache,VIZIER,{"REQUEST":"doQuery","LANG":"ADQL",
                                      "FORMAT":"votable","MAXREC":25001,"QUERY":sql})
                xml = ET.fromstring(raw)
                if any(e.attrib.get("value") in ("ERROR","OVERFLOW") for e in xml.iter()
                       if e.attrib.get("name")=="QUERY_STATUS"):
                    raise ValueError("VSX TAP reported an error or overflow")
                table = Table.read(io.BytesIO(raw),format="votable",use_names_over_ids=True)
                if len(table)>=25001 or not set(VSX_COLS)<=set(table.colnames):
                    raise ValueError("truncated or invalid VSX TAP response")
                rows = [{k:None if np.ma.is_masked(r[k]) else str(r[k]) for k in VSX_COLS} for r in table]
            for r in rows:
                if _number(r.get("RAJ2000")) is None or _number(r.get("DEJ2000")) is None:
                    raise ValueError("missing VSX coordinates")
                unique_rows[r["Name"],r["RAJ2000"],r["DEJ2000"]] = r
            entry.update(receipt=receipt,row_count=len(rows))
            covered.extend([ra,dec,len(batches)] for ra,dec in batch)
            consecutive_failures = 0
        except Exception as exc:
            entry["error"] = str(exc)[:500]
            consecutive_failures += 1
        batches.append(entry)
    return {"centers":covered,"rows":list(unique_rows.values()),"batches":batches,
            "radius_arcsec":VSX_RADIUS,"visits_missing_coordinates":missing}


def seed_catalogs(inventory,cache_dir,previous_seed=None):
    directory = Path(cache_dir)
    directory.mkdir(parents=True,exist_ok=True)
    if previous_seed:
        with gzip.open(previous_seed,"rt") as src:
            snapshot=json.load(src)
        nasa={"pscomppars","toi"}
        if not nasa.issubset(snapshot["sources"]) or not all(any(r[0]==key for r in snapshot["rows"]) for key in nasa):
            raise ValueError("previous seed lacks verified NASA tables")
        if "tess_eb" not in snapshot["sources"]:
            try:
                rows,source=_eb_bulk(directory)
                snapshot["rows"]=[r for r in snapshot["rows"] if r[0]!="tess_eb"]
                snapshot["rows"].extend(("tess_eb",r,_norm("TIC"+str(r["TIC"]))) for r in rows)
                snapshot["sources"]["tess_eb"]=source
                snapshot["errors"].pop("tess_eb",None)
            except Exception as exc:
                snapshot["errors"]["tess_eb"]=str(exc)[:2400]
    else:
        cats=Catalogs(directory)
        snapshot=dict(sources=cats.sources,errors=cats.errors,rows=cats.rows)
    vsx = _vsx_seed(inventory,directory)
    snapshot["vsx"]=vsx
    encoded = json.dumps(snapshot,separators=(",",":"),allow_nan=False).encode()
    compressed = gzip.compress(encoded,mtime=0)
    (directory/"catalog-seed.json.gz").write_bytes(compressed)
    summary = dict(sources=snapshot["sources"],errors=snapshot["errors"],mask_rows=len(snapshot["rows"]),
        previous_seed_sha256=_filehash(previous_seed) if previous_seed else None,
        inventory_sha256=_filehash(inventory),seed_sha256=_hash(compressed),
        vsx={k:v for k,v in vsx.items() if k not in ("centers","rows")},
        vsx_stars=len(vsx["rows"]),vsx_successful_pointings=len(vsx["centers"]))
    (directory/"catalog-receipts.json").write_text(json.dumps(summary,indent=2)+"\n")


def prepare_catalogs(cache_dir):
    seed = os.environ.get("CHEOPS_CATALOG_SEED")
    if seed is not None:
        with gzip.open(seed,"rt") as source:
            return Catalogs(cache_dir,json.load(source))
    return Catalogs(cache_dir)


def _gaia(event, cache):
    ra, dec = event["ra_deg"], event["dec_deg"]
    rows, source = _tap(cache, GAIA,
        "SELECT TOP 3000 source_id,ra,dec,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE " +
        _cone(ra, dec, 125, "ra", "dec"))
    ranked = sorted([(_sep(ra,dec,r["ra"],r["dec"]), r) for r in rows], key=lambda x:x[0])
    result = {"source":source, "cone_arcsec":125, "rows":len(rows)}
    if not ranked or ranked[0][0]>3 or any(_number(r.get("phot_g_mean_mag")) is None for _,r in ranked):
        return dict(result, status="incomplete", reason="ambiguous target or missing Gaia magnitudes")
    if len(ranked)>=3000 or (len(ranked)>1 and ranked[1][0]<3):
        return dict(result, status="incomplete", reason="ambiguous Gaia association or truncated cone")
    target=ranked[0][1]
    fluxratios=[10**(-.4*(float(r["phot_g_mean_mag"])-float(target["phot_g_mean_mag"]))) for _,r in ranked[1:]]
    maxflux=sum(fluxratios)
    neighbors=[dict(source_id=r["source_id"],sep_arcsec=sep,max_fraction=ratio/(1+ratio))
               for (sep,r),ratio in zip(ranked[1:], fluxratios) if ratio/(1+ratio)>=event["depth"]*.5]
    # Gaia G is a proxy for the CHEOPS band. This deliberately cannot replace
    # a pixel/difference-image localization, including unresolved companions.
    return dict(result,status="screened", target_source_id=target["source_id"],
                target_separation_arcsec=ranked[0][0],bright_neighbors=neighbors,
                summed_catalog_flux_fraction=maxflux/(1+maxflux))


def _simbad(cache, query):
    errors=[]
    for url in SIMBAD_MIRRORS:
        try:
            return _tap(cache,url,query)
        except Exception as exc:
            errors.append(f"{url}: {str(exc)[:180]}")
    raise RuntimeError(" | ".join(errors))


def _literature(event, cache):
    ra, dec = event["ra_deg"],event["dec_deg"]
    cone = _cone(ra,dec,3,"b.ra","b.dec")
    sources,errors=[],{}
    def lookup(name,query):
        try:
            rows,source=_simbad(cache,query)
            sources.append(source)
            return rows
        except Exception as exc:
            errors[name]=str(exc)[:500]
            return []
    ids=lookup("SIMBAD identifiers","SELECT TOP 501 b.main_id,i.id FROM basic b "
        "JOIN ident i ON i.oidref=b.oid WHERE "+cone)
    refs=lookup("SIMBAD bibliography","SELECT DISTINCT TOP 1001 r.bibcode,r.title FROM basic b "
        "JOIN has_ref h ON h.oidref=b.oid JOIN ref r ON r.oidbib=h.oidbibref WHERE "+cone)
    aliases = {str(event["star"])}
    aliases.update(str(r["main_id"]) for r in ids)
    for r in ids:
        alias = str(r["id"]).strip()
        if re.match(r"^(HD |HIP |TOI[- ]|GJ |WASP[- ]|KELT[- ]|HAT-|TIC |LHS |LTT |55 Cnc|K2-)",alias):
            aliases.add(alias)
    query_aliases=sorted(aliases, key=lambda n:(n!=str(event["star"]),len(n),n))[:10]
    match = " OR ".join('all:"'+re.sub(r'["\r\n]'," ",a)+'"' for a in query_aliases)
    query = "("+match+") AND (all:transit OR all:eclipse OR all:CHEOPS OR all:occultation)"
    total,papers=0,[]
    try:
        xml, source_arxiv = _get(cache,ARXIV,dict(search_query=query,start=0,max_results=100))
        root=ET.fromstring(xml)
        ns={"a":"http://www.w3.org/2005/Atom","os":"http://a9.com/-/spec/opensearch/1.1/"}
        if root.tag != "{"+ns["a"]+"}feed":
            raise ValueError("arXiv did not return an Atom feed")
        papers=[dict(id=e.findtext("a:id",namespaces=ns),title=e.findtext("a:title",namespaces=ns),
                     abstract=e.findtext("a:summary",namespaces=ns)) for e in root.findall("a:entry",ns)]
        if any("/api/errors" in str(p["id"]) for p in papers):
            raise ValueError("arXiv query error")
        total=int(root.findtext("os:totalResults","0",ns))
        sources.append(source_arxiv)
    except Exception as exc:
        errors["arXiv"]=str(exc)[:400]
    if errors:
        try:
            raw,source=_get(cache,"https://api.crossref.org/works",{"query.bibliographic":
                " ".join(query_aliases[:3])+" transit eclipse CHEOPS","rows":20,"select":"DOI,title"})
            items=json.loads(raw)["message"]["items"]
            sources.append(source)
            refs=refs+[{"doi":v.get("DOI"),"title":" ".join(v.get("title") or []),
                         "match":"Crossref keyword result; verify source and coordinates"} for v in items]
        except Exception as exc:
            errors["Crossref"]=str(exc)[:250]
    return dict(sources=sources,aliases=query_aliases,associated_references=refs,
                arxiv_records=papers,arxiv_total=total,search_errors=errors,
                truncated=len(refs)>=1001 or len(ids)==501 or total>100,
                status="incomplete" if errors else "requires publication and event epoch review")


def _sector_period_search(time, flux, event):
    # Twenty-minute medians preserve resolved eclipses while bounding the cost
    # of the secondary archive check on short-cadence products.
    from astropy.timeseries import BoxLeastSquares
    order = np.argsort(time)
    t, f = np.asarray(time)[order], np.asarray(flux, dtype=float)[order]
    usable = np.isfinite(t) & np.isfinite(f)
    t, f = t[usable], f[usable]
    if len(t) < 60 or np.median(f) <= 0 or np.ptp(t) < 1.5:
        return {"status": "insufficient baseline"}
    f = f / np.median(f)
    index = np.floor((t-t[0])*72).astype(int)
    pieces = np.split(np.arange(len(t)), np.flatnonzero(np.diff(index))+1)
    tb = np.array([np.median(t[k]) for k in pieces])
    fb = np.array([np.median(f[k]) for k in pieces])
    # Remove changes on scales longer than the searched eclipse durations.
    from scipy.ndimage import median_filter
    chunks = np.split(np.arange(len(tb)), np.flatnonzero(np.diff(tb)>.15)+1)
    for chunk in chunks:
        fb[chunk] -= median_filter(fb[chunk], size=min(73,len(chunk)), mode="nearest")
    sigma = 1.4826*np.median(np.abs(fb-np.median(fb)))
    if not sigma > 0:
        return {"status": "no measurable variability"}
    bls = BoxLeastSquares(tb,fb,dy=np.full(len(tb),sigma))
    dgrid = np.array([1,2,3,5])/24
    result = bls.autopower(dgrid,minimum_period=.45,maximum_period=min(15,np.ptp(tb)/2),
                           frequency_factor=3,minimum_n_transit=2)
    best = int(np.nanargmax(result.power))
    per, t0, dur = (float(v[best]) for v in (result.period,result.transit_time,result.duration))
    def relative_phase(times):
        return (times-t0+per/2)%per-per/2
    in_eclipse = np.abs(relative_phase(tb)) < dur/2
    transits = np.rint((tb[in_eclipse]-t0)/per).astype(int)
    cycles = int(len(np.unique(transits)))
    depth,snr = float(result.depth[best]),float(result.depth_snr[best])
    return {"status":"periodic dip found" if cycles>=2 and snr>=9 else "no strong periodic dip",
            "epoch_bjd":t0,"period_days":per,"duration_days":dur,"depth_fraction":depth,
            "white_noise_snr":snr,"sampled_cycles":cycles,
            "cheops_phase_offset_days":float(relative_phase(event["epoch_bjd"])),
            "minimum_period_days":.45,"maximum_period_days":min(15,np.ptp(tb)/2)}


def _difference_image(t, cube, x, y, epoch, dur, period=None):
    # For a single event use the actual time. For repeating TESS events use
    # local in/out samples from the same cycles, reducing pointing drift.
    phase = t-epoch
    if period:
        phase = (phase+period/2)%period-period/2
    inside = np.abs(phase)<dur/2
    outside = (np.abs(phase)>.7*dur)&(np.abs(phase)<min(2*dur,period/2 if period else 2*dur))
    if inside.sum()<3 or outside.sum()<8:
        return {"status":"insufficient in/out pixels"}
    diff = np.nanmedian(cube[outside],axis=0)-np.nanmedian(cube[inside],axis=0)
    finite = np.isfinite(diff)
    if finite.sum()<10:
        return {"status":"insufficient valid pixels"}
    sky = np.nanmedian(diff)
    noise = 1.4826*np.nanmedian(np.abs(diff-sky))
    signal = np.where(finite,np.maximum(diff-sky-3*noise,0),0)
    if signal.sum()<=0 or np.max(signal)<=2*noise:
        return {"status":"unresolved difference image"}
    yy,xx=np.indices(signal.shape)
    dx = float(np.sum(signal*(xx-x))/signal.sum())
    dy = float(np.sum(signal*(yy-y))/signal.sum())
    return {"status":"measured", "offset_pixels":float(math.hypot(dx,dy)),
            "delta_column":dx, "delta_row":dy, "peak_above_sky":float(np.nanmax(diff)-sky),
            "spatial_mad":float(noise),"in_cadences":int(inside.sum()),"out_cadences":int(outside.sum())}


SECTOR_SCREENS = {}


def _compare_flux(time, flux, event, product_id):
    duration, epoch=event["duration_days"],event["epoch_bjd"]
    good=np.isfinite(time)&np.isfinite(flux)
    t, f=np.asarray(time)[good], np.asarray(flux)[good]
    if not len(t):
        return dict(status="no valid cadences")
    dt=np.abs(t-epoch)
    inside, outside=dt<duration/2, (dt>duration*.7)&(dt<duration*2)
    if product_id not in SECTOR_SCREENS:
        if len(SECTOR_SCREENS)>512:
            SECTOR_SCREENS.clear()
        SECTOR_SCREENS[product_id]=_sector_period_search(t,f,event)
    screen=dict(SECTOR_SCREENS[product_id])
    if "period_days" in screen:
        per=screen["period_days"]
        screen["cheops_phase_offset_days"]=float((epoch-screen["epoch_bjd"]+per/2)%per-per/2)
    out=dict(start_bjd=float(t.min()),end_bjd=float(t.max()),n=len(t),in_event=int(inside.sum()),
             periodic_screen=screen)
    if inside.sum()<3 or outside.sum()<8 or not np.any(t[outside]<epoch) or not np.any(t[outside]>epoch):
        return dict(out,status="no usable simultaneous baseline")
    base=float(np.median(f[outside]))
    if not base>0:
        return dict(out,status="nonpositive local flux")
    norm=f/base
    depth=1-float(np.median(norm[inside]))
    med=np.median(norm[outside])
    noise=1.4826*np.median(np.abs(norm[outside]-med))
    return dict(out,status="simultaneous comparison",depth=depth,
                out_of_event_scatter=float(noise),in_event_median_uncertainty_white=float(1.253*noise/math.sqrt(inside.sum())))


def _sector_order(event,sectors,cache):
    # Mission Guest Investigator office schedule, checked against product BJD
    # ranges whenever comparing an actual sector.
    url="https://raw.githubusercontent.com/tessgi/tesswcs/main/src/tesswcs/data/pointings.csv"
    raw,receipt=_get(cache,url,{})
    rows=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    ranges={int(r["Sector"]):(float(r["Start"]),float(r["End"])) for r in rows if r.get("Sector")}
    missing=sorted(set(sectors)-ranges.keys())
    if missing:
        # Supplement newly observed sectors using the mission's orbit schedule.
        from astropy.time import Time
        extra,extra_src=_get(cache,"https://tess.mit.edu/public/files/TESS_orbit_times.csv",{})
        times={}
        for r in csv.DictReader(l for l in extra.decode("utf-8-sig").splitlines() if not l.startswith("#")):
            sec=_number(r.get("Sector"))
            if sec not in missing:
                continue
            pair=[float(Time(r[k].strip().replace("T"," "),scale="utc").tdb.jd)
                  for k in ("Start of Orbit","End of Orbit")]
            times.setdefault(int(sec),[]).extend(pair)
        ranges.update({sec:(min(pair),max(pair)) for sec,pair in times.items()})
        receipt["supplement"]=extra_src
    epoch=event["epoch_bjd"]
    unknown=sorted(set(sectors)-ranges.keys())
    order=sorted((s for s in sectors if s in ranges),key=lambda s:max(ranges[s][0]-epoch,epoch-ranges[s][1],0))
    return order+unknown,receipt,unknown


def _tess(event, cache):
    from astropy import units as u
    from astropy.coordinates import SkyCoord
    from astropy.io import fits
    from astroquery.mast import Catalogs as MastCatalogs, Observations, Tesscut

    pos=SkyCoord(event["ra_deg"],event["dec_deg"],unit="deg")
    for client in (MastCatalogs,Observations,Tesscut):
        client.TIMEOUT=35
    sectors=_route("MAST sector footprints",lambda:Tesscut.get_sectors(coordinates=pos))
    sector_ids=sorted(set(int(s) for s in sectors["sector"])) if len(sectors) else []
    out=dict(sectors_on_silicon=sector_ids,sectors_examined=[],products=[],errors=[])
    if not sector_ids:
        return dict(out,status="no TESS sector footprint")
    order,receipt,unknown=_sector_order(event,sector_ids,cache)
    out.update(sector_schedule=receipt,undated_sectors=unknown,
               sector_policy="Nearest usable sector by date, overlapping sector first")
    tesspath=Path(cache)/"tess"
    tesspath.mkdir(parents=True,exist_ok=True)
    try:
        stars=_route("MAST TIC match",lambda:MastCatalogs.query_region(pos,radius=3*u.arcsec,catalog="TIC"))
        tic=str(stars[0]["ID"]) if len(stars)==1 else None
    except Exception as exc:
        tic=None
        out["tic_lookup_error"]=str(exc)[:350]
    if tic:
        out["tic_id"]=tic
    for sec in order:
        out["sectors_examined"].append(sec)
        sel=[]
        if tic:
            try:
                sel=_route("MAST target light curves",lambda:Observations.query_criteria(
                      target_name=tic,obs_collection="TESS",dataproduct_type="timeseries",sequence_number=sec))
            except Exception as exc:
                out["errors"].append(f"SPOC sector {sec}: {str(exc)[:200]}")
        if len(sel):
            try:
                plist=_route("MAST product lists",lambda:Observations.get_product_list(sel))
            except Exception as exc:
                plist=[]
                out["errors"].append(f"SPOC product list, sector {sec}: {str(exc)[:160]}")
            for r in plist:
                fname=str(r["productFilename"])
                if not fname.endswith("-s_lc.fits") and not fname.endswith("-a_fast-lc.fits"):
                    continue
                uri=str(r["dataURI"])
                path=tesspath/Path(fname).name
                try:
                    if not path.exists():
                        status,msg,_=_route("MAST data file",lambda:Observations.download_file(uri,local_path=str(path)))
                        if status!="COMPLETE":
                            raise ValueError(str(msg))
                    with fits.open(path) as h:
                        hdr, tab=h[1].header,h[1].data
                        if str(h[0].header.get("TICID"))!=tic:
                            raise ValueError("TIC identifier mismatch")
                        ok=np.asarray(tab["QUALITY"])==0
                        t=np.asarray(tab["TIME"][ok])+hdr["BJDREFI"]+hdr.get("BJDREFF",0.)
                        ans=_compare_flux(t,tab["PDCSAP_FLUX"][ok],event,uri)
                        sector=int(h[0].header["SECTOR"])
                        if sector!=sec:
                            continue
                        out["products"].append(dict(uri=uri,sha256=_filehash(path),sector=sector,
                                                    dilution=hdr.get("CROWDSAP"),**ans))
                        if ans.get("n",0)>=60:
                            return dict(out,status="incomplete" if unknown else "compared")
                except Exception as exc:
                    out["errors"].append(f"{uri}: {exc}"[:400])
        # Share one local cutout per star/sector across multiple dip epochs.
        try:
            from astropy.table import Table
            folder=tesspath/f"ffi-{pos.ra.deg:.7f}-{pos.dec.deg:.7f}-{sec}"
            folder.mkdir(exist_ok=True)
            files=list(folder.glob("*.fits"))
            manifest=(Table({"Local Path":[str(p) for p in files]}) if files else
                      _route("MAST FFI download",lambda:Tesscut.download_cutouts(
                      coordinates=pos,size=9,sector=sec,path=str(folder))))
            if not len(manifest):
                raise ValueError("empty cutout response")
            for row in manifest:
                path=Path(str(row["Local Path"]))
                with fits.open(path) as h:
                    tab, hdr=h[1].data,h[1].header
                    quality=np.asarray(tab["QUALITY"])==0
                    cube=np.asarray(tab["FLUX"][quality],dtype=float)
                    if not len(cube):
                        raise ValueError("zero quality-filtered FFI cadences")
                    # Use the target's actual WCS location, not an assumed pixel.
                    from astropy.wcs import WCS
                    x,y=WCS(h[2].header).celestial.world_to_pixel(pos)
                    yy,xx=np.indices(cube.shape[1:])
                    r2=(xx-x)**2+(yy-y)**2
                    aperture, sky=r2<=1.6**2,r2>=3**2
                    if aperture.sum()<4 or sky.sum()<8:
                        raise ValueError("insufficient aperture/sky pixels")
                    f=np.sum(cube[:,aperture],axis=1)-aperture.sum()*np.nanmedian(cube[:,sky],axis=1)
                    t=np.asarray(tab["TIME"][quality])+hdr["BJDREFI"]+hdr.get("BJDREFF",0.)
                    ans=_compare_flux(t,f,event,str(path))
                    ans["cheops_epoch_difference_image"]=_difference_image(t,cube,x,y,event["epoch_bjd"],event["duration_days"])
                    repeated=ans.get("periodic_screen",{})
                    if repeated.get("status")=="periodic dip found":
                        ans["periodic_difference_image"]=_difference_image(t,cube,x,y,repeated["epoch_bjd"],
                                repeated["duration_days"],repeated["period_days"])
                    out["products"].append(dict(file=path.name,source="MAST TESSCut",sha256=_filehash(path),
                                                sector=sec,**ans))
                    if ans.get("n",0)>=60:
                        return dict(out,status="incomplete" if unknown else "compared")
        except Exception as exc:
            out["errors"].append(f"FFI sector {sec}: {exc}"[:400])
    return dict(out,status="incomplete")


def known_overlaps(event,match):
    overlaps=[]
    for ep in match["ephemerides"]:
        p=ep["period_days"]
        n=round((event["epoch_bjd"]-ep["epoch_bjd"])/p) if p else 0
        span=.5*(event["duration_days"]+ep["duration_days"])+3*(ep["epoch_error_days"]+abs(n)*ep["period_error_days"])
        if abs(event["epoch_bjd"]-ep["epoch_bjd"]-n*(p or 0))<=span:
            overlaps.append(ep["label"]+" "+ep["kind"])
    return overlaps


LITERATURE_CACHE={}


def vet_event(event,match,cache_dir,previous=None):
    event={**{k:match.get(k) for k in ("ra_deg","dec_deg","star")},**event}
    checks={}
    old=(previous or {}).get("checks",{})
    for key,action in (("field",_gaia),("tess",_tess),("literature",_literature)):
        frozen=old.get(key,{})
        if frozen.get("status") in ("screened","compared","no TESS sector footprint",
                                    "requires publication and event epoch review"):
            checks[key]=frozen
            continue
        try:
            if key=="literature":
                cache_key=(event["star"],event["ra_deg"],event["dec_deg"])
                if cache_key not in LITERATURE_CACHE:
                    if len(LITERATURE_CACHE)>96:
                        LITERATURE_CACHE.clear()
                    LITERATURE_CACHE[cache_key]=action(event,cache_dir)
                # The caller moves references out of individual event records.
                checks[key]=copy.deepcopy(LITERATURE_CACHE[cache_key])
            else:
                checks[key]=action(event,cache_dir)
        except Exception as exc:
            checks[key]={"status":"incomplete","error":str(exc)[:800]}
    overlaps=known_overlaps(event,match)
    return dict(disposition="catalog ephemeris overlap" if overlaps else "requires review",
                known_ephemeris_overlaps=overlaps, known_objects=match["known_objects"],
                catalog_errors=match["errors"],unresolved=match["unresolved"],
                nearby_variables=match["nearby_variables"],sources=match["sources"],checks=checks,
                required_review=["CHEOPS source localization and reduction agreement",
                                 "TESS coverage and dilution","alias and published epoch reconciliation"])


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation",choices=("seed","repair"))
    parser.add_argument("inventory",type=Path)
    parser.add_argument("cache_dir",type=Path)
    parser.add_argument("--previous-seed",type=Path)
    args=parser.parse_args()
    if (args.operation=="repair") != bool(args.previous_seed):
        parser.error("repair requires --previous-seed; seed does not take it")
    seed_catalogs(args.inventory,args.cache_dir,args.previous_seed)
