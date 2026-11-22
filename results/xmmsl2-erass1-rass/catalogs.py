"""Acquire public catalogues, compare individual slew detections in a common soft band."""
import argparse
import gzip
import hashlib
import itertools
import json
import re
import shutil
import tarfile
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.spatial import cKDTree

NASA = "https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
CDS = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
MPE = "https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/"
BASE = "https://heasarc.gsfc.nasa.gov/FTP/"
F = {
    "ra": "ra raj2000 raicrs radeg ra2000",
    "dec": "dec dej2000 deicrs decdeg dec2000",
}
SCHEMA = {
    "xmm": dict(F, name="name srcname xmmsl2 sourcename iauname designation uniquesrcname",
                obs="obsid observationid observationnumber revolution",
                time="time starttime observationdate dateobs mjdstart mjdobs mjd",
                id="sourcenum detid detectionid uniquedetectionid",
                rate="rateb6", counts="sctsb6", exp="expmapb6",
                pos="radecerr poserr positionerror errorradius",
                flux="fluxb6 fluxband6 b6flux softflux fluxsoft fluxsoftband fluxsb flux6",
                error="fluxb6err fluxb6error fluxband6error b6fluxerr softfluxerror softfluxerr fluxsofterror fluxsofterr fluxsberr efluxb6 eflux6 flux6err",
                ml="detmlb6 detmlband6 b6detml softlikelihood softbandlikelihood softdetlike detmlsoft lhb6",
                ext="extb6 extentb6 b6ext softextent extentsoft",
                issue="verpsusp verpsuspflag suspectflag qualityflag",
                val="valflag validationflag", xtyp="idcategory type category",
                counterpart="ident idident identification idname", oldrass="rassname", iddist="iddist", altident="altident", unique="uniquesrcname",
                veropt="veroptload", verext="verinext", verhalo="verhalo",
                verbg="verhibgnd verhighbgnd", verfalse="verfalse", veredge="vernredg"),
    "ero": dict(F, name="name iauname erass1 srcname sourcename designation",
                pos="radecerr poserr positionerror errorradius",
                flux="mlflux1 mlflux0 flux0 broadflux fluxb flux broadbandflux",
                error="mlfluxerr1 mlfluxerr0 emlflux0 fluxerr0 eflux0 broadfluxerror broadfluxerr fluxerror fluxerr eflux",
                ml="detlike0 detectionlikelihood broadbandlikelihood detlike likelihood",
                ext="ext extent extent0 sourceextent",
                extlike="extlike extentlikelihood"),
    "rass": dict(F, name="name iauname 2rxs irocess rosat sourcename",
                 pos="xerr poserrerr poserr errorradius positionerror",
                 posy="yerr",
                 rate="crates count_rate countrate crate srccountrate rate ct",
                 errrate="crateserr countrateerror countrateerr countrateerr1 ecrate ecrates crateerr erate rateerr",
                 ml="eximl existence_likelihood existlikelihood detectionlikelihood likelihood",
                 flag="sourcequalityflag sflag screeningflag spuriousflag suspectflag",
                 near="nearbysrcdetflag", region="extendedregionflag"),
}


def norm(value):
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def qident(name):
    return '"' + name.strip('"').replace('"', '""') + '"'


class Archive:
    def __init__(self, work, out):
        self.downloads = Path(work) / "downloads"
        self.work = Path(work) / "catalog-cache"
        self.work.mkdir(parents=True, exist_ok=True)
        self.out = Path(out)
        self.out.mkdir(parents=True, exist_ok=True)
        self.records = {}
        self.failures = {}
        self.session = requests.Session()
        self.session.headers["User-Agent"] = "Open-Space-XMMSL2-public-archive-analysis/1.0"

    def fetch(self, url, tag, params=None):
        dest = self.work / tag
        frozen = {
            BASE + "xmm/data/catalogues/xmmsl2_clean.fits.gz": "xmm_clean.fits.gz",
            BASE + "xmm/data/catalogues/xmmsl2_total.fits.gz": "xmm_full.fits.gz",
            MPE + "MerloniA_DR1/eRASS1_Main.v1.2.fits.tar.gz": "ero_main.tar.gz",
            MPE + "MerloniA_DR1/eRASS1_Supp.tar.gz": "ero_supp.tar.gz",
            BASE + "heasarc/dbase/msc/heasarc_rass2rxs.fits.gz": "rass.fits.gz",
        }
        if url in frozen and (self.downloads / frozen[url]).is_file() and not dest.exists():
            dest.hardlink_to(self.downloads / frozen[url])
        if not dest.exists():
            with self.session.get(url, params=params, stream=True, timeout=(12, 125)) as r:
                r.raise_for_status()
                with dest.with_suffix(".part").open("wb") as stream:
                    for block in r.iter_content(2**20):
                        stream.write(block)
            dest.with_suffix(".part").replace(dest)
        if tag not in self.records:
            with dest.open("rb") as stream:
                digest = hashlib.file_digest(stream, "sha256").hexdigest()
            self.records[tag] = {"url": url, "params": params, "sha256": digest,
                                 "bytes": dest.stat().st_size}
            self.save()
        return dest

    def tap(self, url, sql, maxrec=2000000):
        tag = hashlib.sha256((url + sql + str(maxrec)).encode()).hexdigest()[:20] + ".csv"
        params = {"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": ("votable" if url == NASA else "csv"),
                  "QUERY": sql, "MAXREC": maxrec}
        path = self.fetch(url, tag, params)
        with path.open("rb") as stream:
            preview = stream.read(512)
        if preview.lstrip().startswith(b"<"):
            from astropy.io.votable import parse
            vot = parse(str(path))
            for info in list(vot.infos) + [i for r in vot.resources for i in r.infos]:
                if str(info.name).upper() == "QUERY_STATUS" and str(info.value).upper() in {"ERROR", "OVERFLOW"}:
                    raise ValueError(str(info.content or info.value))
            return vot.get_first_table().to_table(use_names_over_ids=True).to_pandas()
        return pd.read_csv(path, low_memory=False, skipinitialspace=True)

    def save(self):
        (self.out / "catalog_inputs.json").write_text(
            json.dumps({"resources": self.records, "unsuccessful_routes": self.failures}, indent=2))

    def note(self, key, exc):
        self.failures[key] = str(exc)[:650]
        self.save()
        print("route:", key, str(exc)[:200], flush=True)


def select_columns(role, names, descriptions):
    mapping = {}
    lookup = {norm(n): n for n in names}
    for key, aliases in SCHEMA[role].items():
        for alias in aliases.split():
            if norm(alias) in lookup:
                mapping[key] = lookup[norm(alias)]
                break
    # HEASARC's normalized names and the release names differ between services.
    for key in ["flux", "error", "ml", "ext"]:
        if key in mapping or role != "xmm":
            continue
        for name in names:
            desc = str(descriptions.get(name, "")).lower()
            n = norm(name)
            soft = "soft" in n or "b6" in n or "band6" in n or "soft band" in desc or "soft-band" in desc
            iserr = "err" in n or "error" in desc or "uncertainty" in desc
            if not soft:
                continue
            if key in {"flux", "error"} and "flux" in n and (key == "error") == iserr:
                mapping[key] = name
            if key == "ml" and ("likelihood" in n or "detml" in n):
                mapping[key] = name
            if key == "ext" and ("extent" in n or "extb6" in n) and not iserr and "lik" not in n:
                mapping[key] = name
    required = {"ra", "dec", "name", "flux", "error"} if role != "rass" else {
        "ra", "dec", "name", "rate", "errrate"}
    if required - mapping.keys():
        raise ValueError(f"missing {required-mapping.keys()}; available {list(names)}")
    return mapping


def get_tap_table(arc, endpoint, table, role):
    # VizieR TAP stores table_name with surrounding quotes on some mirrors.
    path = table.strip('"')
    clean = path.replace("'", "''")
    metadata = arc.tap(endpoint, "SELECT column_name,unit,description FROM TAP_SCHEMA.columns "
                       f"WHERE table_name='{clean}' OR table_name='\"{clean}\"'")
    metadata.columns = metadata.columns.str.lower()
    metadata["column_name"] = metadata.column_name.str.strip('"')
    desc = dict(zip(metadata.column_name, metadata.description))
    units = dict(zip(metadata.column_name, metadata.unit))
    mapping = select_columns(role, list(metadata.column_name), desc)
    fields = list(dict.fromkeys(mapping.values()))
    sqlfrom = " FROM " + qident(path)
    # Some NASA tables are advertised with a schema qualifier.
    if "." in path and "/" not in path:
        sqlfrom = " FROM " + ".".join(qident(s) for s in path.split("."))
    count = int(arc.tap(endpoint, "SELECT COUNT(*) AS n" + sqlfrom).iloc[0, 0])
    if count <= 0:
        raise ValueError("empty catalogue")
    ra = qident(mapping["ra"])
    divisions = max(1, (count + 99999) // 100000)
    chunks = []
    for i in range(divisions):
        where = f" WHERE {ra}>={360*i/divisions:.12f} AND {ra}<{360*(i+1)/divisions:.12f}"
        sql = "SELECT " + ",".join(map(qident, fields)) + sqlfrom + where
        chunks.append(arc.tap(endpoint, sql))
    frame = pd.concat(chunks, ignore_index=True)
    if len(frame) != count:
        raise ValueError(f"catalogue truncation: expected {count}, got {len(frame)}")
    return frame.rename(columns={v: k for k, v in mapping.items()}), mapping, units, desc


def get_tdat(arc, table, role):
    # HEASARC's full database exports avoid TAP result-size limits.
    url = BASE + f"heasarc/dbase/tdat_files/heasarc_{table}.tdat.gz"
    path = arc.fetch(url, table + ".tdat.gz")
    header = []
    with gzip.open(path, "rt", errors="replace") as stream:
        for line in stream:
            if line.strip().upper() == "<DATA>":
                break
            header.append(line.rstrip())
            if len(header) > 20000:
                raise ValueError("invalid TDAT header")
        text = "\n".join(header)
        (arc.out / f"schema_{table}.txt").write_text(text)
        layout = re.search(r"(?im)^\s*line\[1\]\s*=\s*(.*(?:\n[ \t]+[^=\n]+)*)", text)
        names = re.findall(r"[a-zA-Z][\w]*", layout[1]) if layout else []
        fields = re.findall(r"(?im)^\s*field\[([^\]]+)\]\s*=\s*(.+)$", text)
        if not names:
            names = [k if not k.isdigit() else v.split()[0] for k, v in fields]
        desc = {k: v for k, v in fields}
        mapping = select_columns(role, names, desc)
        idx = [names.index(v) for v in mapping.values()]
        rows = []
        for line in stream:
            if line.lstrip().startswith("<"):
                break
            vals = line.rstrip("\r\n").split("|")
            if len(vals) < max(idx) + 1:
                raise ValueError("short TDAT record, invalid layout")
            rows.append([vals[i].strip() for i in idx])
    if not rows:
        raise ValueError("no TDAT rows")
    return pd.DataFrame(rows, columns=list(mapping)), mapping, {}, desc


def get_fits(arc, url, tag, role):
    path = arc.fetch(url, tag)
    local = arc.work / (tag + ".unpacked")
    with path.open("rb") as f:
        compressed = f.read(2) == b"\x1f\x8b"
    if tag.endswith((".tar.gz", ".tgz")):
        if not local.exists():
            with tarfile.open(path) as tar:
                group = [m for m in tar.getmembers() if m.isfile() and m.name.lower().endswith((".fits", ".fit"))]
                if not group:
                    raise ValueError("no FITS file in archive")
                member = max(group, key=lambda m: m.size)
                with tar.extractfile(member) as a, local.open("wb") as b:
                    shutil.copyfileobj(a, b)
        compressed = True
    elif compressed and not local.exists():
        with gzip.open(path, "rb") as a, local.open("wb") as b:
            shutil.copyfileobj(a, b)
    with fits.open(local if compressed else path, memmap=True) as hdus:
        h = next(h for h in hdus if isinstance(h, fits.BinTableHDU))
        units = {c.name: c.unit for c in h.columns}
        desc = {c.name: h.header.comments[f"TTYPE{i+1}"] for i, c in enumerate(h.columns)}
        mapping = select_columns(role, list(units), desc)
        if role == "ero":
            mapping.update({n.lower(): n for n in units if n.upper().startswith("FLAG_SP_")})
        frame = {}
        for key, name in mapping.items():
            arr = np.array(h.data[name])
            if arr.dtype.kind == "S":
                arr = np.char.decode(arr, errors="replace")
            if arr.dtype.byteorder not in {"=", "|"} and not arr.dtype.isnative:
                arr = arr.byteswap().view(arr.dtype.newbyteorder("="))
            frame[key] = arr
        return pd.DataFrame(frame), mapping, units, desc


def discover_tables(arc, endpoint, predicate):
    meta = arc.tap(endpoint, "SELECT table_name,description FROM TAP_SCHEMA.tables WHERE " + predicate)
    names = meta.iloc[:, 0].astype(str).str.strip('"')
    return list(zip(names, meta.iloc[:, 1].fillna("").astype(str)))


def get_catalog(arc, key):
    role = "xmm" if key.startswith("xmmslew") else ("rass" if key == "rass" else "ero")
    options = []
    if role == "xmm":
        # The rolling HEASARC xmmslewcln/xmmslewful are XMMSL3.1 as of June 2026.
        # Only explicit, version-frozen XMMSL2 files qualify here.
        file = "xmmsl2_clean.fits.gz" if key.endswith("cln") else "xmmsl2_total.fits.gz"
        options.append((get_fits, (BASE + "xmm/data/catalogues/" + file, file, role)))
    elif role == "ero":
        file = "eRASS1_Supp.tar.gz" if key == "ero_supp" else "eRASS1_Main.v1.2.fits.tar.gz"
        options.append((get_fits, (MPE + "MerloniA_DR1/" + file, file, role)))
    else:
        options.extend([(get_fits, (BASE + "heasarc/dbase/msc/heasarc_rass2rxs.fits.gz", "rass.fits.gz", role)),
                        (get_tap_table, (NASA, "rass2rxs", role)), (get_tdat, ("rass2rxs", role))])
    # Probe additional archives only as needed.
    def try_options(opts):
        for fn, args in opts:
            try:
                frame, mapping, units, desc = fn(arc, *args)
                if len(frame) == 0:
                    raise ValueError("empty source table")
                info = {"rows": len(frame), "mode": fn.__name__, "location": args[:-1],
                        "columns": {k: {"field": n, "unit": str(units.get(n)),
                                          "description": str(desc.get(n))} for k, n in mapping.items()}}
                print(key, info["rows"], mapping, flush=True)
                return frame, info
            except Exception as exc:
                arc.note(key + " " + str(args[:-1]), exc)
        return None
    result = try_options(options)
    if result or role == "xmm":
        return result
    for endpoint in (CDS, NASA):
        cond = ("table_name LIKE '%682/A34%'" if role == "ero" else "table_name LIKE '%588/A103%'") if endpoint == CDS else (
            "table_name LIKE '%erass%'" if role == "ero" else "table_name LIKE '%2rxs%'")
        try:
            tables = discover_tables(arc, endpoint, cond)
            if role == "ero":
                label = "supp" if key == "ero_supp" else "main"
                tables = [(t, d) for t, d in tables if label in (t + " " + d).lower()
                          and "hard" not in (t + " " + d).lower()]
            else:
                tables.sort(key=lambda v: ("cat2rxs" not in v[0].lower(), len(v[0])))
            result = try_options([(get_tap_table, (endpoint, t, role)) for t, _ in tables])
            if result:
                return result
        except Exception as exc:
            arc.note(key + " discovery " + endpoint, exc)
    return None


# Morrison & McCammon (1983), Table 2. Coefficients apply to energies in keV.
EDGES = np.array([.03,.1,.284,.4,.532,.707,.867,1.303,1.840,2.471,3.210,4.038,7.111,8.331,10.])
COEF = np.array([(17.3,608.1,-2150.),(34.6,267.9,-476.1),(78.1,18.8,4.3),
                 (71.4,66.8,-51.4),(95.5,145.8,-61.1),(308.9,-380.6,294.),
                 (120.6,169.3,-47.7),(141.3,146.8,-31.5),(202.7,104.7,-17.),
                 (342.7,18.7,0.),(352.2,18.7,0.),(433.9,-2.4,.75),
                 (629.,30.9,0.),(701.2,25.2,0.)])


def photons(e, gamma, nh):
    e = np.asarray(e).clip(.03,10.)
    a,b,c = COEF[np.searchsorted(EDGES[1:-1], e, side="right")].T
    return e**(-gamma) * np.exp(-nh * 1e-24 * (a+b*e+c*e**2)/e**3)


def band(lower, upper, gamma, nh):
    e = np.geomspace(lower, upper, 3000)
    # ergs in 1 keV; arbitrary photon-spectrum normalisation cancels.
    return np.trapezoid(e*photons(e,gamma,nh),e) * 1.602176634e-9


def pspc_response(path):
    with fits.open(path) as h:
        m = next(x for x in h if "MATRIX" in list(getattr(getattr(x,"columns",None),"names",[]) or []))
        data = m.data
        lo, hi = data["ENERG_LO"].astype(float), data["ENERG_HI"].astype(float)
        area = []
        for row in data:
            fchan = np.atleast_1d(row["F_CHAN"]).astype(int)[:int(row["N_GRP"])]
            nchan = np.atleast_1d(row["N_CHAN"]).astype(int)[:int(row["N_GRP"])]
            channels = np.concatenate([np.arange(f,f+n) for f,n in zip(fchan,nchan)]) if len(fchan) else np.empty(0)
            vals = np.atleast_1d(row["MATRIX"]).astype(float)[:len(channels)]
            area.append(sum(vals[(channels>=11)&(channels<=235)]))
        area = np.asarray(area)
        if np.max(area)<5.:
            raise ValueError("response lacks a collecting area in cm2")
        return lo, hi, area


def calibrate(arc):
    root = BASE + "caldb/data/rosat/pspc/cpf/matrices/"
    files = ["pspcc_gain1_256.rsp", "pspcc_gain1_256.rsp.gz"]
    # PSPC-C was the detector used in the principal RASS scan epoch.
    try:
        listing = arc.fetch(root, "rosat-response-index.html").read_text()
        files = [s for s in re.findall(r'href=["\']([^"\']+)["\']',listing,re.I)
                 if re.search(r"pspc_?c.*256\.rsp(\.gz)?$",s,re.I)] + files
    except Exception as exc:
        arc.note("PSPC-C listing", exc)
    response = None
    response_url = None
    for file in dict.fromkeys(files):
        try:
            response_url = urljoin(root, file)
            response = pspc_response(arc.fetch(response_url, file.split("/")[-1]))
            break
        except Exception as exc:
            arc.note(file, exc)
    def ros(g,n):
        if response is None:
            return np.nan
        lo,hi,area = response
        denom = np.sum(photons((lo+hi)/2,g,n)*(hi-lo)*area)
        return band(.2,2.,g,n) / denom
    combos = list(itertools.product((1.5,2.,2.5,3.), (1e20,3e20,1e21)))
    c0 = ros(2.,3e20)
    def ero(g,n):
        return band(.2,2.,g,n)/band(.2,2.3,g,n)
    e0 = ero(2.,3e20)
    er = [ero(g,n)/e0 for g,n in combos]
    rr = [ros(g,n)/c0 for g,n in combos] if response is not None else [np.nan]
    def finite(x):
        return float(x) if np.isfinite(x) else None
    model = {"gamma":2., "NH_cm-2":3e20, "absorption":"Morrison & McCammon (1983) Table 2",
             "common_band_keV":[.2,2.],"ero_band_factor":e0,"rass_ecf":finite(c0),
             "rass_response":response_url if response is not None else None,
             "rass_channels":[11,235],
             "notes":"XMM catalogue flux is already 0.2–2 keV; native response assumptions are retained. eROSITA uses released 0.2–2.3 keV flux times the absorbed energy-flux band ratio. RASS uses an on-axis PSPC-C survey response. Point-source calibration and an illustrative grid do not cover arbitrary spectra, absorption, pile-up or spatial changes.",
             "model_grid":{"gamma":[1.5,2.,2.5,3.],"nh":[1e20,3e20,1e21]},
             "xmm_relative":[1.,1.], "ero_relative":[min(er),max(er)],
             "rass_relative":[finite(min(rr)),finite(max(rr))]}
    (arc.out / "spectral_model.json").write_text(json.dumps(model,indent=2))
    return model


def numeric(table,col,default=np.nan):
    if col not in table:
        return pd.Series(default,index=table.index,dtype=float)
    return pd.to_numeric(table[col],errors="coerce")


def truth(series):
    return series.astype(str).str.strip().str.lower().isin(["true","t","1","1.0","y","yes"])


def prepare(table,role,meta,model):
    table = table.copy()
    table["ra"] = numeric(table,"ra")
    table["dec"] = numeric(table,"dec")
    valid = table.ra.between(0,360,inclusive="left") & table.dec.between(-90,90)
    if valid.mean()<.99:
        raise ValueError(f"{role}: over 1% invalid positions; coordinate mapping or unit mismatch")
    table = table.loc[valid].reset_index(drop=True)
    for c in table.select_dtypes(include=["object","string"]):
        if c not in ("ra","dec"):
            table[c] = table[c].map(lambda v: v.strip("\0 ") if isinstance(v,str) else v)
    table["name"] = table.name.astype(str).str.strip()
    # Position floors allow for residual astrometric systematics. The match
    # radius caps bound contamination from low-count or pathological errors.
    floor, fallback = {"xmm":(6.,8.),"ero":(1.5,4.),"rass":(6.,20.)}[role]
    errpos = numeric(table,"pos")
    if role == "rass" and "posy" in table:
        # VizieR XERR/YERR are image pixels, not absolute astrometric errors.
        errpos[:] = np.nan
    table["stat_poserr"] = errpos
    table["pos_sigma"] = np.hypot(errpos.where(errpos.between(0,120), fallback),floor)
    table["poserr"] = table["pos_sigma"]
    if role == "rass":
        scale = model["rass_ecf"] if model["rass_ecf"] is not None else np.nan
        table["flux"] = numeric(table,"rate") * scale
        table["err"] = numeric(table,"errrate") * scale
    else:
        scale = 1. if role=="xmm" else model["ero_band_factor"]
        # Source-specific table units are required if a service scales flux columns.
        fu = meta["columns"]["flux"]["unit"]
        if role=="xmm" and norm(fu)=="10e12ergsscm2":
            # This non-standard legacy string literally reads "(10E-12)ergs/s/cm2".
            # The XMMSL2 release guide defines stored values in 10^-12 cgs.
            # The rate/flux consistency test below independently guards this factor.
            scale *= 1e-12
        elif fu not in ("None","nan",""):
            try:
                scale *= u.Unit(fu).to(u.erg/u.cm**2/u.s)
            except (ValueError, TypeError):
                if not ("erg" in fu.lower() and re.search(r"cm",fu,re.I)):
                    raise ValueError(f"unexpected {role} flux unit {fu}")
        table["flux"] = numeric(table,"flux") * scale
        table["err"] = numeric(table,"error") * scale
    good = (table.flux>0)&(table.err>0)&np.isfinite(table.flux)&np.isfinite(table.err)
    if role == "xmm":
        rate = numeric(table,"rate")
        calibration = (table.flux/rate).where(rate.gt(0)).median()
        if not 1.3e-12 < calibration < 1.6e-12:
            raise ValueError(f"XMMSL2 B6 count rate to flux inconsistent with 1.436e-12: {calibration}")
        good &= numeric(table,"ml").ge(10)
        if "ext" in table:
            good &= numeric(table,"ext").fillna(0).le(.5)
        for field in ("issue","veropt","verext","verhalo","verbg","verfalse","veredge"):
            if field in table:
                good &= ~truth(table[field])
    elif role == "ero":
        good &= numeric(table,"ml").ge(10)
        if "extlike" in table:
            good &= numeric(table,"extlike").fillna(0).lt(6)
        elif "ext" in table:
            good &= numeric(table,"ext").fillna(0).le(0)
        for c in table.columns:
            if c.startswith("flag_sp_"):
                good &= ~truth(table[c])
    else:
        good &= numeric(table,"ml").ge(10)
        if "flag" in table:
            good &= numeric(table,"flag").eq(0)
        if "near" in table:
            good &= numeric(table,"near").eq(0)
        if "region" in table:
            good &= numeric(table,"region").eq(1)
    table["quality"] = good
    return table


def join_mission(slew, target, mission):
    src = SkyCoord(slew.ra.to_numpy()*u.deg,slew.dec.to_numpy()*u.deg)
    dst = SkyCoord(target.ra.to_numpy()*u.deg,target.dec.to_numpy()*u.deg)
    tree = cKDTree(dst.cartesian.xyz.value.T)
    ss,ts = slew.pos_sigma.to_numpy(),target.pos_sigma.to_numpy()
    cap = 40 if mission=="ero" else 75
    chord = 2*np.sin(np.radians(cap/3600)/2)
    def evaluate(at):
        groups = tree.query_ball_point(at.cartesian.xyz.value.T, chord)
        lens = np.fromiter(map(len,groups),int,count=len(groups))
        ai = np.repeat(np.arange(len(groups)),lens)
        if not len(ai):
            return np.full(len(groups),-1),np.full(len(groups),np.nan),np.zeros(len(groups),int)
        bi = np.concatenate([g for g in groups if len(g)]).astype(int)
        sep = at[ai].separation(dst[bi]).arcsec
        keep = sep <= np.minimum(cap,3.44*np.hypot(ss[ai],ts[bi]))
        a,b,d = ai[keep],bi[keep],sep[keep]
        counts = np.bincount(a,minlength=len(groups))
        best = np.full(len(groups),-1)
        dist = np.full(len(groups),np.nan)
        # Last assignment to duplicate indexing is not a stable reduction.
        # Lexicographic order then first index ensures a deterministic nearest.
        order = np.lexsort((b,d,a))
        _, first = np.unique(a[order],return_index=True)
        ix = order[first]
        best[a[ix]],dist[a[ix]]=b[ix],d[ix]
        return best,dist,counts
    idx,dist,n = evaluate(src)
    observed = n>0
    for field in ["name","ra","dec","flux","err","quality","stat_poserr","poserr","pos_sigma","ml","rate","errrate","release"]:
        if field in target:
            values = target[field].iloc[np.maximum(idx,0)].reset_index(drop=True)
            slew[mission+"_"+field] = values.where(observed,False if field=="quality" else np.nan)
    slew[mission+"_sep"] = dist
    slew[mission+"_ambiguous"] = n>1
    slew[mission+"_nmatch"] = n
    # A ROSAT beam can contain distinct XMM entries. Identify mutually
    # competing slew positions while allowing repeated visits to one source.
    matched = pd.DataFrame({"idx":idx[observed],"group":slew.loc[observed,"source_group"].values})
    blends = matched.groupby("idx").group.nunique()
    blended = np.isin(idx,blends[blends>1].index)
    slew[mission+"_ambiguous"] |= blended
    controls = []
    for angle in range(0,360,45):
        shifted = src.directional_offset_by(angle*u.deg,20*u.arcmin)
        si,_,sn = evaluate(shifted)
        high = (np.maximum(slew.flux.to_numpy(),target.flux.to_numpy()[np.maximum(si,0)]) /
                np.minimum(slew.flux.to_numpy(),target.flux.to_numpy()[np.maximum(si,0)]))
        eligible = slew.quality.to_numpy() & target.quality.to_numpy()[np.maximum(si,0)] & (sn==1)
        controls.append({"pa_deg":angle,"shift_arcmin":20,"n":len(slew),"matches":int((sn>0).sum()),
                         "unique":int((sn==1).sum()),"unique_nominal_10x":int(((sn==1)&(high>=10)).sum()),
                         "quality_unique":int(eligible.sum()), "quality_unique_nominal_10x":int((eligible&(high>=10)).sum())})
    mean = float(np.mean([c["matches"] for c in controls]))
    return {"observed":int(observed.sum()),"observed_unambiguous":int((observed&~slew[mission+"_ambiguous"]).sum()),
            "mean_shifted":mean,"estimated_random_fraction":mean/int(observed.sum()) if observed.any() else None,
            "offsets":controls,"outer_radius_arcsec":cap,"adaptive_radius":"min(cap, 3.44 * hypot(position errors))"}


def catalogue_ratios(table):
    n = len(table)
    best = np.full(n,np.nan)
    best_lo = np.full(n,np.nan)
    pair = np.full(n,"",dtype=object)
    for a,b in itertools.combinations(("xmm","ero","rass"),2):
        fa,fb = table[a+"_flux"].to_numpy(float),table[b+"_flux"].to_numpy(float)
        ea,eb = table[a+"_err"].to_numpy(float),table[b+"_err"].to_numpy(float)
        valid = (fa>0)&(fb>0)&(ea>0)&(eb>0)
        with np.errstate(divide="ignore",invalid="ignore"):
            nominal = np.where(valid,np.maximum(fa/fb,fb/fa),np.nan)
            low = np.where(fa>=fb,(fa-1.96*ea)/(fb+1.96*eb),(fb-1.96*eb)/(fa+1.96*ea))
        low = np.where(valid,np.maximum(low,0),np.nan)
        prefix = a + "_" + b
        table[prefix+"_ratio"] = nominal
        table[prefix+"_ratio_lower"] = low
        table[prefix+"_brighter"] = np.where(valid,np.where(fa>=fb,a,b),"")
        improved = valid & (~np.isfinite(best_lo) | (low>best_lo))
        best[improved],best_lo[improved],pair[improved] = nominal[improved],low[improved],prefix
    table["ratio"],table["ratio_lower"],table["pair"] = best,best_lo,pair


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",required=True)
    p.add_argument("--out",required=True)
    args = p.parse_args()
    arc = Archive(args.work,args.out)
    model = calibrate(arc)
    meta, data = {},{}
    for key in ["xmmslewcln","xmmslewful","ero_main","ero_supp","rass"]:
        result = get_catalog(arc,key)
        if result is None:
            raise RuntimeError("A required catalogue is unavailable. See catalog_inputs.json: "+key)
        raw,meta[key] = result
        # Protect against changes to the rolling mission-service aliases.
        expected = {"xmmslewcln":29393, "xmmslewful":72352, "ero_main":930203, "ero_supp":347274}
        if key in expected and len(raw)!=expected[key]:
            raise ValueError(f"Release count mismatch {key} ({len(raw)} versus {expected[key]})")
        role = "xmm" if key.startswith("xmm") else ("rass" if key=="rass" else "ero")
        data[key] = prepare(raw,role,meta[key],model)
        if role=="ero":
            data[key]["release"] = key
        (arc.out / "catalog_schema.json").write_text(json.dumps(meta,indent=2))
        del raw
    clean,full = data.pop("xmmslewcln"),data.pop("xmmslewful")
    # Coordinates + epoch identify detections, independently of each catalogue's
    # row IDs and of source designations shared by multiple slew passages.
    def keys(t):
        return list(zip(t.ra,t.dec,t.get("obs",pd.Series("",index=t.index)).astype(str)))
    clean_keys = set(keys(clean))
    if not clean_keys.issubset(set(keys(full))):
        raise ValueError("XMMSL2 clean coordinates/obsids are not a subset of full")
    full["clean"] = [k in clean_keys for k in keys(full)]
    if int(full.clean.sum())!=len(clean):
        raise ValueError("ambiguous clean/full join, inspect FITS detection identifiers")
    full["quality"] &= full.clean
    full["source_group"] = full.get("unique",full.name).astype(str).str.strip()
    full["source_group"] = full.source_group.where(full.source_group.ne(""),full.name)
    full["detection_uid"] = ["XMMSL2-total-"+str(i) for i in range(1,len(full)+1)]
    controls = {}
    comparisons = {"ero":pd.concat([data["ero_main"],data["ero_supp"]],ignore_index=True),
                   "rass":data["rass"]}
    for mission,cat in comparisons.items():
        controls[mission] = join_mission(full,cat,mission)
    full.rename(columns={c:"xmm_"+c for c in ["rate","counts","exp","id","flux","err","quality","ml","obs","time","stat_poserr","poserr","pos_sigma",
                                              "ext","xtyp","counterpart","iddist","oldrass","altident","issue",
                                              "veropt","verext","verhalo","verbg","verfalse","veredge"] if c in full},inplace=True)
    for mission in ("xmm","ero","rass"):
        full[mission+"_model_min"],full[mission+"_model_max"] = model[mission+"_relative"]
    both = full.ero_name.notna() & full.rass_name.notna()
    er_coord = SkyCoord(full.loc[both,"ero_ra"].to_numpy()*u.deg, full.loc[both,"ero_dec"].to_numpy()*u.deg)
    ro_coord = SkyCoord(full.loc[both,"rass_ra"].to_numpy()*u.deg, full.loc[both,"rass_dec"].to_numpy()*u.deg)
    direct_sep = er_coord.separation(ro_coord).arcsec
    direct_gate = np.minimum(75., 3.44*np.hypot(full.loc[both,"ero_poserr"],full.loc[both,"rass_poserr"]))
    full["ero_rass_sep"] = np.nan
    full.loc[both,"ero_rass_sep"] = direct_sep
    full["ero_rass_consistent"] = False
    full.loc[both,"ero_rass_consistent"] = direct_sep <= direct_gate
    catalogue_ratios(full)
    full.to_csv(arc.out / "catalog_sources.csv.gz",index=False,compression={"method":"gzip","mtime":0},float_format="%.10g")
    (arc.out / "shifted_controls.json").write_text(json.dumps(controls,indent=2))
    stats = {"inputs":{k:v["rows"] for k,v in meta.items()},
             "clean_flagged":int(full.clean.sum()), "quality_slew":int(full.xmm_quality.sum()),
             "detection_rows":len(full),"unique_source_designations":int(full.source_group.nunique()),
             "positive_soft_slew":int((full.xmm_flux>0).sum()),
             "catalogue_only_nominal_10x":int((full[[p+"_ratio" for p in ["xmm_ero","xmm_rass","ero_rass"]]].max(axis=1)>=10).sum()),
             "control_unit":"XMMSL2 detection (includes repeat visits). Offset estimates are approximate and do not measure individual host probabilities."}
    (arc.out / "catalog_statistics.json").write_text(json.dumps(stats,indent=2))
    arc.save()
    print(json.dumps(stats),flush=True)


if __name__ == "__main__":
    main()
