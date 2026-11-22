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
                time="dateobs time starttime observationdate mjdstart mjdobs mjd",
                mjd_min="mjdstart", mjd_max="mjdstop",
                id="sourcenum detid detectionid uniquedetectionid",
                rate="rateb6", counts="sctsb6", exp="expmapb6", background_rate="imagebgrate",
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
                mjd="mjd", mjd_min="mjdmin", mjd_max="mjdmax",
                flux="mlflux1 mlflux0 flux0 broadflux fluxb flux broadbandflux",
                error="mlfluxerr1 mlfluxerr0 emlflux0 fluxerr0 eflux0 broadfluxerror broadfluxerr fluxerror fluxerr eflux",
                ml="detlike0 detectionlikelihood broadbandlikelihood detlike likelihood",
                ext="ext extent extent0 sourceextent",
                extlike="extlike extentlikelihood"),
    "rass": dict(F, name="name iauname 2rxs irocess rosat sourcename",
                 pos="poserror positionalerror poserrerr poserr errorradius positionerror positionuncertainty radecerror xerr",
                 posy="yerr",
                 rate="crates count_rate countrate crate srccountrate rate ct",
                 errrate="crateserr countrateerror countrateerr countrateerr1 ecrate ecrates crateerr erate rateerr",
                 ml="eximl existence_likelihood existlikelihood detectionlikelihood detectlikelihood detlikelihood likelihoodofthedetection likelihood",
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
    if role=="rass":
        for n in names:
            label=norm(n)
            desc=str(descriptions.get(n,"")).lower()
            if "ml" not in mapping and ("likelihood" in label or "likelihood" in desc):
                if ("det" in label or "exist" in label or "detection likelihood" in desc) and "ext" not in label:
                    mapping["ml"]=n
            if "pos" not in mapping and ("err" in label or "uncertainty" in label) and "pos" in label:
                mapping["pos"]=n
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
            mapping.update({n.lower(): n for n in units if n.upper().startswith("FLAG_SP_") or
                            n.upper() in {"FLAG_OPT","FLAG_NO_RADEC_ERR","FLAG_NO_CTS_ERR"}})
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
        options.extend([(get_tap_table, (NASA, "rass2rxs", role)), (get_tdat, ("rass2rxs", role))])
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
    for number in ("ml","rate","errrate","ext","extlike","counts","exp","mjd","mjd_min","mjd_max"):
        if number in table:
            table[number] = numeric(table,number)
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
            if c.startswith("flag_sp_") or c in {"flag_opt","flag_no_radec_err","flag_no_cts_err"}:
                good &= ~(truth(table[c]) | numeric(table,c).gt(0))
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


def join_mission(slew, target, mission, supplemental_states):
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
    for field in ["name","ra","dec","flux","err","quality","stat_poserr","poserr","pos_sigma","ml","rate","errrate","release","mjd","mjd_min","mjd_max","flag_opt","flag_no_cts_err","flag_no_radec_err"]:
        if field in target:
            values = target[field].iloc[np.maximum(idx,0)].reset_index(drop=True)
            slew[mission+"_"+field] = values.where(observed,False if field=="quality" else np.nan)
    slew[mission+"_sep"] = dist
    slew[mission+"_ambiguous"] = n>1
    slew[mission+"_nmatch"] = n
    # A ROSAT beam can contain distinct XMM entries. Identify mutually
    # competing slew positions while allowing repeated visits to one source.
    def ambiguous(indices, multiplicity):
        keep = multiplicity>0
        matches = pd.DataFrame({"idx":indices[keep],
                                "group":slew.loc[keep,"source_group"].values})
        blends = matches.groupby("idx").group.nunique()
        return (multiplicity>1) | np.isin(indices,blends[blends>1].index)

    def summarize(indices, multiplicity, confusion, separation, orientation):
        safe_idx = np.maximum(indices,0)
        fx,ex = slew.flux.to_numpy(),slew.err.to_numpy()
        ft,et = target.flux.to_numpy()[safe_idx],target.err.to_numpy()[safe_idx]
        reliable_target = (multiplicity==1) & ~confusion & target.quality.to_numpy()[safe_idx]
        good = reliable_target & slew.quality.to_numpy()
        valid = (fx>0)&(ft>0)&(ex>0)&(et>0)
        with np.errstate(divide="ignore",invalid="ignore"):
            nominal = np.maximum(fx/ft,ft/fx)
            xbright = (fx-1.96*ex)/(ft+1.96*et)
            tbright = (ft-1.96*et)/(fx+1.96*ex)
        supplemental_states[orientation] = {
            "secure":reliable_target & (numeric(target,"ml").to_numpy()[safe_idx]>=20) &
                     (separation<=(10 if mission=="ero" else 15)),
            "guarded":reliable_target&valid&((xbright>=10)|(tbright>=10))}
        return {"n":len(slew),"matches":int((multiplicity>0).sum()),
                "unique":int((multiplicity==1).sum()),
                "unique_nominal_10x":int(((multiplicity==1)&valid&(nominal>=10)).sum()),
                "quality_unique":int(good.sum()),"quality_unique_nominal_10x":int((good&valid&(nominal>=10)).sum()),
                "guarded_10x":int((good&valid&((xbright>=10)|(tbright>=10))).sum()),
                "guarded_10x_xmm_brighter":int((good&valid&(xbright>=10)).sum()),
                "guarded_10x_comparison_brighter":int((good&valid&(tbright>=10)).sum())}

    confusion = ambiguous(idx,n)
    slew[mission+"_ambiguous"] = confusion
    obs = summarize(idx,n,confusion,dist,"observed")
    # Named calibration / persistent or known-variable targets. Recovery means
    # only that the positional association machinery found the object.
    benchmarks = []
    for label,ra,dec in [("RX J1856.5-3754",284.1463,-37.9084),
                         ("Vela X-1",135.5286,-40.5547),("LMC X-3",84.7355,-64.0843)]:
        fixed = SkyCoord(ra*u.deg,dec*u.deg)
        separation = src.separation(fixed).arcsec
        k = int(np.argmin(separation))
        entry = {"reference":label,"ra":ra,"dec":dec,"search_arcsec":30,
                 "nearest_slew_sep":float(separation[k]),"slew_within_30arcsec":bool(separation[k]<=30)}
        if separation[k]<=30:
            entry.update(xmm_name=slew.name.iloc[k], positive_match=str(target.name.iloc[idx[k]]) if idx[k]>=0 else None,
                         match_sep=float(dist[k]) if idx[k]>=0 else None,
                         match_ambiguous=bool(confusion[k]), slew_row=k)
        benchmarks.append(entry)
    controls = []
    for angle in range(0,360,45):
        shifted = src.directional_offset_by(angle*u.deg,20*u.arcmin)
        si,sd,sn = evaluate(shifted)
        conf = ambiguous(si,sn)
        controls.append(dict(pa_deg=angle,shift_arcmin=20,**summarize(si,sn,conf,sd,angle)))
        if angle==90:
            for item in benchmarks:
                if "slew_row" not in item:
                    continue
                k = item.pop("slew_row")
                item["east_20arcmin_control"] = {
                    "ra":float(shifted[k].ra.deg),"dec":float(shifted[k].dec.deg),
                    "comparison_match":str(target.name.iloc[si[k]]) if si[k]>=0 else None,
                    "sep":float(sd[k]) if si[k]>=0 else None,"ambiguous":bool(conf[k])}
    def mean_shifted(k):
        return float(np.mean([c[k] for c in controls]))
    def fraction(k):
        return mean_shifted(k)/obs[k] if obs[k] else None
    return {"observed":obs["matches"],"observed_unambiguous":int((observed&~confusion).sum()),
            "mean_shifted":mean_shifted("matches"),"estimated_random_fraction":fraction("matches"),
            "observed_selection_counts":obs, "guarded_10x_mean_shifted":mean_shifted("guarded_10x"),
            "guarded_10x_estimated_random_fraction":fraction("guarded_10x"),
            "guarded_10x_shifted_range":[min(c["guarded_10x"] for c in controls),
                                       max(c["guarded_10x"] for c in controls)],
            "offset_dependence":"Offsets share the same sky and source sample. The eight-shift range describes empirical variation and is not a confidence interval.",
            "control_selection":"Same clean point-source/likelihood, positive-error, mutual-confusion and (Fhi-1.96 sigma_hi)/(Flo+1.96 sigma_lo)>=10 cuts as observed XMM versus this catalogue. Does not cover the map or ERO-RASS pair selections.",
            "benchmarks":benchmarks, "offsets":controls,"outer_radius_arcsec":cap,
            "adaptive_radius":"min(cap, 3.44 * hypot(position errors))"}



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
        expected = {"xmmslewcln":29393, "xmmslewful":72352, "ero_main":930203, "ero_supp":347274, "rass":135118}
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
    full["native_quality"] = full.quality.copy()
    flagcols = ["issue","veropt","verext","verhalo","verbg","verfalse","veredge"]
    full["no_warnings"] = True
    full["warning_fields"] = ""
    for field in flagcols:
        if field not in full:
            raise ValueError("Required XMMSL2 warning column absent: " + field)
        flagged = truth(full[field])
        full.loc[flagged,"no_warnings"] = False
        native = meta["xmmslewful"]["columns"][field]["field"]
        full.loc[flagged,"warning_fields"] += native + "|"
    full["warning_fields"] = full.warning_fields.str.rstrip("|")
    full["quality"] &= full.clean
    full["source_group"] = full.get("unique",full.name).astype(str).str.strip()
    full["source_group"] = full.source_group.where(full.source_group.ne(""),full.name)
    full["detection_uid"] = ["XMMSL2-total-"+str(i) for i in range(1,len(full)+1)]
    controls = {}
    comparisons = {"ero":pd.concat([data["ero_main"],data["ero_supp"]],ignore_index=True),
                   "rass":data["rass"]}
    supplemental_states = {m:{} for m in comparisons}
    for mission,cat in comparisons.items():
        controls[mission] = join_mission(full,cat,mission,supplemental_states[mission])
    # An explicitly labelled, higher-risk full-only stratum. Independent
    # visits use distinct OBSIDs, not repeated detections within one slewing exposure.
    positive = full.flux.gt(0) & full.err.gt(0)
    enough_counts = numeric(full,"counts").ge(5)
    adequate_snr = (full.flux / full.err).ge(3)
    strong = full.native_quality & enough_counts & adequate_snr
    def other_visits(use):
        keys = ["source_group","obs"]
        visits = full.loc[use,keys].drop_duplicates()
        n = full.source_group.map(visits.groupby("source_group").size()).fillna(0).astype(int)
        own = pd.MultiIndex.from_frame(full[keys]).isin(pd.MultiIndex.from_frame(visits))
        return n-own.astype(int)
    full["other_positive_slew_obsids"] = other_visits(positive)
    full["other_strong_slew_obsids"] = other_visits(strong)
    for m,cap in (("ero",10),("rass",15)):
        full["close_secure_"+m] = (full[m+"_quality"] & ~full[m+"_ambiguous"] &
                                   full[m+"_ml"].ge(20) & full[m+"_sep"].le(cap))
    both_secure = full.close_secure_ero & full.close_secure_rass
    full["full_only_corroborated"] = (~full.clean & strong &
                                      ((full.other_strong_slew_obsids>0) | both_secure))
    for m in comparisons:
        offset_counts = {}
        for angle in ["observed", *range(0,360,45)]:
            corroboration = (supplemental_states["ero"][angle]["secure"] &
                             supplemental_states["rass"][angle]["secure"])
            sample = (~full.clean & strong & (full.other_strong_slew_obsids.gt(0)|corroboration))
            offset_counts[str(angle)] = int((sample & supplemental_states[m][angle]["guarded"]).sum())
        mean = float(np.mean([offset_counts[str(a)] for a in range(0,360,45)]))
        controls[m]["corroborated_full_only_guarded_10x"] = {"counts":offset_counts,
            "mean_shifted":mean,"estimated_random_fraction":mean/offset_counts["observed"] if offset_counts["observed"] else None,
            "selection":"Recompute BOTH close mission counterparts at each identical 20-arcmin XMM shift before applying the full-only corroboration and guarded-amplitude rules; repeat-slew evidence is intrinsic and fixed."}
    full["full_only_reasons"] = ""
    reasons = {"no_positive_soft_flux_or_error":~positive,
               "soft_ml_below_10_or_missing":~numeric(full,"ml").ge(10),
               "extension_gt_0.5_pixel":numeric(full,"ext").gt(.5),
               "instrument_warning":~full.no_warnings,
               "net_counts_below_5_or_missing":~enough_counts,
               "flux_snr_below_3_or_missing":~adequate_snr,
               "no_strong_independent_slew_or_two_close_missions":~(full.other_strong_slew_obsids.gt(0)|both_secure)}
    for label,mask in reasons.items():
        full.loc[~full.clean & mask,"full_only_reasons"] += label+"|"
    full["full_only_reasons"] = full.full_only_reasons.str.rstrip("|")
    off = ~full.clean
    def count_off(mask):
        return int((off & mask).sum())
    full_only = {"rows":int(off.sum()), "positive_soft_flux":count_off(full.flux.gt(0)),
                 "soft_likelihood_at_least_10":count_off(numeric(full,"ml").ge(10)),
                 "soft_likelihood_at_least_15p5":count_off(numeric(full,"ml").ge(15.5)),
                 "no_instrument_warnings":count_off(full.no_warnings),"native_quality":count_off(full.native_quality),
                 "positive_with_distinct_positive_slew_obsid":count_off(positive & full.other_positive_slew_obsids.gt(0)),
                 "strong_soft_photometry":count_off(strong),"strong_with_distinct_strong_slew_obsid":count_off(strong & full.other_strong_slew_obsids.gt(0)),
                 "close_secure_ero":count_off(full.close_secure_ero), "close_secure_rass":count_off(full.close_secure_rass),
                 "close_secure_both":count_off(both_secure),"corroborated_rows":int(full.full_only_corroborated.sum()),
                 "corroborated_unique_slew_names":int(full.loc[full.full_only_corroborated,"source_group"].nunique()),
                 "failure_counts_nonexclusive":{reason:count_off(mask) for reason,mask in reasons.items()},
                 "corroboration":"Positive unflagged point source with soft ML≥10, ≥5 net counts and flux/error≥3 plus a distinct equally strong slew OBSID with same release name OR secure, unambiguous counterparts in BOTH ERO and 2RXS at ≤10/15 arcsec and detection likelihood≥20. These supplemental rows retain clean=False."}
    full.rename(columns={c:"xmm_"+c for c in ["rate","counts","exp","id","background_rate","full_only_corroborated","native_quality","no_warnings","warning_fields","full_only_reasons","other_positive_slew_obsids","other_strong_slew_obsids","close_secure_ero","close_secure_rass","flux","err","quality","ml","obs","time","mjd_min","mjd_max","stat_poserr","poserr","pos_sigma",
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
    def distribution(source, flux, error, quality):
        f,e = source[flux],source[error]
        good = source[quality].fillna(False).astype(bool)
        positive = f.gt(0)&e.gt(0)
        def median(column,mask):
            return float(column[mask].median()) if mask.any() else None
        return {"rows":len(source),"positive_flux":int(f.gt(0).sum()),"positive_flux_error":int(positive.sum()),
                "quality_rows":int(good.sum()),"median_positive_flux_cgs":median(f,positive),
                "median_positive_error_cgs":median(e,positive),"median_quality_flux_cgs":median(f,good)}
    qa = {"xmm_clean":distribution(clean,"flux","err","quality"),
          "xmm_full":distribution(full,"xmm_flux","xmm_err","xmm_quality")}
    for key in ("ero_main","ero_supp","rass"):
        qa[key] = distribution(data[key],"flux","err","quality")
    for mission in ("ero","rass"):
        matched = full.loc[full.clean & full[mission+"_name"].notna()]
        qa[mission+"_at_clean_slew"] = distribution(matched,mission+"_flux",mission+"_err",mission+"_quality")
        qa[mission+"_at_clean_slew"]["unambiguous_matches"] = int((~matched[mission+"_ambiguous"]).sum())
        qa[mission+"_at_clean_slew"]["distinct_counterpart_names"] = int(matched[mission+"_name"].nunique())
    stats = {"inputs":{k:v["rows"] for k,v in meta.items()}, "quality_distribution":qa, "full_only_audit":full_only,
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
