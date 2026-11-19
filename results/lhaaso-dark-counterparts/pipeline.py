"""Public-catalogue acquisition for the 1LHAASO dark-source experiment.

This script runs on the *GitHub Actions runner*, not the sandbox. Failed catalogues
produce explicit pending receipts; a missing catalogue is never a non-detection.
The original 1LHAASO CDS CSV bytes are frozen for replay.
"""
from __future__ import annotations

import datetime as dt
import gzip
import hashlib
import io
import json
import os
import re
import tarfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
CAT = OUT / "catalogs"
FROZEN = OUT / "frozen"
for p in (OUT, CAT, FROZEN):
    p.mkdir(parents=True, exist_ok=True)
RECEIPTS: dict = {}
SESSION = requests.Session()
SESSION.headers.update({"User-Agent": "Open-Space-public-astronomy-audit/1.0 (citation and hashes)"})
TAP = [
    "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync",
    "https://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync",
]

def utc():
    return dt.datetime.now(dt.timezone.utc).isoformat()

def digest(data):
    return hashlib.sha256(data).hexdigest()

def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, sort_keys=True, default=str) + "\n")

def get(url, *, params=None, timeout=130, max_bytes=150_000_000):
    last = None
    for attempt in range(3):
        try:
            with SESSION.get(url, params=params, timeout=(20, timeout), stream=True) as r:
                r.raise_for_status()
                chunks, n = [], 0
                for b in r.iter_content(1 << 20):
                    if b:
                        chunks.append(b)
                        n += len(b)
                        if n > max_bytes:
                            raise ValueError("catalogue exceeds configured byte cap")
                data = b"".join(chunks)
                if not data or b"Too many requests were received" in data[:200]:
                    raise ValueError("blank/rate-limited catalogue")
                return data, r.url
        except Exception as exc:
            last = repr(exc)
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"{url}: {last}")

def csv_table(raw):
    if raw.lstrip().lower().startswith(b"<!doctype") or raw.lstrip().startswith(b"<html"):
        raise ValueError("TAP returned HTML rather than a catalogue")
    f = pd.read_csv(io.BytesIO(raw), comment="#", low_memory=False)
    f.columns = [x.strip() for x in f.columns]
    if len(f.columns) <= 2:
        raise ValueError(f"unexpected TAP CSV columns: {list(f.columns)}")
    return f

def vizier(key, table, min_rows=1):
    """Prefer independent live CfA/IDIA VizieR mirrors while Strasbourg is down."""
    urls = ["https://vizier.cfa.harvard.edu/viz-bin/asu-csv",
            "https://vizier.idia.ac.za/viz-bin/asu-csv"]
    errors=[]
    for host in urls:
        try:
            response=SESSION.get(host,params={"-source":table,"-out.all":"1",
                                              "-out.max":"20000"},
                                 timeout=(9,35))
            response.raise_for_status()
            raw=response.content
            if len(raw)>15_000_000:
                raise ValueError("mirror response too large; use targeted query")
            tab=csv_table(raw)
            if len(tab)<min_rows:
                raise ValueError(f"mirror row count {len(tab)} < {min_rows}")
            RECEIPTS[key]={"status":"loaded","table":table,"url":response.url,
                "mirror":"CfA/IDIA", "retrieved_utc":utc(),
                "raw_sha256":digest(raw),"bytes":len(raw),"rows":len(tab),
                "columns":list(tab.columns)}
            if key=="lhaaso":
                (FROZEN/"1lhaaso_cao2024_vizier.csv").write_bytes(raw)
                RECEIPTS[key]["frozen_path"]="frozen/1lhaaso_cao2024_vizier.csv"
                RECEIPTS[key]["input_variant"]="vizier_2024_cds_mirror"
            return tab
        except Exception as exc: errors.append(repr(exc)[:400])
    RECEIPTS[key]={"status":"pending","table":table,"errors":errors}
    return None


def github_preprint_fallback():
    """Author-associated Gammapy transcription of Cao's arXiv 2023 table.

    This is NOT VizieR 2024 data. Keep its distinct provenance and a pending
    primary-release comparison, even if counts reproduce 90/43.
    """
    key="lhaaso"
    url="https://raw.githubusercontent.com/gammapy/gammapy-data/main/catalogs/lhaaso_cat.csv"
    raw=file_url("gammapy_2023_preprint_table",url,max_bytes=500_000)
    if raw is None:return None
    try:
        original=pd.read_csv(io.BytesIO(raw),dtype=str,keep_default_na=True)
        if len(original)!=180 or len(original.columns)<11:
            raise ValueError(f"unexpected transcription shape {original.shape}")
        name=original.iloc[:,0].astype("string").str.strip().replace("",pd.NA).ffill()
        name=(name.str.replace("$+$","+",regex=False)
                  .str.replace("$-$","-",regex=False)
                  .str.replace("*","",regex=False)
                  .str.replace("1LHAASO","1LHAASO ",regex=False)
                  .str.replace(r"\s+"," ",regex=True).str.strip())
        result=pd.DataFrame({
          "1LHAASO":name,
          "f_1LHAASO":name.str.contains(r"u$",regex=True).map({True:"u",False:""}),
          "Comp":original.iloc[:,1].astype(str).str.replace("*","",regex=False).str.strip().str.upper(),
          "RAJ2000":original.iloc[:,2].map(number),
          "DEJ2000":original.iloc[:,3].map(number),
          "ePos":original.iloc[:,4].map(number),
          "r39":original.iloc[:,5].map(number),
          "TS":original.iloc[:,6].map(number),
          "N0":original.iloc[:,7].map(number),
          "Index":original.iloc[:,8].map(number),
          "TS100":original.iloc[:,9].map(number),
          "Assoc":original.iloc[:,10].str.strip().replace("",pd.NA)
        })
        if result["1LHAASO"].nunique()!=90:
            raise ValueError(f"wrong unique names: {result['1LHAASO'].nunique()}")
        (FROZEN/"1lhaaso_gammapy_2023_source.csv").write_bytes(raw)
        RECEIPTS[key]={
            "status":"loaded","url":url,"retrieved_utc":utc(),
            "input_variant":"gammapy_2023_preprint_transcription",
            "source_original_cds_verified":False,
            "source_original_cds_status":"pending_during_CDS_outage",
            "frozen_path":"frozen/1lhaaso_gammapy_2023_source.csv",
            "raw_sha256":digest(raw),"bytes":len(raw),"rows":len(original),
            "unique_names":result["1LHAASO"].nunique(),
            "columns":list(original.columns),
            "note":"Cao arXiv 2023 table transcribed by Gammapy, not 2024 VizieR original; no original-byte equivalence claimed."
        }
        return result
    except Exception as exc:
        RECEIPTS[key]={"status":"invalid_schema","source":url,"error":repr(exc)}
        return None

def file_url(key, url, *, max_bytes=150_000_000):
    try:
        raw, resolved = get(url, max_bytes=max_bytes)
        RECEIPTS[key] = {
            "status": "loaded", "url": resolved, "retrieved_utc": utc(),
            "raw_sha256": digest(raw), "bytes": len(raw)
        }
        return raw
    except Exception as exc:
        RECEIPTS[key] = {"status": "pending", "url": url, "error": str(exc)[:650]}
        return None

def col(df, *names, default=np.nan):
    for n in names:
        if n in df.columns:
            return df[n]
    return pd.Series([default] * len(df), index=df.index)

def number(x):
    try:
        if x is None or pd.isna(x): return float("nan")
        if isinstance(x, (int, float)): return float(x)
        m = re.search(r"[-+]?\d*\.?\d+(?:[Ee][-+]?\d+)?", str(x))
        return float(m.group()) if m else float("nan")
    except Exception:
        return float("nan")

def galactic_to_equatorial(l, b):
    v = SkyCoord(l=np.asarray(l)*u.deg, b=np.asarray(b)*u.deg, frame="galactic").icrs
    return np.asarray(v.ra.deg), np.asarray(v.dec.deg)

def coords(df, *, ra=("RAJ2000", "RAdeg", "RA_ICRS", "RAJ", "RAJD"), dec=("DEJ2000", "DEdeg", "DE_ICRS", "DEJ", "DECJD")):
    a, d = col(df, *ra), col(df, *dec)
    av = pd.to_numeric(a, errors="coerce").to_numpy()
    dv = pd.to_numeric(d, errors="coerce").to_numpy()
    if np.isfinite(av).sum() and np.isfinite(dv).sum():
        return av, dv
    # Some catalogue fields are sexagesimal, but never parse numerical degrees as hours.
    try:
        sky = SkyCoord(a.astype(str).to_list(), d.astype(str).to_list(),
                       unit=(u.hourangle, u.deg), frame="icrs")
        return np.asarray(sky.ra.deg), np.asarray(sky.dec.deg)
    except Exception:
        return av, dv

def records(items, fname):
    frame = pd.DataFrame(items)
    frame.to_csv(CAT / fname, index=False, float_format="%.9g")
    return frame

def base_record(name, ra, dec, size=0.0, **kw):
    if not np.isfinite(ra) or not np.isfinite(dec) or abs(dec) > 90:
        return None
    c = SkyCoord(ra*u.deg, dec*u.deg, frame="icrs").galactic
    return {"name": str(name).strip(), "ra": float(ra), "dec": float(dec),
            "l": float(c.l.deg), "b": float(c.b.deg), "radius_deg": float(size or 0),
            **kw}

def lhaaso():
    tab = vizier("lhaaso", "J/ApJS/271/25/catalog", min_rows=170)
    if tab is None:
        tab = github_preprint_fallback()
    if tab is None:
        return None
    if "1LHAASO" not in tab or "Comp" not in tab:
        RECEIPTS["lhaaso"]["status"] = "invalid_schema"
        tab = github_preprint_fallback()
        if tab is None: return None
    src = []
    for name, rows in tab.groupby("1LHAASO", dropna=True, sort=True):
        rows = rows.copy()
        rows["_ts"] = pd.to_numeric(col(rows, "TS"), errors="coerce")
        rows["_ra"] = pd.to_numeric(col(rows, "RAJ2000"), errors="coerce")
        rows["_dec"] = pd.to_numeric(col(rows, "DEJ2000"), errors="coerce")
        measured = rows[rows["_ra"].notna() & rows["_dec"].notna()]
        if measured.empty:
            raise AssertionError(f"no measured component: {name}")
        # Source name follows the most-significant instrument component.
        winner = measured.sort_values("_ts", ascending=False).iloc[0]
        assoc = sorted(set(x for x in rows["Assoc"].dropna().astype(str).str.strip()
                           if x and x.lower() != "nan")) if "Assoc" in rows else []
        flags = "".join(rows.get("f_1LHAASO", pd.Series(dtype=str)).fillna("").astype(str))
        flag_u = "u" in flags or str(name).rstrip().endswith("u")
        km = rows[(rows.Comp.astype(str).str.upper() == "KM2A") & rows["_ra"].notna()]
        k = km.sort_values("_ts", ascending=False).iloc[0] if len(km) else None
        val = base_record(str(name).strip() + ("u" if flag_u and not str(name).endswith("u") else ""),
                          number(winner["RAJ2000"]), number(winner["DEJ2000"]),
                          size=0,
                          r39_deg=max(0.0, number(winner.get("r39", 0))) if np.isfinite(number(winner.get("r39", 0))) else 0,
                          r95_deg=max(0.0, number(winner.get("ePos", 0))) if np.isfinite(number(winner.get("ePos", 0))) else 0,
                          component=str(winner["Comp"]),
                          ts=number(winner.get("TS", np.nan)),
                          km2a_ts100=number(k.get("TS100", np.nan)) if k is not None else float("nan"),
                          km2a_n0=number(k.get("N0", np.nan)) if k is not None else float("nan"),
                          km2a_index=number(k.get("Index", np.nan)) if k is not None else float("nan"),
                          uhe=flag_u, known_tev=bool(assoc),
                          tev_associations=";".join(assoc),
                          n_components=len(rows))
        if val is not None:
            src.append(val)
    s = records(src, "1lhaaso_sources.csv")
    # Fail closed if the frozen input is not the published 90/43.
    assert len(s) == 90, f"distinct source count {len(s)} != 90"
    assert int(s.uhe.sum()) == 43, f"UHE source count {int(s.uhe.sum())} != 43"
    RECEIPTS["lhaaso"]["unique_sources"] = len(s)
    RECEIPTS["lhaaso"]["uhe_sources"] = int(s.uhe.sum())
    RECEIPTS["lhaaso"]["published_tev_unassociated"] = int((~s.known_tev).sum())
    return s

def fermi():
    url = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/14yr_catalog/gll_psc_v35.fit"
    raw = file_url("fermi_4fgl_dr4_v35", url)
    if raw is None: return
    try:
        with fits.open(io.BytesIO(raw), memmap=False) as hdus:
            t = pd.DataFrame(np.array(hdus[1].data).byteswap().view(np.ndarray)) if False else hdus[1].data
            result = []
            for r in t:
                name = r["Source_Name"]
                if isinstance(name, bytes): name = name.decode(errors="replace")
                x = base_record(name, float(r["RAJ2000"]), float(r["DEJ2000"]),
                                size=0,
                                cls=str(r["CLASS1"]).strip() if "CLASS1" in t.names else "",
                                assoc=str(r["ASSOC1"]).strip() if "ASSOC1" in t.names else "",
                                data_release=int(r["DataRelease"]) if "DataRelease" in t.names else None)
                if x: result.append(x)
        if len(result) < 7000: raise ValueError(f"4FGL contains only {len(result)} valid sources")
        records(result, "fermi_dr4.csv")
        RECEIPTS["fermi_4fgl_dr4_v35"]["sources"] = len(result)
    except Exception as exc:
        RECEIPTS["fermi_4fgl_dr4_v35"]["status"] = "invalid_schema"
        RECEIPTS["fermi_4fgl_dr4_v35"]["error"] = repr(exc)

def snr_radius(size):
    # Green/SNR candidates: angular sizes often in arcminutes, e.g. '30x20'.
    if pd.isna(size): return 0.0
    nums = re.findall(r"\d+(?:\.\d+)?", str(size))
    return min(2.0, max([float(v) for v in nums], default=0.0) / 120)

def green():
    """Official Cambridge Green 2024 Oct summary, not a VizieR proxy."""
    from bs4 import BeautifulSoup
    url="https://www.mrao.cam.ac.uk/projects/surveys/snrs/snrs.data.html"
    raw=file_url("green_2024",url,max_bytes=1_000_000)
    if raw is None:return
    try:
        soup=BeautifulSoup(raw,"html.parser")
        text=soup.get_text()
        pat=r"(?m)^\s*(\d{1,3}\.\d+)\s+([+-]\d+\.\d+)\s+(\d{2})\s+(\d{2})\s+(\d{2})\s+([+-]\d{2})\s+(\d{2})\s+([0-9.x?]+)"
        result=[]
        for m in re.finditer(pat,text):
            l=float(m.group(1)); b=float(m.group(2))
            size=m.group(8)
            ra,dec=galactic_to_equatorial(l,b)
            x=base_record("G"+m.group(1)+m.group(2),float(ra),float(dec),
                          size=snr_radius(size),catalogue="Green2024Oct",
                          size_original_arcmin=size)
            if x:result.append(x)
        if len(result)!=310:
            raise ValueError(f"Green 2024 page parser recovered {len(result)} rather than 310 SNRs")
        records(result,"green_snrs.csv")
        RECEIPTS["green_2024"]["sources"]=len(result)
        RECEIPTS["green_2024"]["version"]="2024 October, summary of all 310 Green SNRs"
    except Exception as exc:
        RECEIPTS["green_2024"]["status"]="invalid_schema"
        RECEIPTS["green_2024"]["error"]=repr(exc)


def meerkat_snrcat():
    # These are *radio-SNR candidates*, not a 510k-row compact-source census.
    for key,table,fn,minrow in [
        ("smgps_new_snr_candidates","J/A+A/693/A247/meerkatc","smgps_candidates.csv",230),
        ("snrcat_smgps_overlap","J/A+A/693/A247/ksnrcat","smgps_snrcat.csv",25),
    ]:
        t=vizier(key,table,min_rows=minrow)
        if t is None: continue
        rows=[]
        for _,r in t.iterrows():
            l=number(r.get("GLON",np.nan)); b=number(r.get("GLAT",np.nan))
            if not np.isfinite(l) or not np.isfinite(b): continue
            ra,dec=galactic_to_equatorial(l,b)
            x=base_record(str(r.get("Name",r.get("SNR",""))),float(ra),float(dec),
                          size=number(r.get("Radius",r.get("Rad",0)))/60 if np.isfinite(number(r.get("Radius",r.get("Rad",0)))) else 0,
                          catalogue=key)
            if x: rows.append(x)
        records(rows,fn)
        RECEIPTS[key]["normalized"] = len(rows)

def erass_hard():
    """Official eROSITA-DE DR1 hard 2.3-5 keV FITS; only German half sky."""
    url="https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Hard.tar.gz"
    raw=file_url("erass1_dr1_hard",url,max_bytes=15_000_000)
    if raw is None:return
    try:
        with tarfile.open(fileobj=io.BytesIO(raw),mode="r:gz") as tar:
            names=[n for n in tar.getmembers() if n.name.endswith(".fits")]
            if not names:raise ValueError("No .fits member in official hard archive")
            fitsbytes=tar.extractfile(names[0]).read()
        with fits.open(io.BytesIO(fitsbytes),memmap=False) as h:
            t=h[1].data
            out=[]
            for r in t:
                name=r["IAUNAME"]
                if isinstance(name,bytes):name=name.decode(errors="replace")
                x=base_record(name,float(r["RA"]),float(r["DEC"]),
                              size=max(0,float(r["EXT"])/3600) if np.isfinite(r["EXT"]) else 0,
                              catalogue="eRASS1_Hard_DR1",flux_2_3_5=float(r["ML_FLUX_0"]) if np.isfinite(r["ML_FLUX_0"]) else float("nan"))
                if x:out.append(x)
        if not 5000<=len(out)<=7000:raise ValueError(f"unexpected hard catalogue count {len(out)}")
        records(out,"erass1_hard.csv")
        RECEIPTS["erass1_dr1_hard"].update({"sources":len(out),"fits_sha256":digest(fitsbytes),
         "scope":"eROSITA-DE public western half sky only; approximate l=180..360 Galactic"})
    except Exception as exc:
        RECEIPTS["erass1_dr1_hard"]["status"]="invalid_schema"
        RECEIPTS["erass1_dr1_hard"]["error"]=repr(exc)


def gaia_young():
    t=vizier("gaia_dr3_hunt_ref", "J/A+A/673/A114/clusters", min_rows=7000)
    if t is None: return
    ra,dec=coords(t,ra=("RA_ICRS","RAJ2000"),dec=("DE_ICRS","DEJ2000"))
    result=[]
    for i,(a,d) in enumerate(zip(ra,dec)):
        r=t.iloc[i]
        age=number(r.get("logAge50",np.nan))
        typ=str(r.get("Type","")).strip()
        # Open clusters only. Young is <10 Myr (logAge<7), never all Gaia stars.
        if typ != "o" or not np.isfinite(age) or age>=7: continue
        x=base_record(str(r.get("Name","")),a,d,size=max(0,number(r.get("rt",0))) if np.isfinite(number(r.get("rt",0))) else 0,
                      log_age=age,distance_pc=number(r.get("dist50",np.nan)))
        if x: result.append(x)
    records(result,"gaia_young_clusters.csv")
    RECEIPTS["gaia_dr3_hunt_ref"]["young_open_clusters"] = len(result)

def atnf():
    """Parse the published ATNF psrcat.db directly, freezing its archive digest.

    The ATNF database records are key/value lines separated by '@'; all derived
    Edot/characteristic ages use conventional I=10^45 cgs assumptions. Missing
    DM distances remain unknown rather than filled with fabricated distances.
    """
    url="https://www.atnf.csiro.au/research/pulsar/psrcat/downloads/psrcat_pkg.tar.gz"
    blob=file_url("atnf_psrcat_v2_latest",url,max_bytes=40_000_000)
    if blob is None:return
    try:
        with tarfile.open(fileobj=io.BytesIO(blob),mode="r:gz") as tf:
            paths=[m for m in tf.getmembers() if m.isfile() and m.name.endswith("/psrcat.db")]
            if len(paths)!=1:raise ValueError(f"database archive paths {len(paths)}")
            db=tf.extractfile(paths[0]).read()
        from astropy.coordinates import BarycentricTrueEcliptic
        sources=[]
        current={}
        def flush():
            if not current:return
            name=current.get("PSRJ",current.get("JNAME",current.get("NAME","")))
            if not name:return
            ra,dec=float("nan"),float("nan")
            if "RAJD" in current and "DECJD" in current:
                ra=number(current["RAJD"]);dec=number(current["DECJD"])
            elif "RAJ" in current and "DECJ" in current:
                try:
                    sky=SkyCoord(current["RAJ"],current["DECJ"],unit=(u.hourangle,u.deg),frame="icrs")
                    ra=float(sky.ra.deg);dec=float(sky.dec.deg)
                except (ValueError,TypeError):pass
            elif "ELONG" in current and "ELAT" in current:
                try:
                    sky=SkyCoord(lon=number(current["ELONG"])*u.deg,
                                 lat=number(current["ELAT"])*u.deg,
                                 frame=BarycentricTrueEcliptic).icrs
                    ra=float(sky.ra.deg);dec=float(sky.dec.deg)
                except (ValueError,TypeError):pass
            p0=number(current.get("P0"))
            p1=number(current.get("P1"))
            f0=number(current.get("F0"))
            f1=number(current.get("F1"))
            if not np.isfinite(p0) and np.isfinite(f0) and f0>0:p0=1/f0
            if not np.isfinite(p1) and np.isfinite(f1) and np.isfinite(f0) and f0>0:
                p1=-f1/(f0*f0)
            edot=number(current.get("EDOT"))
            if not np.isfinite(edot) and np.isfinite(p0) and p0>0 and np.isfinite(p1) and p1>0:
                edot=4*np.pi**2*1e45*p1/p0**3
            age=number(current.get("AGE"))
            if not np.isfinite(age) and np.isfinite(p0) and p0>0 and np.isfinite(p1) and p1>0:
                age=p0/(2*p1)/(365.25*86400)
            d=float("nan"); ds="no_reported_distance"
            for candidate in ("DIST_A","DIST","DIST_DM","DIST_DM1"):
                x=number(current.get(candidate))
                if np.isfinite(x) and x>0:
                    d=x;ds=candidate;break
            if not np.isfinite(d):
                px=number(current.get("PX"))
                if np.isfinite(px) and px>0:
                    d=1/px;ds="inverse_parallax_unbiased_not_justified"
            val=base_record(name,ra,dec,edot_erg_s=edot,
                            distance_kpc=d,distance_source=ds,
                            age_yr=age,period_s=p0,period_derivative=p1,
                            age_model="characteristic_P_2Pdot",
                            edot_model="4pi2_I_Pdot_over_P3_I1e45")
            if val:sources.append(val)
        for rawline in db.decode("ascii",errors="replace").splitlines():
            line=rawline.strip()
            if not line or line.startswith("#"):continue
            if line.startswith("@"):
                flush();current={}
                continue
            parts=line.split()
            if len(parts)>=2 and re.match(r"^[A-Z_][A-Z0-9_]*$",parts[0]):
                current[parts[0]]=parts[1]
        flush()
        if not (3000<=len(sources)<=7000):
            raise ValueError(f"ATNF db parser yielded {len(sources)} positions, expected 3000-7000")
        df=records(sources,"atnf_pulsars.csv")
        if not df["name"].is_unique:
            raise ValueError("duplicate PSRJ identifiers")
        first=db[:2500].decode("ascii",errors="replace")
        version=re.search(r"(?i)(?:version|ver)\s*:?\s*v?(\d+\.\d+(?:\.\d+)?)",first)
        RECEIPTS["atnf_psrcat_v2_latest"].update({
            "psrcat_db_sha256":digest(db),
            "version":version.group(1) if version else "not_in_first_header",
            "pulsars_with_coordinates":len(sources),
            "pulsars_with_edot":int(np.isfinite(df.edot_erg_s).sum()),
            "pulsars_with_distance":int(np.isfinite(df.distance_kpc).sum()),
            "parser":"direct published ATNF psrcat.db record decoder, no cached psrqpy calculations"})
    except Exception as exc:
        RECEIPTS["atnf_psrcat_v2_latest"]["status"]="invalid_schema"
        RECEIPTS["atnf_psrcat_v2_latest"]["error"]=repr(exc)


def hess():
    # HGPS 2018 is a legacy veto catalogue; do not describe it as post-2024.
    url="https://hess.in2p3.fr/wp-content/uploads/2025/09/hgps_catalog_v1.fits.gz.zip"
    raw=file_url("hess_hgps_legacy",url,max_bytes=20_000_000)
    if raw is None:return
    try:
        import zipfile
        with zipfile.ZipFile(io.BytesIO(raw)) as z:
            names=[n for n in z.namelist() if n.endswith(".fits.gz")]
            if not names: raise ValueError("no FITS GZ in archive")
            zbytes=z.read(names[0])
            arr=zbytes if zbytes.startswith(b"SIMPLE") else gzip.decompress(zbytes)
        with fits.open(io.BytesIO(arr)) as h:
            lines=[]
            for hd in h[1:]:
                if not getattr(hd,"columns",None): continue
                for r in hd.data:
                    names=hd.columns.names
                    if "RA" in names and "DEC" in names:
                        a=number(r["RA"]);d=number(r["DEC"])
                    elif "GLON" in names and "GLAT" in names:
                        a,d=galactic_to_equatorial(number(r["GLON"]),number(r["GLAT"]))
                    else: continue
                    name=str(r["Source_Name"]) if "Source_Name" in names else str(r["Name"]) if "Name" in names else ""
                    x=base_record(name,float(a),float(d),catalogue="HESS-HGPS")
                    if x:lines.append(x)
                if len(lines)>50:break
        if len(lines)<50: raise ValueError(f"HGPS schema yielded {len(lines)} positions")
        records(lines,"hess_hgps.csv")
        RECEIPTS["hess_hgps_legacy"]["sources"]=len(lines)
    except Exception as exc:
        RECEIPTS["hess_hgps_legacy"]["status"]="invalid_schema"
        RECEIPTS["hess_hgps_legacy"]["error"]=str(exc)

def pending_only():
    # These require a validated survey footprint and a machine-readable
    # counterpart-specific sample, not a blanket sky non-detection.
    RECEIPTS["erass1_dr1_main_soft"]={
      "status":"pending", "reason":"930203-source catalogue not ingested; eRASS1-hard-only results do not rule out soft counterparts; half-sky rights mask required.",
      "url":"https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/"}
    RECEIPTS["meerkat_smgps_compact_510599"]={
      "status":"pending", "reason":"Dense 510599-source radio compact catalogue needs footprint/extended-source classification; SMGPS SNR candidate table is NOT this catalogue.",
      "url":"https://doi.org/10.48479/h16d-2z55"}
    RECEIPTS["askap_emu"]={
      "status":"pending", "reason":"No validated release-dependent mosaic footprint or uniform source-selection query in this run; do not claim non-detections."}
    RECEIPTS["snrcat_full"]={
      "status":"pending", "reason":"SMGPS SNRcat overlap subset is NOT the full live SNRcat; external SNRcat coverage/version requires independent freeze.",
      "url":"https://www.samarsafi-harb.com/snrcat"}
    RECEIPTS["hawc_latest"]={
      "status":"pending", "reason":"Newest HAWC release and machine-readable authoritative catalogue not frozen; older 3HWC is not a post-2024 discovery.",
      "url":"https://data.hawc-observatory.org/datasets/3hwc-survey/catalog.php"}

def main():
    s=lhaaso()
    if s is None:
        pending_only()
        write_json(OUT/"receipts.json",RECEIPTS)
        write_json(OUT/"run_status.json",{"status":"primary_catalogue_unavailable","time_utc":utc(),"claims_allowed":False})
        raise SystemExit("Failed to freeze/validate the primary Cao 2024 catalogue; no census will be claimed.")
    stages=[
      ("fermi",fermi),("atnf",atnf),("green",green),
      ("smgps SNR",meerkat_snrcat),("eROSITA hard",erass_hard),
      ("Gaia young clusters",gaia_young),("HESS legacy",hess)
    ]
    for label,fn in stages:
        try:
            fn()
        except Exception as exc:
            RECEIPTS["stage_"+label]={"status":"pending","error":repr(exc)}
    pending_only()
    write_json(OUT/"receipts.json",RECEIPTS)
    write_json(OUT/"run_status.json",{
        "status":"primary_validated_counterpart_catalogues_partial",
        "time_utc":utc(), "claims_allowed":False,
        "loaded":sorted(k for k,v in RECEIPTS.items() if v.get("status")=="loaded"),
        "pending":sorted(k for k,v in RECEIPTS.items() if v.get("status")!="loaded")})
    print("PRIMARY SOURCES:",len(s),"UHE:",s.uhe.sum())
    print("CATALOGUE RECEIPTS:",{k:v.get("status") for k,v in RECEIPTS.items()})

if __name__=="__main__": main()
