"""Coordinate and redshift crosschecks. All archive calls are made by the Actions runner."""
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import io
import json
from itertools import combinations
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import zipfile
from xml.etree import ElementTree

from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky
from astropy.io import fits
from bs4 import BeautifulSoup
import numpy as np
import pandas as pd
import requests

LENS_COLS = ["catalog", "name", "ra", "dec", "classification", "benchmark"]
PUBLIC_COLS = ["row_id", "catalog", "name", "era", "edec", "classification", "z",
               "z_err", "warning", "distance_arcsec", "redshift_ref"]
PUBLIC = ["Milliquas v8", "DESI DR1", "SDSS DR19"]
TAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
GLQ = "https://research.ast.cam.ac.uk/lensedquasars/"
LENS_URL = {
    "GraL X": "J/A+A/707/A345/tablea1",
    "Gaia FPR 381": "J/A+A/685/A130/table3",
    "Lemon IV": "J/MNRAS/520/3305/table1",
    "DESI Hsu spectroscopic pairs": "https://zenodo.org/api/records/17153224",
    "DESI DR1 strong lens VAC": "https://data.desi.lbl.gov/public/dr1/vac/dr1/strong-lensing/v1.0/desi-sl-vac-v1.fits",
    "DESI He": "J/A+A/672/A123/table1",
    "Lenscat": "https://raw.githubusercontent.com/lenscat/lenscat/main/lenscat/data/catalog.csv",
    "DESI Jing quasar pairs": "https://api.github.com/repos/astroliang/DESI_DR1_QP",
    "GraL III shortlist": "https://orbi.uliege.be/bitstream/2268/239586/1/Delchambre_et%20al_III.pdf",
    "GraL VIII radio census": "https://uhra.herts.ac.uk/id/eprint/11027/1/stad4002.pdf",
    "DESI Dawes": "J/ApJS/269/61/table2",
    "DULAG": "J/A+A/692/A154/superset",
    "GraL VI asterisms": "J/ApJ/921/42/table6",
    "GLQ doubles": GLQ + "doubles.html",
    "GLQ quads": GLQ + "quads.html",
    "GLQ binaries": GLQ + "binaries.html",
    "GLQ NIQs": GLQ + "niqs.html",
    "CASTLES": "https://lweb.cfa.harvard.edu/castles/noimages.html",
    "GraL IX appendix": "https://uhra.herts.ac.uk/id/eprint/25641/2/aa51690-24.pdf",
    "Euclid Q1": "https://zenodo.org/records/15025832/files/q1_discovery_engine_lens_catalog.csv?download=1",
    "Euclid DR1 ECLIPSE": "https://zenodo.org/records/19654904/files/ECLIPSE_object_catalogue.zip?download=1",
}
AUDIT = []
HEADERS = {"User-Agent": "Open-Space-Gaia-FPR/1.0 (public astronomical catalog cross-match)"}


def _download(cat, url, params=None, method="GET", files=None, timeout=55, headers=None):
    result = requests.request(method, url, params=params if method == "GET" else None,
                              data=params if method == "POST" else None, files=files,
                              headers={**HEADERS, **(headers or {})}, timeout=(12, timeout))
    result.raise_for_status()
    raw = result.content
    AUDIT.append(dict(catalog=cat, url=result.url, query=str(params or ""),
                      bytes=len(raw), sha256=hashlib.sha256(raw).hexdigest()))
    return raw


def _column(frame, names, required=True):
    # Ignore case and separators, retain explicit column preference.
    mapping = {re.sub(r"[^a-z0-9]", "", str(k).lower()): k for k in frame}
    for n in names.split():
        k = mapping.get(re.sub(r"[^a-z0-9]", "", n.lower()))
        if k is not None:
            return frame[k]
    if required:
        raise ValueError(f"Missing {names} among {list(frame)[:90]}")
    return pd.Series("", index=frame.index)


def _eq(frame):
    a = _column(frame, "RAJ2000 RAICRS RAdeg RA_obj RA _RA ra_center right_ascension rightascension")
    b = _column(frame, "DEJ2000 DEICRS DEdeg DECdeg DEC_obj DEC _DE dec_center declination")
    ra = pd.to_numeric(a, errors="coerce")
    dec = pd.to_numeric(b, errors="coerce")
    for i in frame.index[(ra.isna() | dec.isna()) & a.notna() & b.notna()]:
        p = SkyCoord(str(a[i]).strip(), str(b[i]).strip(), unit=(u.hourangle, u.deg))
        ra.loc[i], dec.loc[i] = p.ra.deg, p.dec.deg
    return ra, dec


def _tap(cat, query, endpoint=TAP):
    arguments = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=100000,
                     QUERY=query)
    raw = _download(cat, endpoint, arguments, timeout=100)
    if raw.lstrip().startswith(b"<"):
        raise ValueError(f"TAP error: {raw[:400]!r}")
    return pd.read_csv(io.BytesIO(raw), dtype=str)


def _lens_rows(cat, frame):
    if cat == "GraL X":
        frame = pd.DataFrame({"name": _column(frame, "Name").str.strip(),
                              "ra": _column(frame, "RAcdeg"),
                              "dec": _column(frame, "DEcdeg")}).drop_duplicates("name")
        if len(frame) != 364:
            raise ValueError(f"Expected 364 distinct named GraL X lenses; found {len(frame)}")
        frame["classification"] = "spectroscopically confirmed"
    else:
        ra, dec = _eq(frame)
        names = _column(frame, "Name design_id Unique_Identifier Euclid_Designation "
                        "designation EUCLID id_str object_id IntID Source GaiaDR3 GraL ID recno",
                        required=False).astype(str).str.strip()
        classification = _column(frame, "Class classification IdentClass Grade expert_grade "
                                   "Final_Grade Rank grading Type", required=False)
        frame = pd.DataFrame(dict(name=names, ra=ra, dec=dec, classification=classification))
        missing = frame.name.isin(["", "nan"])
        frame.loc[missing, "name"] = [f"J{a:.6f}{b:+.6f}" for a, b in
                                     zip(frame.loc[missing, "ra"], frame.loc[missing, "dec"])]
    frame["catalog"] = cat
    frame["benchmark"] = cat == "GraL X"
    frame["ra"] = pd.to_numeric(frame.ra)
    frame["dec"] = pd.to_numeric(frame.dec)
    if not len(frame) or not (frame.ra.between(0, 360) & frame.dec.between(-90, 90)).all():
        raise ValueError("Empty table or invalid celestial coordinates")
    return frame[LENS_COLS].drop_duplicates()


def _table_from_html(cat, raw):
    tables = BeautifulSoup(raw, "html.parser").find_all("table")
    all_rows = []
    for table in tables:
        for tr in table.find_all("tr"):
            cells = [c.get_text(" ", strip=True) for c in tr.find_all(["td", "th"], recursive=False)]
            if len(cells) < 4:
                continue
            if cat == "CASTLES":
                # CASTLES uses hh:mm:ss and signed dd:mm:ss and has no fixed leading
                # thumbnail column across page versions.
                idx = next((i for i in range(len(cells)-1) if
                            re.fullmatch(r"\d{1,2}:\d\d:[\d.]+", cells[i]) and
                            re.fullmatch(r"[+−-]?\d{1,2}:\d\d:[\d.]+", cells[i+1])), None)
                if idx is None:
                    continue
                name = next((s for s in cells[:idx] if re.search(r"[A-Za-z].*\d", s)), "")
                all_rows.append([name, cells[idx], cells[idx+1].replace("−", "-"), "CASTLES entry"])
            else:
                try:
                    ai = next(i for i in range(1, min(5, len(cells)-1)) if
                              re.search(r"\d", cells[i-1]) and
                              re.fullmatch(r"\d+\.\d+", cells[i]))
                    all_rows.append([cells[ai-1], float(cells[ai]), float(cells[ai+1]), cat])
                except (StopIteration, ValueError):
                    continue
    return pd.DataFrame(all_rows, columns=["Name", "RA", "Dec", "Class"])


def _gral_names(names):
    rows = []
    for s in sorted(set(names)):
        # GraL also uses packed coordinate designations without the decimal point.
        parsed = re.sub(r"(\d{6})(\d{3})([+-]\d{6})(\d{3})", r"\1.\2\3.\4", s)
        m = re.search(r"(\d{2})(\d{2})(\d{2}(?:\.\d+)?)([+-])(\d{2})(\d{2})(\d{2}(?:\.\d+)?)", parsed)
        if m:
            h, mi, sec, sign, d, dm, ds = m.groups()
            ra = (float(h) + float(mi)/60 + float(sec)/3600)*15
            dec = (1 if sign == "+" else -1)*(float(d) + float(dm)/60 + float(ds)/3600)
            rows.append([s, ra, dec])
    return pd.DataFrame(rows, columns=["Name", "RA", "Dec"])


def _static_cds(cat, table):
    """Read a byte-described CDS file from the official Harvard mirror.

    Preserve fixed-width identifier columns and construct decimal coordinates from
    the published sexagesimal parts only when that is the table's native encoding.
    """
    from astropy.table import Table
    prefix, filename = table.rsplit("/", 1)
    errors = []
    for host in ["https://vizier.cfa.harvard.edu/ftp/cats/",
                 "https://cdsarc.cds.unistra.fr/ftp/"]:
        try:
            root = host + prefix + "/"
            readme = _download(cat, root+"ReadMe", timeout=20)
            with tempfile.TemporaryDirectory() as tmp:
                dictionary = Path(tmp)/"ReadMe"
                dictionary.write_bytes(readme)
                file = Path(tmp)/(filename+".dat")
                try:
                    data = _download(cat, root+file.name, timeout=30)
                    if data[:100].lower().lstrip().startswith((b"<",b"<!")):
                        raise ValueError("Received HTML instead of CDS data")
                except (requests.RequestException, ValueError):
                    data = gzip.decompress(_download(cat, root+file.name+".gz", timeout=40))
                file.write_bytes(data)
                tab = Table.read(file, format="ascii.cds", readme=str(dictionary))
                df = tab.to_pandas()
            if {"RAh","RAm","RAs","DE-","DEd","DEm","DEs"}.issubset(df):
                df["RAJ2000"] = (df.RAh + df.RAm/60 + df.RAs/3600)*15
                df["DEJ2000"] = (df.DEd + df.DEm/60 + df.DEs/3600) * np.where(df["DE-"]=="-", -1, 1)
            elif cat == "DESI Dawes":
                coordinates = df.Name.astype(str).str.extract(r"(\d{1,3}\.\d+)([+-]\d{1,2}\.\d+)")
                df["RAJ2000"], df["DEJ2000"] = coordinates[0], coordinates[1]
            if cat != "Gaia FPR 381":
                _lens_rows(cat, df)
            return df
        except Exception as exc:
            errors.append(str(exc))
    raise ValueError("CDS static file unavailable: "+" / ".join(errors))


def _load_one(cat, url):
    if url.startswith("J/"):
        try:
            frame = _static_cds(cat, url)
        except Exception:
            try:
                frame = _tap(cat, f'SELECT * FROM "{url}"')
            except Exception:
                frame = None
        if frame is None:
            # Independent CfA VizieR mirror; request synthetic decimal coordinates
            # alongside source fields, avoiding dependence on TAP upload support.
            from astropy.table import Table
            raw = _download(cat, "https://vizier.cfa.harvard.edu/viz-bin/votable",
                            {"-source":url, "-out.all":1, "-out.add":"_RAJ2000,_DEJ2000",
                             "-out.max":"unlimited"})
            frame = Table.read(io.BytesIO(raw), format="votable", use_names_over_ids=True).to_pandas()
            for col in frame.select_dtypes("object"):
                frame[col] = frame[col].map(lambda x: x.decode().strip() if isinstance(x, bytes) else x)
        if cat == "Gaia FPR 381":
            frame = _gral_names(_column(frame, "Name").str.strip())
            if len(frame) != 381:
                raise ValueError(f"Expected 381 FPR leads, found {len(frame)}")
    else:
        raw = _download(cat, url)
        if cat.startswith("GLQ") or cat == "CASTLES":
            frame = _table_from_html(cat, raw)
        elif cat == "DESI DR1 strong lens VAC":
            publisher = _download(cat, url.rsplit("/",1)[0] + "/dr1_vac_dr1_strong-lensing_v1.0.sha256sum")
            if hashlib.sha256(raw).hexdigest().encode() not in publisher:
                raise ValueError("DESI release SHA256 mismatch")
            with fits.open(io.BytesIO(raw)) as hdus:
                tab = hdus["DESI_SL_CATALOG"].data
                frame = pd.DataFrame(dict(Name=tab["TARGETID"].astype(str), RA=tab["TARGET_RA"].astype(float),
                                          Dec=tab["TARGET_DEC"].astype(float), Class="DESI strong-lens target"))
        elif cat == "DESI Hsu spectroscopic pairs":
            rec = json.loads(raw)
            frames = []
            for f in rec["files"]:
                filename = f["key"].lower()
                if f["size"] > 25_000_000 or filename.endswith((".pdf", ".jpg", ".png")):
                    continue
                contents = _download(cat, f["links"]["self"])
                members = [(filename, contents)]
                if zipfile.is_zipfile(io.BytesIO(contents)):
                    with zipfile.ZipFile(io.BytesIO(contents)) as archive:
                        members = [(item.filename.lower(), archive.read(item)) for item in archive.infolist()
                                   if not item.is_dir() and item.file_size < 25_000_000 and
                                   not item.filename.startswith("__MACOSX")]
                for entry_name, contents in members:
                    if contents.startswith(b"\x1f\x8b"):
                        contents = gzip.decompress(contents)
                    ids = re.findall(rb"DESI[- ](\d{1,3}\.\d{4})([+-]\d{1,2}\.\d{4})", contents)
                    if len(ids) > 3:
                        frames.append(pd.DataFrame([(f"DESI-{a.decode()}{b.decode()}",float(a),float(b))
                                                   for a,b in set(ids)], columns=["Name","RA","DEC"]))
                if len(frames) and sum(len(t) for t in frames) >= 2364:
                    break
                try:
                    if filename.endswith(".csv"):
                        source = pd.read_csv(io.BytesIO(contents), comment="#")
                    else:
                        from astropy.table import Table
                        source = Table.read(io.BytesIO(contents), format="fits" if filename.endswith(".fits") else "ascii.cds").to_pandas()
                except Exception:
                    continue
                try:
                    names = _column(source, "Name Designation DESI System", required=False)
                    ra = _column(source, "RA_d RA_lens lensra radeg RA RAJ2000")
                    dec = _column(source, "DEC_d DEC_lens lensdec DEdeg DECdeg DEC DEJ2000")
                    part = pd.DataFrame(dict(Name=names, RA=ra, DEC=dec))
                except ValueError:
                    # Hsu's DESI names encode decimal RA and signed Dec.
                    ids = re.findall(r"DESI[- ](\d{1,3}\.\d{4})([+-]\d{1,2}\.\d{4})",
                                     source.to_csv(index=False))
                    part = pd.DataFrame([(f"DESI-{a}{b}", float(a), float(b)) for a, b in set(ids)],
                                        columns=["Name","RA","DEC"])
                frames.append(part)
            if not frames:
                raise ValueError("Published record did not contain a machine-readable lens list")
            frame = pd.concat(frames, ignore_index=True).drop_duplicates()
        elif cat == "DESI Jing quasar pairs":
            meta = json.loads(raw)
            commit = json.loads(_download(cat, f"{url}/commits/{meta['default_branch']}"))
            tree = json.loads(_download(cat, commit["commit"]["tree"]["url"], {"recursive":"1"}))
            rows = []
            for entry in tree["tree"]:
                filename = entry["path"]
                if not filename.lower().endswith((".csv", ".fits", ".fits.gz")) or entry["size"] > 35_000_000:
                    continue
                contents = _download(cat, f"https://raw.githubusercontent.com/astroliang/DESI_DR1_QP/{commit['sha']}/{filename}")
                if filename.lower().endswith(".csv"):
                    source = pd.read_csv(io.BytesIO(contents), dtype=str)
                else:
                    from astropy.table import Table
                    source = Table.read(io.BytesIO(contents), format="fits").to_pandas()
                if "RA_A" not in source:
                    continue
                for side in ("A", "B"):
                    rows.append(pd.DataFrame(dict(Name=_column(source, "System SysCenterName"),
                                RA=source["RA_"+side], Dec=source["DEC_"+side],
                                Class=_column(source, "Classify Class classification", required=False))))
            if not rows:
                raise ValueError("No author table with RA_A/DEC_A and RA_B/DEC_B found")
            frame = pd.concat(rows, ignore_index=True)
        elif cat in {"GraL IX appendix", "GraL III shortlist", "GraL VIII radio census"}:
            if not raw.startswith(b"%PDF"):
                raise ValueError("Response is not a PDF")
            pages = ["-f", "12", "-l", "13"] if cat == "GraL IX appendix" else []
            data = subprocess.run(["pdftotext", "-layout", *pages, "-", "-"], input=raw,
                                  check=True, capture_output=True).stdout.decode().replace("−", "-")
            if cat == "GraL IX appendix":
                names = re.findall(r"DR3Gaia\d{6}\.\d{3}[+-]\d{6}\.\d{2}", data)
            else:
                names = re.findall(r"(?<!\d)(?:\d{9}[+-]\d{9}|\d{6}\.\d+[+-]\d{6}(?:\.\d+)?)(?!\d)", data)
            frame = _gral_names(names)
        else:
            if cat == "Euclid Q1" and hashlib.md5(raw).hexdigest() != "21f9540189a8ea11a9165a80c341858b":
                raise ValueError("Zenodo publisher MD5 mismatch")
            if cat == "Euclid DR1 ECLIPSE":
                if hashlib.md5(raw).hexdigest() != "56433fecbdc459046ffe9eca44bb4dcf":
                    raise ValueError("Zenodo publisher MD5 mismatch")
                with zipfile.ZipFile(io.BytesIO(raw)) as z:
                    item = next(n for n in z.namelist() if n.endswith("ECLIPSE_object_catalogue.csv")
                                and not n.startswith("__MACOSX"))
                    raw = z.read(item)
            frame = pd.read_csv(io.BytesIO(raw), low_memory=False)
    minrows = {"GLQ doubles": 160, "GLQ quads": 55, "GLQ binaries": 18,
               "GLQ NIQs": 1, "CASTLES": 90, "GraL IX appendix": 48,
               "Euclid Q1": 497, "Euclid DR1 ECLIPSE": 26662, "Lemon IV": 179,
               "DESI Hsu spectroscopic pairs": 2364, "Lenscat": 32000, "DESI Jing quasar pairs": 2000,
               "GraL III shortlist": 32, "GraL VIII radio census": 10,
               "DESI He": 971, "DESI Dawes": 875, "DULAG": 5286}
    if len(frame) < minrows.get(cat, 1):
        raise ValueError(f"Too few parsed input rows: {len(frame)}")
    return _lens_rows(cat, frame)


def _write_audit(out):
    path = Path(out) / "catalog_input_checksums.csv"
    old = pd.read_csv(path, keep_default_na=False) if path.exists() and path.stat().st_size > 1 else pd.DataFrame()
    combined = pd.concat([old, pd.DataFrame(AUDIT)], ignore_index=True).fillna("")
    combined.drop_duplicates().to_csv(path, index=False)


def load_lenses(outdir):
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    frames, status = [], []
    with ThreadPoolExecutor(max_workers=4) as pool:
        submitted = [(cat, pool.submit(_load_one, cat, path)) for cat, path in LENS_URL.items()]
        for cat, fut in submitted:
            try:
                df = fut.result()
                frames.append(df)
                status.append([cat, True, len(df), ""])
            except Exception as exc:
                status.append([cat, False, 0, str(exc)[:500]])
            print("Publication catalogue", *status[-1], flush=True)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=LENS_COLS)
    combined.to_csv(out/"known_lenses.csv", index=False)
    pd.DataFrame(status, columns=["catalog", "ok", "rows", "error"]).to_csv(out/"catalog_status.csv", index=False)
    _write_audit(out)
    return combined


def _milliquas(positions):
    urls = [
        "https://vizier.cfa.harvard.edu/ftp/cats/VII/294/catalog.dat.gz",
        "https://cdsarc.cds.unistra.fr/ftp/VII/294/catalog.dat.gz",
    ]
    for url in urls:
        try:
            raw = _download(PUBLIC[0], url, timeout=65)
            # Native CDS VII/294 byte offsets. All 1,021,800 rows are required.
            columns = ["era","edec","name","classification","z","redshift_ref"]
            chunks = pd.read_fwf(gzip.GzipFile(fileobj=io.BytesIO(raw)), header=None,
                                 colspecs=[(0,11),(12,23),(25,50),(51,55),(76,82),(90,96)],
                                 names=columns, chunksize=50000, encoding="ascii")
            loc = SkyCoord(positions.ra.to_numpy()*u.deg, positions.dec.to_numpy()*u.deg)
            total, rows = 0, []
            for chunk in chunks:
                cat = SkyCoord(chunk.era.to_numpy()*u.deg, chunk.edec.to_numpy()*u.deg)
                i, j, _, _ = search_around_sky(loc, cat, 2*u.arcsec)
                if len(i):
                    part = chunk.iloc[j].copy()
                    part["row_id"] = positions.row_id.to_numpy()[i]
                    part["catalog"] = PUBLIC[0]
                    rows.append(part)
                total += len(chunk)
            if total != 1021800:
                raise ValueError(f"Milliquas row count {total} != 1021800")
            return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=PUBLIC_COLS)
        except (ValueError,OSError,requests.RequestException):
            continue
    return _public_remote(PUBLIC[0], positions)


def _desi_noirlab(positions):
    rows = []
    for pos in range(0, len(positions), 35):
        chunk = positions.iloc[pos:pos+35]
        entries = ",".join(f"({int(r.row_id)},{r.ra:.10f},{r.dec:.10f})" for r in chunk.itertuples())
        query = (f"SELECT p.row_id, z.targetid as name, z.target_ra as era, z.target_dec as edec, "
                 f"z.spectype as classification, z.z, z.zerr as z_err, z.zwarn as warning "
                 f"FROM (VALUES {entries}) AS p(row_id,ra,dec) JOIN desi_dr1.zpix z "
                 f"ON q3c_radial_query(z.target_ra,z.target_dec,p.ra,p.dec,2.0/3600.0)")
        # Public, fixed anonymous identifier from dl.authClient.ANON_TOKEN.
        # https://datalab.noirlab.edu/docs/api/_modules/dl/authClient.html
        raw = _download(PUBLIC[1], "https://datalab.noirlab.edu/query/query",
                        {"sql":query, "ofmt":"csv", "async":"False", "out":"", "drop":"False",
                         "profile":"default"}, timeout=95,
                        headers={"X-DL-AuthToken":"anonymous.0.0.anon_access",
                                 "X-DL-TimeoutRequest":"90", "Content-Type":"text/ascii"})
        group = pd.read_csv(io.BytesIO(raw), dtype=str)
        if not {"row_id","name","era","edec","classification","z","z_err","warning"}.issubset(group):
            raise ValueError(f"DESI NOIRLab query error {raw[:300]!r}")
        group["catalog"] = PUBLIC[1]
        rows.append(group)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame(columns=PUBLIC_COLS)


def _desi_vizier(positions):
    """Crossmatch at the independent CfA mirror; verify a real control in each batch."""
    url = "https://vizier.cfa.harvard.edu/viz-bin/votable"
    columns = "TargetID,RAICRS,DEICRS,OType,z,e_z,ZWARN"
    source = "V/161/zcatdr1"

    def parse(raw, allow_limit=False):
        root = ElementTree.fromstring(raw)
        for info in root.findall(".//{*}INFO"):
            value = str(info.get("value", "")).upper()
            if value == "ERROR" or (value == "OVERFLOW" and not allow_limit) or str(info.get("name", "")).lower() == "error":
                raise ValueError(f"DESI CfA VizieR: {value} {info.text}")
        # Some other VizieR responses have only metadata or descriptions.
        for table in root.findall(".//{*}TABLE"):
            if table.get("name", "") != source:
                continue
            names = [field.get("name") for field in table.findall("{*}FIELD")]
            if not set(columns.split(",")).issubset(names):
                raise ValueError(f"DESI CfA VizieR missing fields: {names}")
            data = table.find("{*}DATA/{*}TABLEDATA")
            if data is None:
                raise ValueError("DESI CfA VizieR response omitted TABLEDATA")
            records = [[v.text or "" for v in row.findall("{*}TD")] for row in data.findall("{*}TR")]
            return pd.DataFrame(records, columns=names)
        raise ValueError("DESI CfA VizieR response omitted zcatdr1")

    base = {"-source":source, "-out":columns, "-out.max":"2", "-out.form":"VOTable"}
    reference = parse(_download(PUBLIC[1],url,base,timeout=45), allow_limit=True)
    if reference.empty:
        raise ValueError("Missing DESI CfA validation source")
    seed = reference.iloc[0]
    anchor = f"{float(seed.RAICRS):.9f}{float(seed.DEICRS):+.9f}"
    frames = []
    for start in range(0,len(positions),120):
        chunk = positions.iloc[start:start+120]
        sky = [f"{r.ra:.9f}{r.dec:+.9f}" for r in chunk.itertuples()] + [anchor]
        # The VizieR CGI position-list request used by astroquery.vizier.
        # A real final control position prevents a false empty result.
        directives = [f"-source={source}",f"-out={columns}","-out.add=_q","-out.max=10000",
                      "-c.rs=2","-oc.form=d","-c=<<====Positions", *sky,"====Positions"]
        raw = _download(PUBLIC[1],url,"\n".join(directives),"POST",timeout=75)
        rows = parse(raw)
        query_ids = pd.to_numeric(rows["_q"], errors="raise").astype(int)
        if not query_ids.between(1,len(sky)).all() or len(rows)>=10000:
            raise ValueError("DESI CfA VizieR query identifiers or row limit invalid")
        control = rows[query_ids==len(sky)]
        if str(seed.TargetID).strip() not in set(control.TargetID.str.strip()):
            raise ValueError("DESI CfA VizieR did not recover its positive control")
        data = rows[query_ids < len(sky)].copy()
        data["row_id"] = chunk.row_id.to_numpy()[query_ids[query_ids < len(sky)].to_numpy()-1]
        data["catalog"] = PUBLIC[1]
        frames.append(data.rename(columns=dict(TargetID="name",RAICRS="era",DEICRS="edec",
                                               OType="classification",e_z="z_err",ZWARN="warning")))
    return pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=PUBLIC_COLS)


def _desi(positions):
    errors = []
    for label, query in [("NOIRLab zpix",_desi_noirlab), ("CfA VizieR",_desi_vizier),
                         ("Strasbourg VizieR",lambda pos: _public_remote(PUBLIC[1],pos))]:
        try:
            return query(positions)
        except Exception as exc:
            detail = f"{type(exc).__name__}: {str(exc).split('?')[0][:300]}"
            errors.append(f"{label} {detail}")
            AUDIT.append(dict(catalog=PUBLIC[1], service=label, retrieval_error=detail))
            print("DESI redshift service failed",label,detail,flush=True)
    raise ValueError(" | ".join(errors))


def _public_remote(cat, positions):
    """TAP upload crossmatch, no TOP and explicit overflow detection via row cap."""
    if cat == PUBLIC[0]:
        table = '"VII/294/catalog"'
        columns = ['"RAJ2000"', '"DEJ2000"', '"Name"', '"Type"', '"z"', '"rz"']
        aliases = ["era", "edec", "name", "classification", "z", "redshift_ref"]
    else:
        table = '"V/161/zcatdr1"'
        columns = ['"RAICRS"', '"DEICRS"', '"TargetID"', '"OType"', '"z"', '"e_z"', '"ZWARN"']
        aliases = ["era", "edec", "name", "classification", "z", "z_err", "warning"]
    from astropy.table import Table
    results = []
    for start in range(0, len(positions), 1500):
        part = positions.iloc[start:start+1500]
        payload = io.BytesIO()
        Table.from_pandas(part[["row_id", "ra", "dec"]]).write(payload, format="votable")
        # VizieR's TAP service supports upload only via multipart/VOTable.
        fields = ", ".join(f"c.{c} AS {a}" for c, a in zip(columns, aliases))
        query = (f"SELECT p.row_id, {fields} FROM TAP_UPLOAD.positions AS p JOIN {table} AS c "
                 f"ON 1=CONTAINS(POINT('ICRS',c.{columns[0]},c.{columns[1]}),"
                 "CIRCLE('ICRS',p.ra,p.dec,0.000555555556))")
        args = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=100000,
                    QUERY=query, UPLOAD="positions,param:positions")
        try:
            raw = _download(cat, TAP, args, "POST",
                            {"positions": ("positions.xml", payload.getvalue(), "application/x-votable+xml")},
                            timeout=110)
            df = pd.read_csv(io.BytesIO(raw), dtype=str)
            df.columns = df.columns.str.lower()
            if not set(["row_id"]+aliases).issubset(df) or len(df) >= 100000:
                raise ValueError("TAP schema or overflow error")
        except Exception:
            # Independent CDS XMatch server. Match distance is recomputed below.
            args2 = dict(request="xmatch", distMaxArcsec=2, selection="all", responseFormat="csv",
                         cat2="vizier:"+table.strip('"'), colRA1="input_ra", colDec1="input_dec")
            csv = part[["row_id", "ra", "dec"]].rename(columns={"ra":"input_ra","dec":"input_dec"}).to_csv(index=False)
            raw = _download(cat, "https://cdsxmatch.u-strasbg.fr/xmatch/api/v1/sync", args2, "POST",
                            {"cat1": ("positions.csv", csv.encode(), "text/csv")}, timeout=110)
            matched = pd.read_csv(io.BytesIO(raw), dtype=str)
            df = pd.DataFrame({"row_id":_column(matched, "row_id")})
            for c, a in zip(columns, aliases):
                df[a] = _column(matched, c.strip('"'))
        results.append(df)
    if not results:
        return pd.DataFrame(columns=PUBLIC_COLS)
    df = pd.concat(results, ignore_index=True)
    df["catalog"] = cat
    return df


def _public_sdss(positions):
    # DR19 SDSS-V BOSS daily coadd catalogue. Legacy SDSS I-IV are additionally
    # represented in Milliquas; that overlap is not a complete DR17 spAll scan.
    name = "spAll-lite-v6_1_3.fits.gz"
    path = "/sas/dr19/spectro/boss/redux/v6_1_3/" + name
    results = []
    target = SkyCoord(positions.ra.to_numpy()*u.deg, positions.dec.to_numpy()*u.deg)
    with tempfile.TemporaryDirectory() as tmp:
        archive = Path(tmp)/name
        for host in ["https://dr19.sdss.org", "https://data.sdss.org", "https://dev-mirror.sdss.org"]:
            try:
                checksum = hashlib.sha256()
                with requests.get(host+path, headers=HEADERS, stream=True, timeout=(12, 90)) as r:
                    r.raise_for_status()
                    with archive.open("wb") as dest:
                        for buf in r.iter_content(1024*1024):
                            checksum.update(buf)
                            dest.write(buf)
                    with gzip.open(archive, "rb") as f:
                        if f.read(6) != b"SIMPLE":
                            raise ValueError("Response is not a gzipped FITS file")
                    AUDIT.append(dict(catalog=PUBLIC[2], url=r.url, bytes=archive.stat().st_size,
                                      sha256=checksum.hexdigest()))
                break
            except (requests.RequestException, OSError, ValueError):
                if host == "https://dev-mirror.sdss.org":
                    raise
        file = Path(tmp)/"spectra.fits"
        with gzip.open(archive) as src, file.open("wb") as dest:
            shutil.copyfileobj(src, dest, length=1024*1024)
        with fits.open(file, memmap=True) as hdus:
            data = hdus[1].data
            if len(data) < 100000 or not {"ZWARNING", "Z", "CLASS"}.issubset(data.names):
                raise ValueError("Unexpected DR19 spAll-lite table")
            columns = set(data.names)
            # Catalog coordinates locate the science target. Prefer measured fiber
            # pointing where available for associating the spectrum.
            for start in range(0, len(data), 40000):
                chunk = data[start:start+40000]
                ra = np.asarray(chunk["RACAT" if "RACAT" in columns else "PLUG_RA"], dtype=float)
                dec = np.asarray(chunk["DECCAT" if "DECCAT" in columns else "PLUG_DEC"], dtype=float)
                if {"FIBER_RA", "FIBER_DEC"}.issubset(columns):
                    fr, fd = np.asarray(chunk["FIBER_RA"]), np.asarray(chunk["FIBER_DEC"])
                    valid = (fr>=0) & (fr<360) & (np.abs(fd)<=90)
                    ra, dec = np.where(valid, fr, ra), np.where(valid, fd, dec)
                valid = np.flatnonzero((ra>=0) & (ra<360) & (np.abs(dec)<=90))
                if not len(valid):
                    continue
                i, j, _, _ = search_around_sky(target, SkyCoord(ra[valid]*u.deg, dec[valid]*u.deg), 2*u.arcsec)
                if not len(i):
                    continue
                j = valid[j]
                rows = chunk[j]
                identifier = rows["CATALOGID"] if "CATALOGID" in columns else rows["SDSS_ID"]
                results.append(pd.DataFrame(dict(
                    row_id=positions.row_id.to_numpy()[i], catalog=PUBLIC[2],
                    name=identifier.astype(str), classification=np.char.strip(rows["CLASS"].astype(str)),
                    era=ra[j], edec=dec[j], z=rows["Z"].astype(float),
                    z_err=rows["Z_ERR"].astype(float), warning=rows["ZWARNING"].astype(int))))
    return pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=PUBLIC_COLS)


def _write_matches(components, lenses, out):
    registered = []
    if len(components) and len(lenses):
        c = SkyCoord(components.ra.to_numpy()*u.deg, components.dec.to_numpy()*u.deg)
        known = SkyCoord(lenses.ra.to_numpy()*u.deg, lenses.dec.to_numpy()*u.deg)
        i, j, sep, _ = search_around_sky(c, known, 10*u.arcsec)
        registered = pd.concat([components.iloc[i][["system_id", "row_id"]].reset_index(drop=True),
                                lenses.iloc[j].reset_index(drop=True)], axis=1)
        registered["distance_arcsec"] = sep.arcsec
    if not isinstance(registered, pd.DataFrame):
        registered = pd.DataFrame(columns=["system_id", "row_id"]+LENS_COLS+["distance_arcsec"])
    registered.to_csv(out/"literature_matches.csv", index=False)
    return set(registered.system_id)


def annotate(components, outdir):
    """One row per field. Close fibers support a system anchor, not two resolved spectra."""
    out = Path(outdir)
    components = components.copy().reset_index(drop=True)
    components["row_id"] = np.arange(len(components))
    if "component_id" not in components:
        components["component_id"] = (components["image_id"] if "image_id" in components else
                                      components.groupby("system_id").cumcount()+1)
    known = pd.read_csv(out/"known_lenses.csv") if (out/"known_lenses.csv").exists() else load_lenses(out)
    prior = _write_matches(components, known, out)
    status = pd.read_csv(out/"catalog_status.csv")
    status = status[~status.catalog.isin(PUBLIC)]
    loaded = status.loc[status.ok.astype(str).str.lower()=="true", "catalog"]
    lens_ok = set(LENS_URL).issubset(loaded)
    frames = []
    pub_status = []
    if len(components):
        with ThreadPoolExecutor(max_workers=3) as pool:
            tasks = [pool.submit(fn, components) for fn in [_milliquas, _desi, _public_sdss]]
            for cat, task in zip(PUBLIC, tasks):
                try:
                    frame = task.result()
                    frames.append(frame)
                    pub_status.append([cat, True, len(frame), ""])
                except Exception as exc:
                    pub_status.append([cat, False, 0, str(exc)[:500]])
                print("Public quasar catalogue", *pub_status[-1], flush=True)
    else:
        pub_status = [[cat, False, 0, "No surviving component positions to query"] for cat in PUBLIC]
    pub_ok = all(row[1] for row in pub_status)
    pd.concat([status, pd.DataFrame(pub_status, columns=status.columns)]).to_csv(out/"catalog_status.csv", index=False)
    _write_audit(out)
    matches = (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()).reindex(columns=PUBLIC_COLS)
    matches["row_id"] = pd.to_numeric(matches.row_id).astype(int)
    matches = matches.merge(components[["row_id", "component_id", "system_id", "ra", "dec"]], on="row_id")
    for col in ["era", "edec", "z", "z_err", "warning"]:
        matches[col] = pd.to_numeric(matches[col], errors="coerce")
    matches["classification"] = matches.classification.astype(str).str.strip()
    if len(matches):
        matches["distance_arcsec"] = SkyCoord(matches.era.to_numpy()*u.deg, matches.edec.to_numpy()*u.deg).separation(
            SkyCoord(matches.ra.to_numpy()*u.deg, matches.dec.to_numpy()*u.deg)).arcsec
    matches["nearest_component"] = False
    matches["isolated_fiber"] = False
    for (sysid, cat), group in matches.groupby(["system_id", "catalog"]):
        for _, obs in group.groupby(["name", "era", "edec"], dropna=False):
            nearest = obs.loc[obs.distance_arcsec.idxmin()]
            good = obs[(obs.row_id == nearest.row_id) & (obs.distance_arcsec <= 1.0)].index
            matches.loc[good, "nearest_component"] = True
            neighbors = components[(components.system_id == sysid) & (components.row_id != nearest.row_id)]
            if cat != PUBLIC[0] and len(good) and len(neighbors):
                s = SkyCoord(nearest.era*u.deg, nearest.edec*u.deg).separation(
                    SkyCoord(neighbors.ra.to_numpy()*u.deg, neighbors.dec.to_numpy()*u.deg)).arcsec.min()
                matches.loc[good, "isolated_fiber"] = s >= 2.5 and nearest.distance_arcsec <= 0.65
    matches.to_csv(out/"external_matches.csv", index=False)
    result = []
    for system in components.system_id.unique():
        group = matches[(matches.system_id==system) & matches.nearest_component]
        spectra = group[(group.catalog!=PUBLIC[0]) & (group.warning==0)]
        qso = spectra[(spectra.classification=="QSO") & (spectra.z>0)]
        resolved = qso[qso.isolated_fiber]
        n = resolved.row_id.nunique()
        redshift_differences = []
        for a, b in combinations(resolved.itertuples(), 2):
            if a.row_id == b.row_id:
                continue
            dz = abs(a.z-b.z)
            sigma = np.hypot(a.z_err, b.z_err)
            discrepant = dz > max(0.01*(1+(a.z+b.z)/2), 5*sigma if np.isfinite(sigma) else 0)
            redshift_differences.append(discrepant)
        result.append([system, system in prior, lens_ok, pub_ok, bool(len(qso) or
                       ((group.catalog==PUBLIC[0]) & group.classification.str.startswith("Q")).any()),
                       bool((spectra.classification=="STAR").any()),
                       any(redshift_differences), n, bool(n>=2 and not any(redshift_differences))])
    return pd.DataFrame(result, columns=["system_id", "prior_match", "lens_catalogues_complete",
                        "public_catalogues_complete", "qso_evidence", "star_evidence",
                        "discordant_redshifts", "resolved_qso_spectra", "compatible_redshifts"])
