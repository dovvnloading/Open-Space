"""Read the complete corrected RVBank epoch table, or the ESO Phase 3 RV catalogue."""
from runtime import require_survey_workflow

import argparse
import hashlib
import json
import http.client
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

CDS = "https://cdsarc.cds.unistra.fr/ftp"
VIZTAP = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync"
ESOTAP = "https://archive.eso.org/tap_cat/sync"
PRIMARY = "J/A+A/683/A125/table4"
SPORES = "J/AJ/170/343/table5"
FIBRE_CHANGE = 2457169.0  # Trifonov et al. Fig. 5, during the May 2015 interruption.
COLUMNS = ["star", "bjd", "rv", "err", "era", "fwhm", "bis", "crx", "dlw", "halpha",
           "sindex", "rhkp", "epoch_id"]
# CDS one-based, inclusive byte positions. No whitespace splitting of sentinel fields.
RV_FIELDS = {
    "star": ("Target", 1, 14), "ra_deg": ("RAdeg", 16, 35),
    "dec_deg": ("DEdeg", 37, 57), "bjd": ("BJD", 59, 73),
    "rv": ("DRVmlcnzp", 75, 87), "err": ("e_DRVmlcnzp", 89, 98),
    "crx": ("CRX", 206, 219), "dlw": ("dLW", 235, 246),
    "halpha": ("Halpha", 259, 267), "flag": ("Flag", 318, 322),
    "fwhm": ("FWHMDRS", 324, 335), "bis": ("BIS", 351, 371),
    "epoch_id": ("TimeId", 535, 557), "rhkp": ("RHKp", 662, 682),
}
S_FIELDS = {"hd": ("HD", 1, 8), "bjd": ("BJD", 10, 23), "sindex": ("S-ind", 44, 50),
            "inst": ("Inst", 140, 150)}


def download(url, path, manifest):
    request = Request(url, headers={"User-Agent": "Open-Space HARPS catalog research"})
    for attempt in range(1, 4):
        digest = hashlib.sha256()
        received = 0
        entry = dict(url=url, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                     file=path.name, attempt=attempt)
        try:
            print(f"Retrieving {path.name}, attempt {attempt}", flush=True)
            with urlopen(request, timeout=120) as stream, path.open("wb") as sink:
                expected = stream.headers.get("Content-Length")
                while block := stream.read(1024 * 1024):
                    sink.write(block)
                    digest.update(block)
                    received += len(block)
                if expected and received != int(expected):
                    raise OSError(f"short response: {received} of {expected} bytes")
            manifest.append(dict(**entry, status="downloaded", bytes=received,
                                 sha256=digest.hexdigest(), error=""))
            return path
        except (HTTPError, URLError, OSError, http.client.HTTPException) as exc:
            manifest.append(dict(**entry, status="failed", bytes=received,
                                 sha256="", error=str(exc)))
            path.unlink(missing_ok=True)
            if isinstance(exc, HTTPError) and exc.code not in (408, 425, 429, 500, 502, 503, 504):
                raise
            if attempt == 3:
                raise
            time.sleep(2 ** attempt)


def tap(endpoint, query, limit, file, manifest):
    params = {"request": "doQuery", "lang": "ADQL", "format": "csv",
              "maxrec": limit, "query": query}
    download(endpoint + "?" + urlencode(params), file, manifest)
    return pd.read_csv(file, low_memory=False)


def catalog(table, fields, expected, work, manifest, errors):
    folder, name = table.rsplit("/", 1)
    stem = folder.replace("/", "-") + "-" + name
    for suffix in (".dat.gz", ".dat", ".tap.csv"):
        file = work / (stem + suffix)
        try:
            if suffix == ".tap.csv":
                data = tap(VIZTAP, f'SELECT * FROM "{table}"', expected + 10, file, manifest)
                data.columns = data.columns.str.strip().str.lower()
                # VizieR may standardise CDS coordinate labels in TAP tables.
                ra = next((key for key in ("radeg", "raj2000", "ra_icrs") if key in data), "")
                dec = next((key for key in ("dedeg", "dej2000", "de_icrs") if key in data), "")
                if ra and dec:
                    data = data.rename(columns={ra: "radeg", dec: "dedeg"})
                data = data[[f[0].lower() for f in fields.values()]].copy()
                data.columns = list(fields)
            else:
                url = f"{CDS}/{folder}/{name}{suffix}"
                download(url, file, manifest)
                data = pd.read_fwf(file, names=list(fields),
                                   colspecs=[(a-1, b) for _, a, b in fields.values()],
                                   na_values=["-", "nan"], dtype=str)
            if len(data) != expected:
                raise ValueError(f"expected {expected} rows, received {len(data)}")
            return data
        except Exception as exc:
            if manifest and manifest[-1]["file"] == file.name and manifest[-1]["status"] == "downloaded":
                manifest[-1].update(status="unusable_response", error=str(exc))
            errors.append(f"{table} {suffix}: {exc}")
    raise RuntimeError(f"All retrieval routes for {table} failed")


def eso_catalog(work, manifest):
    metadata = tap(ESOTAP, "SELECT table_name FROM TAP_SCHEMA.columns "
                   "WHERE column_name='drs_ccf_rvc'", 100, work/"eso-tables.csv", manifest)
    names = metadata.iloc[:, 0].dropna().astype(str).unique()
    # Select the documented release by schema and published number of rows.
    for name in sorted(names):
        if not re.fullmatch(r"[a-zA-Z0-9_.]+", name):
            continue
        counts = tap(ESOTAP, f"SELECT COUNT(*) AS n FROM {name}", 2,
                     work/"eso-count.csv", manifest)
        if int(counts.iloc[0, 0]) != 289843:
            continue
        cols = ["adp_dpid", "main_id_simbad", "ra_simbad", "dec_simbad", "drs_bjd",
                "drs_ccf_rvc", "drs_dvrms", "drs_ccf_fwhm", "drs_bis_span",
                "ins_mode", "drs_mask"]
        data = tap(ESOTAP, f"SELECT {','.join(cols)} FROM {name}", 300000,
                   work/"eso-phase3.csv", manifest)
        if len(data) != 289843:
            raise ValueError("truncated ESO catalogue")
        data.columns = [c.strip().lower() for c in data.columns]
        data = data[cols].copy()
        data.columns = ["epoch_id", "star", "ra_deg", "dec_deg", "bjd", "rv", "err",
                        "fwhm", "bis", "mode", "mask"]
        data["rv"] = pd.to_numeric(data.rv, errors="coerce") * 1000
        return data
    raise ValueError("ESO Phase 3 HARPS DR1 table could not be identified")


def add_sindex(data, work, manifest, errors):
    """Require a unique source and target epoch within 2 seconds of each other."""
    aux = catalog(SPORES, S_FIELDS, 153902, work, manifest, errors)
    # Henry Draper numbers can be padded or have a stellar-component suffix.
    # A source HD123A must not match HD123 or HD123B.
    def hd_key(names, require_prefix):
        pattern = r"^HD0*([1-9][0-9]*)([A-Z]*)$" if require_prefix else r"^(?:HD)?0*([1-9][0-9]*)([A-Z]*)$"
        compact = names.astype("string").str.upper().str.replace(r"\s+", "", regex=True)
        parts = compact.str.extract(pattern)
        return "HD" + parts[0] + parts[1]

    aux["star"] = hd_key(aux.hd, require_prefix=False)
    aux["bjd"] = pd.to_numeric(aux.bjd, errors="coerce")
    aux["sindex"] = pd.to_numeric(aux.sindex, errors="coerce")
    aux = aux[aux.inst.str.strip().str.upper().isin(["HARPS", "HARPS-PRE", "HARPS-POST"]) &
              aux.sindex.gt(0) & aux.bjd.notna()]
    name = hd_key(data.star, require_prefix=True)
    for star, group in aux.groupby("star"):
        left = data.loc[name.eq(star) & data.bjd.notna(), "bjd"].sort_values()
        if left.empty:
            continue
        right = group.sort_values("bjd")
        times = right.bjd.to_numpy()
        limit = 2 / 86400
        lo = np.searchsorted(times, left.to_numpy()-limit)
        hi = np.searchsorted(times, left.to_numpy()+limit, side="right")
        unique = hi - lo == 1
        used = pd.Series(lo[unique]).duplicated(keep=False).to_numpy()
        data.loc[left.index[unique][~used], "sindex"] = right.sindex.to_numpy()[lo[unique][~used]]


def classify(data, source):
    data["star"] = data.star.astype("string").str.strip().replace(
        {"": pd.NA, "nan": pd.NA, "-": pd.NA, "None": pd.NA, "UNKNOWN": pd.NA})
    data["epoch_id"] = data.epoch_id.astype("string").str.strip()
    for col in ["ra_deg", "dec_deg", "bjd", "rv", "err", "fwhm", "bis", "crx", "dlw",
                "halpha", "rhkp", "flag"]:
        if col in data:
            data[col] = pd.to_numeric(data[col], errors="coerce").replace([np.inf,-np.inf],np.nan)
        else:
            data[col] = np.nan
    data["sindex"] = np.nan
    data["era"] = np.where(data.bjd < FIBRE_CHANGE, "PRE", "POST")
    data["disposition"] = "used"

    def reject(mask, reason):
        data.loc[mask & data.disposition.eq("used"), "disposition"] = reason

    reject(data.star.isna(), "missing_canonical_id")
    reject(~data.bjd.between(2452700, 2465000), "missing_or_invalid_bjd")
    reject(data.rv.isna() | data.err.isna() | data.err.le(0), "invalid_rv_or_error")
    if source == PRIMARY:
        reject(data.flag.ne(0), "serval_quality_flag")
    else:
        reject(~data["mode"].str.strip().eq("HAM"), "not_HAM")
        reject(data["mask"].isna(), "missing_ccf_mask")
        mask_counts = (data[data.disposition.eq("used")].groupby(["star","mask"]).size()
                       .rename("n").reset_index().sort_values(["n", "mask"], ascending=[False, True]))
        dominant = mask_counts.drop_duplicates("star").set_index("star")["mask"]
        reject(data["mask"].ne(data.star.map(dominant)), "minority_ccf_mask")
    # A repeated target/time pair cannot contribute extra precision.
    keep = data.loc[data.disposition.eq("used")].sort_values("err", kind="stable")
    duplicates = keep.index[keep.duplicated(["star","bjd"])]
    data.loc[duplicates, "disposition"] = "duplicate_star_bjd"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require_survey_workflow()
    work, out = args.work, args.out
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    manifest, errors = [], []
    release = PRIMARY
    try:
        try:
            data = catalog(PRIMARY, RV_FIELDS, 252615, work, manifest, errors)
        except RuntimeError:
            release = "ESO Phase 3 HARPS RV catalogue DR1"
            data = eso_catalog(work, manifest)
        classify(data, release)
        if not data.disposition.eq("used").any():
            raise ValueError("no usable primary RV measurements")
        try:
            add_sindex(data, work, manifest, errors)
        except Exception as exc:
            data["sindex"] = np.nan
            errors.append(f"S-index retrieval or matching unavailable: {exc}")
        data["n_epochs"] = data.disposition.eq("used").astype(int)
        targets = data.groupby("star", sort=True).agg(
            ra_deg=("ra_deg", "median"), dec_deg=("dec_deg", "median"),
            n_epochs=("n_epochs", "sum"), n_input=("n_epochs", "size")).reset_index()
        targets["source"] = release
        targets.to_csv(out/"targets.csv", index=False)
        compression = {"method": "gzip", "mtime": 0, "compresslevel": 6}
        data[["star", "epoch_id", "bjd", "disposition"]].to_csv(
            out/"epoch-ledger.csv.gz", index=False, compression=compression)
        used = data.loc[data.n_epochs.eq(1), COLUMNS].sort_values(["star", "bjd"])
        used.to_csv(work/"epochs.csv.gz", index=False, compression=compression)
        info = dict(source=release, input_rows=len(data), targets=len(targets),
                    used_epochs=len(used), dispositions=data.disposition.value_counts().to_dict(),
                    sindex_epochs=int(used.sindex.notna().sum()),
                    first_bjd=float(used.bjd.min()), last_bjd=float(used.bjd.max()))
    except Exception as exc:
        errors.append(str(exc))
        info = dict(source=release, status="ingest_failed")
        raise
    finally:
        pd.DataFrame(manifest, columns=["url", "retrieved_utc", "file", "attempt", "status", "bytes", "sha256", "error"]).to_csv(
            out/"inputs.csv", index=False)
        (out/"ingest.json").write_text(json.dumps(dict(**info, warnings=errors), indent=2)+"\n")
        print(json.dumps(dict(**info, warnings=errors)), flush=True)


if __name__ == "__main__":
    main()
