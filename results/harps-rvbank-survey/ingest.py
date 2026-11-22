"""Read the complete corrected RVBank epoch table, or the ESO Phase 3 RV catalogue."""
from activity_supplement import add_published_sindex
from runtime import require_survey_workflow

import argparse
import csv
import zipfile
import hashlib
import json
import http.client
from html import unescape
import re
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote, urlencode

from astropy.table import Table
from urllib.request import Request, urlopen

import numpy as np
import pandas as pd

CDS_STATIC = ("https://vizier.cfa.harvard.edu/ftp/cats",
              "https://cdsarc.cds.unistra.fr/ftp/cats",
              "https://cdsarc.cds.unistra.fr/ftp")
AUTHOR_COMMIT = "09626a652ddc93df8771ce3d818ee04a292f7585"
AUTHOR_BLOB = "bf5e0f4ceed9e8ce07e9412cc33646b45cb6c56a"
AUTHOR_URL = (f"https://raw.githubusercontent.com/3fon3fonov/HARPS_RVBank/{AUTHOR_COMMIT}/"
              "HARPS_RVBank_ver02.csv.zip")
AUTHOR = "Perdelwitz et al. author HARPS-RVBank Ver.02b, 2024-03-03"
ESO = "ESO Phase 3 HARPS RV catalogue DR1"
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


def download(url, path, manifest, timeout=60, attempts=2):
    request = Request(url, headers={"User-Agent": "Open-Space HARPS catalog research"})
    for attempt in range(1, attempts + 1):
        digest = hashlib.sha256()
        received = 0
        entry = dict(url=url, retrieved_utc=datetime.now(timezone.utc).isoformat(),
                     file=path.name, attempt=attempt)
        try:
            print(f"Retrieving {path.name}, attempt {attempt}", flush=True)
            with urlopen(request, timeout=timeout) as stream, path.open("wb") as sink:
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
            if attempt == attempts:
                raise
            time.sleep(2 ** attempt)


def tap(endpoint, query, limit, file, manifest, timeout=120, attempts=2):
    params = {"request": "doQuery", "lang": "ADQL", "format": "csv",
              "maxrec": limit, "query": query}
    download(endpoint + "?" + urlencode(params), file, manifest, timeout, attempts)
    return pd.read_csv(file, low_memory=False)


def field_projection(frame, fields):
    frame.columns = frame.columns.str.strip().str.lower()
    # Mirrors may standardize equatorial column names. sexagesimal values
    # must never be interpreted as decimal degrees.
    aliases = {"radeg": ("radeg", "_raj2000", "ra_icrs", "raj2000"),
               "dedeg": ("dedeg", "_dej2000", "de_icrs", "dej2000")}
    result = {}
    for key, (label, _, _) in fields.items():
        choices = aliases.get(label.lower(), (label.lower(),))
        found = next((col for col in choices if col in frame), None)
        if found is None:
            raise ValueError(f"Missing {label}. Available columns: {list(frame)}")
        values = frame[found]
        if key in ("ra_deg", "dec_deg") and found != label.lower():
            numeric = pd.to_numeric(values, errors="coerce")
            if values.notna().any() and not numeric.notna().any():
                raise ValueError(f"Cannot interpret {found} as degrees")
        result[key] = values
    return pd.DataFrame(result)


ASU_HOSTS = ("vizier.cfa.harvard.edu", "vizier.idia.ac.za")


def asu_recno_limit(table, work, manifest, errors, audit):
    """Ask each mirror for its highest record, independently of the paper's count."""
    results = {}
    for host in ASU_HOSTS:
        try:
            # VizieR-3 contains only catalogue metadata and the constraint form.
            endpoint = f"https://{host}/viz-bin/VizieR-3?" + urlencode({"-source": table})
            html = download(endpoint, work/f"table4.metadata.{host}.html", manifest,
                            timeout=15, attempts=1).read_text(errors="replace")
            plain = " ".join(unescape(re.sub(r"<[^>]+>", " ", html)).split())
            candidates = re.findall(re.escape(table) + r"[^()]{0,200}\(([\d, ]+) rows\)", plain, re.I)
            audit.setdefault("vizier_metadata_row_counts", {})[host] = sorted(
                {int(n.replace(",", "").replace(" ", "")) for n in candidates})
        except Exception as exc:
            errors.append(f"{host} CDS table metadata unavailable: {exc}")
        params = {"-source": table, "-out": "recno", "-out.max": 1, "-sort": "-recno",
                  "-out.meta": "h"}
        path = work / f"table4.last-record.{host}.tsv"
        try:
            download(f"https://{host}/viz-bin/asu-tsv?" + urlencode(params),
                     path, manifest, timeout=20, attempts=1)
            frame = pd.read_csv(path, comment="#", sep="\t", dtype=str)
            frame.columns = frame.columns.str.strip().str.lower()
            nums = pd.to_numeric(frame["recno"], errors="coerce").dropna()
            if len(nums) != 1 or nums.iloc[0] <= 1 or nums.iloc[0] % 1:
                raise ValueError("invalid descending record index")
            limit = int(nums.iloc[0])
            # The inclusive lower bound guarantees a nonempty response. Any
            # record above limit makes this two-record boundary check fail.
            params.update({"recno": f"{limit}..9999999", "-sort": "recno", "-out.max": 2})
            boundary = work / f"table4.boundary.{host}.tsv"
            download(f"https://{host}/viz-bin/asu-tsv?" + urlencode(params),
                     boundary, manifest, timeout=20, attempts=1)
            frame = pd.read_csv(boundary, comment="#", sep="\t", dtype=str)
            frame.columns = frame.columns.str.strip().str.lower()
            boundary_records = pd.to_numeric(frame["recno"], errors="coerce").dropna()
            if boundary_records.tolist() != [limit]:
                raise ValueError(f"ascending boundary disagrees with terminal record {limit}")
            results[host] = limit
        except Exception as exc:
            errors.append(f"{host} terminal record unavailable: {exc}")
    return results


def asu_catalog(table, fields, expected, work, manifest, errors, hosts=ASU_HOSTS, retain_prefix=False):
    """Request bounded, gap-checked pages from independent VizieR mirrors."""
    parts = []
    stem = table.replace("/", "-")
    for first in range(1, expected + 1, 25000):
        last = min(first + 24999, expected)
        prefix, prefix_receipt = None, None
        params = {"-source": table, "-out": "**", "-out.max": last - first + 2,
                  "-out.meta": "h", "-sort": "recno", "recno": f"{first}..{last}"}
        for host in hosts:
            url = f"https://{host}/viz-bin/asu-tsv?" + urlencode(params)
            path = work / f"{stem}.{host}.{first}-{last}.tsv"
            try:
                download(url, path, manifest, timeout=40, attempts=1)
                raw = pd.read_csv(path, comment="#", sep="\t", dtype=str, quoting=csv.QUOTE_NONE)
                raw.columns = raw.columns.str.strip().str.lower()
                recno = pd.to_numeric(raw["recno"], errors="coerce")
                # ASU can print unit and separator rows below its column names.
                data_rows = recno.notna()
                numbers = np.sort(recno[data_rows].to_numpy())
                np.testing.assert_array_equal(numbers, np.arange(first, first + len(numbers)))
                part = field_projection(raw[data_rows].copy(), fields)
                part["recno"] = recno[data_rows].astype(int)
                if len(numbers) != last - first + 1:
                    if 0 < len(numbers) < last-first+1:
                        if prefix is None or len(part) > len(prefix):
                            prefix, prefix_receipt = part, manifest[-1]
                    raise ValueError(f"requested {first}..{last}, received {len(numbers)} contiguous rows")
                parts.append(part)
                break
            except Exception as exc:
                if manifest and manifest[-1]["file"] == path.name and manifest[-1]["status"] == "downloaded":
                    manifest[-1].update(status="unusable_response", error=str(exc)[:500])
                errors.append(f"{table} {host} rows {first}..{last}: {str(exc)[:500]}")
        else:
            if retain_prefix and (parts or prefix is not None):
                if prefix is not None:
                    prefix_receipt["status"] = "used_verified_prefix"
                    parts.append(prefix)
                errors.append(f"Retained only the contiguous available prefix of {table}")
                break
            raise RuntimeError(f"All mirror routes for {table} rows {first}..{last} failed")
    return pd.concat(parts, ignore_index=True)


def catalog(table, fields, expected, work, manifest, errors, audit=None, subsets=None):
    """CDS fixed-width mirrors, checked ASU pages, then one bounded TAP attempt."""
    stem = table.replace("/", "-")
    for host in CDS_STATIC:
        for suffix in (".dat.gz", ".dat"):
            path = work / (stem + "." + host.split("/")[2] + suffix)
            url = f"{host}/" + quote(table+suffix, safe="/")
            try:
                download(url, path, manifest, timeout=20, attempts=1)
                data = pd.read_fwf(path, names=list(fields),
                                   colspecs=[(a-1, b) for _, a, b in fields.values()],
                                   na_values=["-", "nan"], dtype=str)
                if len(data) != expected:
                    raise ValueError(f"expected {expected} rows, received {len(data)}")
                return data
            except Exception as exc:
                if manifest[-1]["status"] == "downloaded":
                    manifest[-1].update(status="unusable_response", error=str(exc)[:500])
                errors.append(f"{url}: {str(exc)[:500]}")
    if table == PRIMARY:
        # CDS's FTP-to-FITS converter can expose the archive original even if
        # the independently indexed VizieR database contains fewer rows.
        for host in ("cdsarc.cds.unistra.fr", "vizier.cfa.harvard.edu"):
            path = work/f"table4.original.{host}.fits"
            url = f"https://{host}/viz-bin/nph-Cat/fits?" + quote(table+".dat", safe="/")
            try:
                download(url, path, manifest, timeout=25, attempts=1)
                original = Table.read(path, format="fits", hdu=1)
                if len(original) != expected:
                    raise ValueError(f"expected {expected} rows, received {len(original)}")
                data = field_projection(original.to_pandas(), fields)
                dates = pd.to_numeric(data.bjd, errors="coerce")
                if not dates.between(2452700, 2465000).any():
                    raise ValueError("no plausible absolute BJDs in original-file FITS")
                for col in ("star", "epoch_id"):
                    if data[col].map(lambda x: isinstance(x, bytes)).any():
                        data[col] = data[col].map(lambda x: x.decode("ascii") if isinstance(x, bytes) else x)
                return data
            except Exception as exc:
                if manifest[-1]["status"] == "downloaded":
                    manifest[-1].update(status="unusable_response", error=str(exc)[:500])
                errors.append(f"{url}: {str(exc)[:500]}")
    try:
        if audit is None:
            return asu_catalog(table, fields, expected, work, manifest, errors)
        limits = asu_recno_limit(table, work, manifest, errors, audit)
        audit["mirror_last_recno"] = limits
        # If two mirrors disagree, independently attempt the largest complete
        # reconstruction first. Each page must cover its exact inclusive range.
        for last in sorted(set(limits.values()) or {expected}, reverse=True):
            mirrors = tuple(h for h in ASU_HOSTS if not limits or limits.get(h) == last)
            if not 0 < last <= expected:
                errors.append(f"Unexpected maximum record {last} in {mirrors}")
                continue
            try:
                rows = asu_catalog(table, fields, last, work, manifest, errors, mirrors, retain_prefix=True)
            except RuntimeError as exc:
                errors.append(str(exc))
                continue
            if len(rows) == expected:
                return rows
            subsets.append(rows)
            audit["verified_contiguous_table_rows"] = len(rows)
            audit["verified_to_reported_endpoint"] = bool(limits and len(rows) == last)
            break
    except RuntimeError as exc:
        errors.append(str(exc))
    path = work / (stem + ".tap.csv")
    try:
        data = tap(VIZTAP, f'SELECT * FROM "{table}"', expected + 10, path, manifest,
                   timeout=30, attempts=1)
        if len(data) != expected:
            raise ValueError(f"expected {expected} rows, received {len(data)}")
        return field_projection(data, fields)
    except Exception as exc:
        if manifest[-1]["status"] == "downloaded":
            manifest[-1].update(status="unusable_response", error=str(exc)[:500])
        errors.append(f"{table} TAP: {str(exc)[:500]}")
        raise RuntimeError(f"All corrected CDS routes for {table} failed") from exc


def author_rvbank(work, manifest):
    """Parse the author-published 2024 v2b archive using its CSV header."""
    path = download(AUTHOR_URL, work/"RVBank-author-v02b.csv.zip", manifest)
    # Verify the public Git object checksum, including its Git object header.
    blob = hashlib.sha1(b"blob " + str(path.stat().st_size).encode() + b"\0")
    with path.open("rb") as raw:
        while buf := raw.read(1024 * 1024):
            blob.update(buf)
    if blob.hexdigest() != AUTHOR_BLOB:
        raise ValueError(f"Archive differs from the file in author commit {AUTHOR_COMMIT}")
    names = {"star":"target", "ra_deg":"ra", "dec_deg":"dec", "bjd":"bjd",
             "rv":"rv_mlc_nzp", "err":"e_rv_mlc_nzp", "crx":"crx", "dlw":"dlw",
             "halpha":"halpha", "flag":"f_rv", "fwhm":"fwhm_drs", "bis":"bis",
             "epoch_id":"timeid", "rhkp":"rhkp"}
    parts = []
    with zipfile.ZipFile(path) as archive:
        members = [n for n in archive.namelist() if n.lower().endswith(".csv")
                   and not any(c.startswith((".", "__MACOSX")) for c in n.split("/"))]
        for member in sorted(members):
            with archive.open(member) as stream:
                header = stream.readline().decode("utf-8-sig").strip().lstrip("#").strip()
            dialect = csv.Sniffer().sniff(header, delimiters=",;\t|")
            columns = next(csv.reader([header], dialect=dialect, skipinitialspace=True))
            cols = [c.strip().lower() for c in columns]
            if not set(names.values()).issubset(cols):
                raise ValueError(f"Unexpected author CSV header in {member}: {cols}")
            with archive.open(member) as stream:
                data = pd.read_csv(stream, names=cols, skiprows=1, sep=dialect.delimiter,
                                   quotechar=dialect.quotechar or '"', skipinitialspace=True,
                                   usecols=list(names.values()), dtype=str,
                                   na_values=["-", "None", "NULL", "null", "NaN", r"\N"],
                                   encoding="latin-1")
            parts.append(data.rename(columns={v:k for k,v in names.items()}))
    if not parts:
        raise ValueError("Author release contains no CSV tables")
    data = pd.concat(parts, ignore_index=True)
    # Full coverage is guaranteed by the pinned blob, CRC-checked members,
    # and parsing every CSV. Record the actual count against the paper's total.
    if data.empty or not pd.to_numeric(data.bjd, errors="coerce").between(2452700,2465000).any():
        raise ValueError("Author release has no rows with a valid absolute BJD")
    return data


def correction_crossmatch(author, table):
    """Match two releases of the same measurements without sky-only associations."""
    def keys(frame):
        star = frame.star.astype("string").str.strip().str.upper().str.replace(r"\s+", "", regex=True)
        star = star.replace({"": pd.NA, "-": pd.NA, "NAN": pd.NA, "UNKNOWN": pd.NA, "NONE": pd.NA})
        stamp = frame.epoch_id.astype("string").str.extract(
            r"(\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2}\.\d{3})(?!\d)", expand=False)
        stamp = stamp.str.replace(" ", "T", regex=False)
        return pd.DataFrame({"target": star, "stamp": stamp,
                             "bjd": pd.to_numeric(frame.bjd, errors="coerce")})
    a, b = keys(author), keys(table)
    ledger = table[["recno", "star", "epoch_id", "bjd"]].copy()
    ledger["author_index"] = pd.Series(pd.NA, index=ledger.index, dtype="Int64")
    ledger["match"] = "no_unique_identity_and_epoch"
    # A known time ID on both sides must agree. For the fallback permit at
    # most 0.022 s BJD separation and require reciprocal uniqueness.
    groups = a.groupby("target", sort=False).groups
    for target, right in b.groupby("target", sort=False):
        if target not in groups:
            continue
        left = a.loc[groups[target]].dropna(subset=["bjd"]).sort_values("bjd")
        right = right.dropna(subset=["bjd"])
        if left.empty or right.empty:
            continue
        unique_right = right.stamp.notna() & ~right.stamp.duplicated(keep=False)
        stamp_to_row = pd.Series(left.index, index=left.stamp)
        stamp_to_row = stamp_to_row[stamp_to_row.index.notna() & ~stamp_to_row.index.duplicated(keep=False)]
        exact = right.loc[unique_right, "stamp"].map(stamp_to_row).dropna().astype(int)
        if len(exact):
            bjd_ok = np.abs(left.bjd.reindex(exact.to_numpy()).to_numpy()
                            - right.bjd.loc[exact.index].to_numpy()) <= 2/86400
            ids = exact.index[bjd_ok]
            ledger.loc[ids, "author_index"] = exact.loc[ids]
            ledger.loc[ids, "match"] = "target_timestamp_and_bjd"
        pending = right.loc[ledger.loc[right.index, "author_index"].isna()]
        if pending.empty:
            continue
        times = left.bjd.to_numpy()
        lo = np.searchsorted(times, pending.bjd.to_numpy() - 2.5e-7)
        hi = np.searchsorted(times, pending.bjd.to_numpy() + 2.5e-7, side="right")
        one = hi - lo == 1
        ids = pending.index[one]
        suggested = left.iloc[lo[one]]
        old_stamp = suggested.stamp.reset_index(drop=True)
        new_stamp = pending.loc[ids, "stamp"].reset_index(drop=True)
        agrees = (old_stamp.isna() | new_stamp.isna() | old_stamp.eq(new_stamp)).fillna(False).to_numpy(dtype=bool)
        ledger.loc[ids[agrees], "author_index"] = suggested.index[agrees].to_numpy()
        ledger.loc[ids[agrees], "match"] = "unique_target_bjd"
    mapped = ledger.author_index.dropna().astype(int)
    if len(mapped):
        a_pos = author.loc[mapped.to_numpy(), ["ra_deg", "dec_deg"]].apply(pd.to_numeric, errors="coerce")
        b_pos = table.loc[mapped.index, ["ra_deg", "dec_deg"]].apply(pd.to_numeric, errors="coerce")
        ra1, dec1 = np.deg2rad(a_pos.to_numpy()).T
        ra2, dec2 = np.deg2rad(b_pos.to_numpy()).T
        h = np.sin((dec2-dec1)/2)**2 + np.cos(dec1)*np.cos(dec2)*np.sin((ra2-ra1)/2)**2
        mismatch = np.rad2deg(2*np.arcsin(np.sqrt(h.clip(0,1)))) > 2/3600
        ids = mapped.index[mismatch]
        ledger.loc[ids, "match"] = "coordinate_conflict"
        ledger.loc[ids, "author_index"] = pd.NA
    collision = ledger.author_index.notna() & ledger.author_index.duplicated(keep=False)
    ledger.loc[collision, "match"] = "ambiguous_author_exposure"
    ledger.loc[collision, "author_index"] = pd.NA
    return ledger


def apply_partial_table(author, table, audit, out):
    ledger = correction_crossmatch(author, table)
    matched = ledger.dropna(subset=["author_index"])
    dest = matched.author_index.to_numpy(dtype=int)
    table_rows = matched.index
    old = author.loc[dest].reset_index(drop=True)
    new = table.loc[table_rows].reset_index(drop=True)
    gzip = {"method": "gzip", "mtime": 0}
    ledger.to_csv(out/"correction-crossmatch.csv.gz", index=False, compression=gzip)
    audit.update(author_rows=len(author), matched_unique_epochs=len(matched),
                 author_epochs_without_match=len(author)-len(matched),
                 unmatched_table_epochs=len(table)-len(matched),
                 crossmatch_counts=ledger.match.value_counts().to_dict(),
                 fields={})
    changes = []
    # Half-unit rounding allowance for fixed decimal columns. E13.6 columns
    # can differ by up to 5e-7 in relative terms from a full-precision export.
    absolute = {"rv": 0, "err": 5.01e-7, "fwhm": 5.01e-7, "bis": 0,
                "halpha": 5.01e-7, "crx": 0, "dlw": 0, "rhkp": 1e-15,
                "flag": 0, "ra_deg": 1e-10, "dec_deg": 1e-10, "bjd": 5.01e-8}
    groups = [old.star.astype("string").str.strip(),
              pd.to_numeric(old.bjd, errors="coerce").ge(FIBRE_CHANGE)]
    for col, atol in absolute.items():
        x = pd.to_numeric(old[col], errors="coerce").replace([np.inf, -np.inf], np.nan)
        y = pd.to_numeric(new[col], errors="coerce").replace([np.inf, -np.inf], np.nan)
        both = x.notna() & y.notna()
        differ = ~(x.eq(y) | (x.isna() & y.isna()))
        rtol = 0 if col in ("flag", "ra_deg", "dec_deg", "bjd") else 5.01e-7
        agree = pd.Series(np.isclose(x, y, atol=atol, rtol=rtol), index=x.index) & both
        fraction = float(agree.sum()/both.sum()) if both.any() else None
        # Match scientific scales empirically across repeated exposures, then
        # demand confirmation within the same star and instrumental zero point.
        global_ok = bool(both.sum() >= 20 and fraction >= .95
                         and (col == "flag" or x[agree].nunique() >= 3))
        n_local = both.astype(int).groupby(groups).transform("sum")
        local = agree.astype(int).groupby(groups).transform("sum")
        same_scale = global_ok & (local >= 3) & (local >= .95 * n_local)
        patch = differ & ~agree & y.notna() & same_scale
        if col in ("bjd", "ra_deg", "dec_deg"):
            patch[:] = False
        # Missing corrected values cannot by themselves document a retraction.
        action = pd.Series("withheld_scale_or_identity", index=x.index)
        action.loc[y.isna()] = "retained_author_cds_missing"
        action.loc[agree] = "kept_author_within_print_precision"
        action.loc[patch] = "used_cds"
        author.loc[dest[patch.to_numpy()], col] = new.loc[patch, col].to_numpy()
        audit["fields"][col] = dict(finite_pairs=int(both.sum()), compatible_fraction=fraction,
                                    global_scale_pass=global_ok, finite_different=int((differ & both).sum()),
                                    beyond_precision_or_missing=int((differ & ~agree).sum()),
                                    applied=int(patch.sum()))
        changed = new.loc[differ, ["recno", "star", "epoch_id", "bjd"]].copy()
        changed["field"] = col
        changed["author_value"] = x[differ]
        changed["table4_value"] = y[differ]
        changed["action"] = action[differ]
        changes.append(changed)
    pd.concat(changes, ignore_index=True).to_csv(out/"correction-changes.csv.gz", index=False, compression=gzip)
    audit["fields_applied"] = sum(item["applied"] for item in audit["fields"].values())
    audit["status"] = "full_author_base_compared_to_partial_cds_asu"
    return author


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


def hd_key(names, require_prefix):
    # Preserve stellar-component suffixes while accepting zero padding and
    # separators between the catalogue prefix, the numeric ID and a component.
    compact = names.astype("string").str.upper().str.strip()
    compact = compact.str.replace(r"(?<=HD)[_\s-]+", "", regex=True)
    compact = compact.str.replace(r"(?<=\d)[_\s-]+(?=[A-Z]$)", "", regex=True)
    pattern = (r"^HD0*([1-9][0-9]*)([A-Z]*)$" if require_prefix else
               r"^(?:HD)?0*([1-9][0-9]*)([A-Z]*)$")
    parts = compact.str.extract(pattern)
    return "HD" + parts[0] + parts[1]


def spores_positions(work, manifest, errors):
    """Resolve alternate survey designations from SPORES's own SIMBAD positions."""
    for host in ASU_HOSTS:
        try:
            params = {"-source": "J/AJ/170/343/table3", "-out": "recno,HD,_RA,_DE",
                      "-out.max": 142, "-out.meta": "h", "-oc.form": "d", "-sort": "recno"}
            path = work/f"spores-identifications.{host}.tsv"
            download(f"https://{host}/viz-bin/asu-tsv?"+urlencode(params), path, manifest,
                     timeout=20, attempts=1)
            frame = pd.read_csv(path, sep="\t", dtype=str, comment="#", quoting=csv.QUOTE_NONE)
            frame.columns = frame.columns.str.strip().str.lower()
            recno = pd.to_numeric(frame["recno"], errors="coerce")
            np.testing.assert_array_equal(np.sort(recno.dropna()), np.arange(1, 142))
            pos = frame.loc[recno.notna(), ["hd", "_ra", "_de"]].copy()
            pos["hd"] = hd_key(pos.hd, False)
            for col in ("_ra", "_de"):
                pos[col] = pd.to_numeric(pos[col], errors="coerce")
            pos = pos[pos["_ra"].between(0, 360) & pos["_de"].between(-90, 90) & pos.hd.notna()]
            if pos.empty:
                raise ValueError("No verified decimal-degree SPORES positions")
            return pos[~pos.hd.duplicated(keep=False)].set_index("hd")
        except Exception as exc:
            errors.append(f"SPORES table3 position crossmatch at {host}: {str(exc)[:300]}")
    return pd.DataFrame(columns=["_ra", "_de"])


def add_sindex(data, work, manifest, errors, out):
    """Require a unique source and target epoch within two seconds."""
    info = {"status": "retrieving"}
    try:
        aux = catalog(SPORES, S_FIELDS, 153902, work, manifest, errors)
        aux["star"] = hd_key(aux.hd, require_prefix=False)
        aux["bjd"] = pd.to_numeric(aux.bjd, errors="coerce")
        aux["sindex"] = pd.to_numeric(aux.sindex, errors="coerce")
        # Preserve raw labels for checking the release's actual vocabulary.
        aux["inst"] = aux.inst.astype("string").str.strip().str.upper()
        aux["positive_s"] = aux.sindex.gt(0) & np.isfinite(aux.sindex)
        counts = aux.groupby("inst", dropna=False).agg(
            rows=("bjd","size"), valid_bjd=("bjd","count"),
            positive_s=("positive_s","sum"), valid_hd=("star","count"))
        info["by_instrument"] = json.loads(counts.reset_index().to_json(orient="records"))
        labels = aux.inst.str.replace(r"[-_\s]+", "", regex=True)
        harps = aux.loc[labels.isin(["HARPS", "HARPSPRE", "HARPSPOST"])].copy()
        era_conflict = ((labels.eq("HARPSPRE") & aux.bjd.ge(FIBRE_CHANGE)) |
                        (labels.eq("HARPSPOST") & aux.bjd.lt(FIBRE_CHANGE)))
        info["harps_conflicting_pre_post_labels"] = int(era_conflict.sum())
        info.update(input_rows=len(aux), harps_rows=len(harps),
                    harps_positive_s=int(harps.positive_s.sum()))
        harps = harps[~era_conflict.loc[harps.index] & harps.positive_s &
                      harps.bjd.between(2452700,2465000) & harps.star.notna()]
        stars = data.groupby("star")[["ra_deg","dec_deg"]].median()
        names = hd_key(stars.index.to_series(), require_prefix=True)
        common = set(harps.star.dropna()) & set(names.dropna())
        info.update(harps_valid_s_bjd_hd=len(harps), direct_hd_targets=len(common))
        if len(harps):
            info["harps_valid_s_bjd_range"] = [float(harps.bjd.min()), float(harps.bjd.max())]
        alt = {}
        if len(set(harps.star)-common):
            pos = spores_positions(work, manifest, errors)
            # A source must be the only survey target within two arcseconds.
            # Also require that only one SPORES source can map to that target.
            for hd, row in pos.iterrows():
                dra = (stars.ra_deg-row["_ra"]+180) % 360-180
                dist = np.hypot(dra*np.cos(np.deg2rad(row["_de"])), stars.dec_deg-row["_de"])
                near = stars.index[dist.le(2/3600)]
                if len(near) == 1:
                    alt.setdefault(near[0], []).append(hd)
            for star, ids in alt.items():
                if len(ids) == 1 and pd.isna(names.at[star]) and ids[0] not in common:
                    names.at[star] = ids[0]
        info["position_verified_alternate_targets"] = int(names.notna().sum()
            - hd_key(stars.index.to_series(), True).notna().sum())
        name = data.star.map(names)
        nearest, matched = [], []
        for star, group in harps.groupby("star"):
            left = data.loc[name.eq(star) & data.bjd.notna(), "bjd"].sort_values()
            if left.empty:
                continue
            right = group.sort_values("bjd")
            times = right.bjd.to_numpy()
            limit = 2 / 86400
            t = left.to_numpy()
            lo = np.searchsorted(times, t-limit)
            hi = np.searchsorted(times, t+limit, side="right")
            # Both directions need exactly one neighbouring exposure.
            reverse_lo = np.searchsorted(t, times-limit)
            reverse_hi = np.searchsorted(t, times+limit, side="right")
            unique = (hi-lo == 1) & ((reverse_hi-reverse_lo)[lo.clip(0, len(times)-1)] == 1)
            data.loc[left.index[unique], "sindex"] = right.sindex.to_numpy()[lo[unique]]
            m = np.searchsorted(t, times)
            nearest.extend(np.minimum(np.abs(times-t[(m-1).clip(0,len(t)-1)]),
                                      np.abs(times-t[m.clip(0,len(t)-1)]))*86400)
            matched.extend(np.abs(t[unique]-times[lo[unique]])*86400)
        info.update(status="completed", mutual_matches=len(matched), nearest_for_common_sources=len(nearest))
        if nearest:
            info["unthresholded_nearest_rvbank_seconds_q0_q25_q50_q90"] = np.percentile(nearest,[0,25,50,90]).tolist()
            info["nearest_within_2s_10s_90s"] = [int((np.asarray(nearest)<=x).sum()) for x in (2,10,90)]
        if matched:
            info["max_accepted_distance_seconds"] = float(max(matched))
    except Exception as exc:
        info.update(status="failed", error=str(exc))
        raise
    finally:
        (out/"sindex-coverage.json").write_text(json.dumps(info, indent=2)+"\n")


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
    if source in (PRIMARY, AUTHOR):
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


def select_archive(source, work, manifest, errors, out, audit):
    if source == "auto":
        subsets = []
        try:
            primary = catalog(PRIMARY, RV_FIELDS, 252615, work, manifest, errors, audit, subsets)
            audit.update(status="full_cds_table", input_rows=len(primary))
            return PRIMARY, primary
        except RuntimeError:
            try:
                author = author_rvbank(work, manifest)
            except Exception as exc:
                if manifest and manifest[-1]["url"] == AUTHOR_URL and manifest[-1]["status"] == "downloaded":
                    manifest[-1].update(status="unusable_response", error=str(exc)[:500])
                errors.append(f"{AUTHOR}: {exc}")
            else:
                if len(author) != 252615:
                    raise ValueError(f"Expected 252615 author exposures, received {len(author)}")
                if subsets:
                    return AUTHOR, apply_partial_table(author, subsets[0], audit, out)
                audit.update(status="author_only_no_verified_cds_coverage", author_rows=len(author))
                return AUTHOR, author
    return ESO, eso_catalog(work, manifest)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--source", choices=("auto", "eso"), default="auto",
                        help="Use ESO directly for the time-separated RV holdout")
    args = parser.parse_args()
    require_survey_workflow()
    work, out = args.work, args.out
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    manifest, errors = [], []
    release = ESO if args.source == "eso" else "auto"
    audit = dict(published_table4_records=252615, catalogue=PRIMARY, status="unavailable")
    try:
        release, data = select_archive(args.source, work, manifest, errors, out, audit)
        print(f"Selected {release}: {len(data)} input rows", flush=True)
        classify(data, release)
        if not data.disposition.eq("used").any():
            raise ValueError("no usable primary RV measurements")
        if args.source != "eso":
            try:
                add_sindex(data, work, manifest, errors, out)
            except Exception as exc:
                data["sindex"] = np.nan
                errors.append(f"S-index retrieval or matching unavailable: {exc}")
        if release in (PRIMARY, AUTHOR):
            try:
                add_published_sindex(data, work, out, manifest, errors, catalog)
            except Exception as exc:
                errors.append(f"Ibáñez Bustos S-index supplement unavailable: {exc}")
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
        version = {PRIMARY: "CDS corrected 2026-04-20", AUTHOR: "author Ver.02b 2024-03-03",
                   ESO: "ESO Phase 3 DR1 2023-12-08"}[release]
        info = dict(source=release, release_version=version, cds_corrected_table=(release == PRIMARY),
                    correction_status=audit["status"] if args.source == "auto" else None,
                    partial_cds_fields_applied=audit.get("fields_applied", 0),
                    author_commit=AUTHOR_COMMIT if release == AUTHOR else None,
                    published_input_rows=289843 if release == ESO else 252615,
                    input_rows=len(data), targets=len(targets),
                    used_epochs=len(used), dispositions=data.disposition.value_counts().to_dict(),
                    sindex_epochs=int(used.sindex.notna().sum()),
                    first_bjd=float(used.bjd.min()), last_bjd=float(used.bjd.max()))
    except Exception as exc:
        errors.append(str(exc))
        info = dict(source=release, status="ingest_failed")
        raise
    finally:
        if args.source == "auto":
            (out/"correction-audit.json").write_text(json.dumps(audit, indent=2)+"\n")
        pd.DataFrame(manifest, columns=["url", "retrieved_utc", "file", "attempt", "status", "bytes", "sha256", "error"]).to_csv(
            out/"inputs.csv", index=False)
        (out/"ingest.json").write_text(json.dumps(dict(**info, warnings=errors), indent=2)+"\n")
        print(json.dumps(dict(**info, warnings=errors)), flush=True)


if __name__ == "__main__":
    main()
