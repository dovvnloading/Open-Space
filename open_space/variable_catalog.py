"""Crossmatch variability evidence against VSX, Gaia DR3, TESS-SVC and TARS."""

import argparse
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing, contextmanager
import csv
import fcntl
import hashlib
import io
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile
from urllib.error import URLError
from urllib.parse import urlencode

import zipfile

from .archive_index import retrieve, record_response
from .tess_archive import atomic_bytes, download
from .tess_eb_catalog import separation

VIZIER = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
MAX_ROWS = 1000
CONE_TABLES = {
    "vsx": ("B/vsx/vsx", "OID,Name,V,RAJ2000,DEJ2000,Type,l_Period,Period,u_Period"),
    "gaia_summary": ("I/358/varisum", "Source,RA_ICRS,DE_ICRS"),
    "gaia_classification": ("I/358/vclassre", "Source,RA_ICRS,DE_ICRS,Class,ClassSc"),
}
SVC_FILE = "hlsp_tess-svc_tess_lcf_all-s0001-s0026_tess_v1.0_cat.csv"
TARS_FILE = "tars_table_2.csv.zip"
# Hash, expected row count and columns retained in the local compact index.
RELEASES = {
    "tess_svc": (
        "https://archive.stsci.edu/hlsps/tess-svc/" + SVC_FILE,
        "6c0ba9aa391b007f9c803f3d847c387adfcf7a7b65c68142fd2440ec1936d7c4",
        84046,
        SVC_FILE,
        "tess_id,ra,dec,Solution,Sector,period_var_1,period_var_uncert_1,period_var_2,period_var_uncert_2",
    ),
    "tars": (
        "https://zenodo.org/records/19917941/files/tars_table_2.csv.zip?download=1",
        "5f774c42d440516c3a679d288fe60f4877d6cd695019e93127bbe704b38d4df6",
        1046317,
        TARS_FILE,
        "TICID,ra,dec,dr3_source_id,adopted_period,adopted_period_unc,flag_multiple_periods,"
        "flag_possible_binary,final_n_contams,n_secs,flag_doubled_period",
    ),
}


def _finite(value):
    if value is None or not str(value).strip():
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def period_comparison(measured, published):
    """Closest of P/2, P, 2P, reported explicitly without inferring an identification."""
    if measured is None or published is None or published <= 0:
        return None
    factor = min((0.5, 1.0, 2.0), key=lambda f: abs(math.log(measured / (f * published))))
    return {
        "measured_over_published": measured / published,
        "closest_factor": factor,
        "fractional_offset_from_factor": measured / (factor * published) - 1,
    }


def _parse_vizier(raw, table, columns):
    lines = raw.decode("utf-8").splitlines()
    marker = "#Name: " + table
    if marker not in lines or any(
        line.lower().startswith(("#error", "#warning", "#info\terror", "#info\twarning"))
        for line in lines
    ):
        raise ValueError(f"VizieR query failed or is missing table {table}")
    body = lines[lines.index(marker) + 1 :]
    names = [line.split("\t")[1] for line in body if line.startswith("#Column\t")]
    if names != columns.split(","):
        raise ValueError(f"VizieR schema changed for {table}")
    # Empty ASU tables contain validated column descriptions and no tabular body.
    data = [line for line in body if line.strip() and not line.startswith("#")]
    rows = []
    if data:
        if data[0].split("\t") != names or len(data) < 3 or not data[2].startswith("---"):
            raise ValueError("Invalid VizieR tabular header")
        for line in data[3:]:
            cells = [v.strip() for v in line.split("\t")]
            rows.append(dict(zip(names, cells, strict=True)))
    dates = [
        line.split("=", 1)[1].split("\t")[0] for line in lines if line.startswith("#INFO\trequest_date=")
    ]
    if len(dates) != 1:
        raise ValueError("VizieR query date missing")
    return rows, {
        "table": table,
        "description": next((l[8:] for l in body if l.startswith("#Title: ")), None),
        "query_date_utc": dates[0],
    }


def _cone(key, ra, dec, radius, period, cache):
    table, columns = CONE_TABLES[key]
    url = (
        VIZIER
        + "?"
        + urlencode(
            {
                "-source": table,
                "-c": f"{ra:.8f} {dec:+.8f}",
                "-c.rs": f"{radius:g}",
                "-out": columns,
                "-out.max": MAX_ROWS,
            }
        )
    )
    raw, path = retrieve(url, cache, "text/tab-separated-values")
    rows, receipt = _parse_vizier(raw, table, columns)
    matches = []
    for row in rows:
        if key == "vsx":
            x, y = float(row["RAJ2000"]), float(row["DEJ2000"])
            published = _finite(row["Period"])
            match = {
                "source_id": row["OID"],
                "name": row["Name"],
                "variability_flag": row["V"],
                "class": row["Type"] or None,
                "period_days": published,
                "period_limit_flag": row["l_Period"] or None,
                "period_uncertainty_flag": row["u_Period"] or None,
                "period_comparison": period_comparison(period, published),
            }
        else:
            x, y = float(row["RA_ICRS"]), float(row["DE_ICRS"])
            match = {"source_id": row["Source"]}
            if key == "gaia_classification":
                match.update({"class": row["Class"], "class_score": _finite(row["ClassSc"])})
        matches.append(
            {**match, "ra_deg": x, "dec_deg": y, "separation_arcsec": separation(ra, dec, x, y)}
        )
    return {
        "status": "possibly_truncated" if len(matches) >= MAX_ROWS else "ok",
        "matches": sorted(matches, key=lambda row: row["separation_arcsec"]),
        "receipt": {**receipt, **record_response(raw, path, url)},
    }


@contextmanager
def _catalog_rows(path):
    if path.name.endswith(".zip"):
        with zipfile.ZipFile(path) as archive, archive.open("tars_table_2.csv") as zipped:
            with io.TextIOWrapper(zipped, newline="", encoding="utf-8") as handle:
                yield csv.reader(handle)
    else:
        with path.open(newline="", encoding="utf-8") as handle:
            yield csv.reader(handle)


def _index(key, cache):
    url, digest, count, filename, columns = RELEASES[key]
    index = cache / f"{key}-{digest}-v1.sqlite"
    cache.mkdir(parents=True, exist_ok=True)
    # Workers can safely share a partially downloaded/rebuilt cache. A crashed
    # builder never publishes its SQLite file and the OS releases the lock.
    with (cache / (key + ".lock")).open("a+b") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if index.exists():
            return index
        source = cache / filename
        raw = source.read_bytes() if source.exists() else download(url, 145_000_000)
        if hashlib.sha256(raw).hexdigest() != digest:
            raise ValueError(f"{key} source checksum differs from the pinned release")
        if not source.exists():
            atomic_bytes(source, raw)
        del raw
        names = columns.split(",")
        with tempfile.NamedTemporaryFile(dir=cache, delete=False) as tmp:
            staged = Path(tmp.name)
        try:
            with closing(sqlite3.connect(staged)) as conn, _catalog_rows(source) as rows:
                header = next(rows)
                offsets = [header.index(field) for field in names]
                conn.execute("PRAGMA journal_mode=OFF")
                conn.execute("PRAGMA synchronous=OFF")
                conn.execute("CREATE TABLE sources (" + ",".join(names) + ")")
                conn.executemany(
                    "INSERT INTO sources VALUES (" + ",".join("?" for _ in names) + ")",
                    ([row[i] for i in offsets] for row in rows),
                )
                if conn.execute("SELECT count(*) FROM sources").fetchone()[0] != count:
                    raise ValueError(f"{key} row count differs from pinned release")
                conn.execute(f"CREATE INDEX tic ON sources({names[0]})")
                conn.commit()
            os.replace(staged, index)
        finally:
            staged.unlink(missing_ok=True)
    return index


def _tic_catalog(key, tic, ra, dec, period, cache):
    url, sha, count, name, columns = RELEASES[key]
    path = _index(key, cache)
    with closing(sqlite3.connect(path)) as conn:
        conn.row_factory = sqlite3.Row
        rows = list(conn.execute(f"SELECT * FROM sources WHERE {columns.split(',')[0]}=?", (str(tic),)))
    matches = []
    for row in rows:
        if key == "tess_svc":
            metadata = {"sector": row["Sector"], "solution": row["Solution"]}
            periods = [
                (_finite(row[f"period_var_{i}"]), _finite(row[f"period_var_uncert_{i}"])) for i in (1, 2)
            ]
        else:
            metadata = {
                "gaia_dr3_id": row["dr3_source_id"],
                "sectors": int(float(row["n_secs"])),
                "possible_binary": row["flag_possible_binary"] == "True",
                "multiple_periods": row["flag_multiple_periods"] == "True",
                "doubled_period": row["flag_doubled_period"] == "True",
                "contaminants": _finite(row["final_n_contams"]),
            }
            periods = [(_finite(row["adopted_period"]), _finite(row["adopted_period_unc"]))]
        p = [
            {
                "period_days": value,
                "uncertainty_days": err,
                "period_comparison": period_comparison(period, value),
            }
            for value, err in periods
            if value is not None
        ]
        x, y = float(row["ra"]), float(row["dec"])
        matches.append(
            {
                "tic_id": tic,
                "same_id": True,
                "ra_deg": x,
                "dec_deg": y,
                "separation_arcsec": separation(ra, dec, x, y),
                "periods": p,
                **metadata,
            }
        )
    return {
        "status": "ok",
        "matches": matches,
        "receipt": {
            "url": url,
            "sha256": sha,
            "rows": count,
            "selection": "exact TIC ID",
        },
    }


def search(ra_deg, dec_deg, cache_dir, radius_arcsec=20, *, tic_id=None, period_days=None):
    """Return each archive outcome separately. An empty successful query has matches=[]."""
    if not (0 <= ra_deg < 360 and -90 <= dec_deg <= 90 and 0 < radius_arcsec <= 180):
        raise ValueError("Invalid cone position or radius (max 180 arcsec)")
    if period_days is not None and not (math.isfinite(period_days) and period_days > 0):
        raise ValueError("Period must be finite and positive")
    if tic_id is not None and (type(tic_id) is not int or tic_id < 1):
        raise ValueError("TIC ID must be a positive integer")
    cache = Path(cache_dir) / "variable-catalog"
    result = {"ra_deg": ra_deg, "dec_deg": dec_deg, "radius_arcsec": radius_arcsec, "tic_id": tic_id}
    for key in (*CONE_TABLES, *RELEASES):
        try:
            if key in CONE_TABLES:
                result[key] = _cone(key, ra_deg, dec_deg, radius_arcsec, period_days, cache)
            elif tic_id is None:
                result[key] = {"status": "requires_tic_id", "matches": []}
            else:
                result[key] = _tic_catalog(key, tic_id, ra_deg, dec_deg, period_days, cache)
        except (OSError, URLError, ValueError, sqlite3.Error, zipfile.BadZipFile) as exc:
            result[key] = {"status": "error", "error": str(exc), "matches": []}
    result["status"] = (
        "ok"
        if all(result[k]["status"] in {"ok", "requires_tic_id"} for k in (*CONE_TABLES, *RELEASES))
        else "incomplete"
    )
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("evidence", type=Path)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--radius-arcsec", type=float, default=42)
    p.add_argument("--workers", type=int, choices=range(1, 9), default=3)
    args = p.parse_args()
    raw = args.evidence.read_bytes()
    all_targets = json.loads(raw)["targets"]
    targets = sorted(
        (t for t in all_targets if t["status"] == "two_sector_coherent_variability"),
        key=lambda t: t["tic_id"],
    )

    def check(t):
        evidence = search(
            t["ra_deg"],
            t["dec_deg"],
            args.cache,
            args.radius_arcsec,
            tic_id=t["tic_id"],
            period_days=t["period_days"],
        )
        return {"role": t["role"], "measured_period_days": t["period_days"], **evidence}

    with ThreadPoolExecutor(args.workers) as pool:
        results = list(pool.map(check, targets))
    report = {
        "input_sha256": hashlib.sha256(raw).hexdigest(),
        "targets_screened": len(all_targets),
        "coherent_targets": len(targets),
        "targets": results,
    }
    atomic_bytes(args.output, (json.dumps(report, indent=2, allow_nan=False) + "\n").encode())
    if any(t["status"] != "ok" for t in results):
        raise SystemExit(
            "Incomplete catalog response. The output records each affected query. Rerun to retry."
        )


if __name__ == "__main__":
    main()
