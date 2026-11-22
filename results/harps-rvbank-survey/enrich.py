"""Recover Mount Wilson S-indices on exact HARPS exposures for otherwise viable RV signals."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.parse import quote

import numpy as np
import pandas as pd
import requests

from runtime import require_survey_workflow
from search import search


# Public endpoint and filter schema used by astro-dace/dace-query 3.0.
ENDPOINT = "https://spectroscopy-webapp-pg.obsuksprd2.unige.ch/target/{}/timeseries/radial-velocities"
FILTERS = {"instrument_group": {"equals": ["HARPS"]}, "drs_qc": {"is": True},
           "status_public": {"is": True}, "is_latest_drs": {"is": True},
           "source_product_file_ext": {"equals": ["POSTDRS_A"]}}
LIMIT = 40000
FIELDS = {"file_rootname", "rjd", "ins_mode", "instrument_name", "spectro_smw"}
COLUMNS = ["star", "epoch_id", "bjd", "file_rootname", "rjd", "ins_mode",
           "instrument_name", "spectro_smw"]


def decompress_columns(data):
    """Read DACE's column vectors and expand run-length encoded values."""
    columns = {}
    for item in data["parameters"]:
        key = item["variableName"]
        if key in FIELDS:
            values = next((item[k] for k in ("doubleValues", "floatValues", "stringValues")
                           if k in item), [])
            repeat = item.get("occurrences")
            columns[key] = np.repeat(values, repeat) if repeat else values
    if not columns:
        return pd.DataFrame(columns=sorted(FIELDS))
    if not FIELDS.issubset(columns):
        raise ValueError("DACE response lacks " + ", ".join(sorted(FIELDS - columns.keys())))
    return pd.DataFrame(columns)


def fetch(session, star, receipt, cache):
    url = ENDPOINT.format(quote(star, safe=""))
    params = dict(limit=LIMIT, filters=json.dumps(FILTERS, sort_keys=True), sort="{}")
    receipt["query_url"] = requests.Request("GET", url, params=params).prepare().url
    for attempt in (1, 2):
        receipt.update(utc=datetime.now(timezone.utc).isoformat(), attempts=attempt)
        try:
            r = session.get(url, params=params, timeout=(10, 60))
            receipt.update(http_status=r.status_code, response_bytes=len(r.content),
                           sha256=hashlib.sha256(r.content).hexdigest())
            (cache / receipt["sha256"]).write_bytes(r.content)
            r.raise_for_status()
            table = decompress_columns(r.json())
            receipt["dace_rows"] = len(table)
            if len(table) >= LIMIT:
                raise ValueError("DACE response reached the query limit")
            return table
        except (requests.RequestException, KeyError, ValueError) as exc:
            receipt["error"] = str(exc)[:350]
            receipt["status"] = "DACE unavailable"
            if (isinstance(exc, requests.HTTPError) and r.status_code < 500
                    and r.status_code != 429) or attempt == 2:
                return None
            time.sleep(2)


def match_exposures(star, epoch, dace):
    # TimeId in RVBank is the instrument timestamp (ISO date with milliseconds).
    # Require it in the raw HARPS frame name, with an independent BJD test.
    stamp = r"(\d{4}-\d\d-\d\d[T ]\d\d:\d\d:\d\d\.\d{3})(?!\d)"
    left = epoch[["epoch_id", "bjd", "era"]].copy()
    left["stamp"] = left.epoch_id.astype("string").str.extract(stamp, expand=False).str.replace(" ", "T")
    right = dace.copy()
    right["stamp"] = right.file_rootname.astype("string").str.extract(
        "HARPS[.]" + stamp, expand=False).str.replace(" ", "T")
    right["spectro_smw"] = pd.to_numeric(right.spectro_smw, errors="coerce")
    right["rjd"] = pd.to_numeric(right.rjd, errors="coerce")
    right = right[right.ins_mode.str.upper().isin(["HAM", "HARPS"])
                  & right.instrument_name.str.upper().isin(["HARPS", "HARPS03", "HARPS15"])
                  & np.isfinite(right.spectro_smw) & right.spectro_smw.gt(0)]
    left = left.dropna(subset=["stamp"]).drop_duplicates("stamp", keep=False)
    right = right.dropna(subset=["stamp"]).drop_duplicates("stamp", keep=False)
    paired = left.merge(right, on="stamp", validate="one_to_one")
    tdiff = abs((paired.bjd - 2400000.) - paired.rjd)
    names = paired.instrument_name.str.upper()
    era_ok = ~(names.eq("HARPS03") & paired.era.ne("PRE") |
               names.eq("HARPS15") & paired.era.ne("POST"))
    paired = paired[(tdiff <= 2 / 86400) & era_ok].copy()
    paired["star"] = star
    return paired[COLUMNS]


def validate_repeat(old, new, cells, original):
    old = old.sort_values("signal_index")
    new = new.sort_values("signal_index")
    np.testing.assert_array_equal(old.signal_index, new.signal_index)
    for key in ("period_d", "k_ms", "fap", "bootstrap_exceedances"):
        np.testing.assert_allclose(old[key], new[key], rtol=1e-8, atol=1e-10)
    dims = ["period_d", "k_ms", "eccentricity", "trials", "recovered"]
    a = cells.sort_values(dims[:3])[dims].to_numpy(dtype=float)
    b = original.sort_values(dims[:3])[dims].to_numpy(dtype=float)
    np.testing.assert_array_equal(a, b)


def main(work, out):
    require_survey_workflow()
    cache = work / "catalog-cache"
    cache.mkdir(parents=True, exist_ok=True)
    signals = pd.read_csv(out / "signals.csv", dtype={"star": str}, keep_default_na=False)
    signals["initial_activity_json"] = signals.activity_json
    signals["initial_flags"] = signals["flags"]
    signals["initial_rv_status"] = signals.rv_status
    signals["sindex_source"] = "SPORES if covered"
    wanted = signals.loc[signals["flags"].eq("missing_sindex"), "star"].unique()
    if len(wanted):
        epochs = pd.read_csv(work / "epochs.csv.gz", dtype={"star": str})
        injections = pd.read_csv(out / "injections.csv", dtype={"star": str})
        n_targets = len(pd.read_csv(out / "targets.csv", usecols=["star"]))
    headers = ["star", "status", "utc", "query_url", "attempts", "http_status",
               "response_bytes", "sha256", "dace_rows", "matched_epochs", "nights", "error"]
    matches = []
    with requests.Session() as session, (out / "sindex-receipts.csv").open("w", newline="") as handle:
        session.headers.update({"Accept":"application/json", "User-Agent":"Open-Space-HARPS-survey/1.0"})
        log = csv.DictWriter(handle, fieldnames=headers)
        log.writeheader()
        for star in sorted(wanted):
            receipt = {"star": star}
            table = fetch(session, star, receipt, cache)
            if table is not None:
                frame = epochs[epochs.star.eq(star)].copy()
                paired = match_exposures(star, frame, table)
                nights = len(np.unique(np.floor(paired.bjd)))
                receipt.update(matched_epochs=len(paired), nights=nights,
                               status="insufficient matched DACE series")
                matches.append(paired)
                if nights >= 10 and paired.spectro_smw.nunique() > 1:
                    # Analyze the homogeneous DACE S time series alone. Mixing with
                    # SPORES could introduce a systematic offset between reductions.
                    frame["sindex"] = frame.epoch_id.map(paired.set_index("epoch_id").spectro_smw)
                    _, repeated, cells = search(star, frame, n_targets)
                    new = pd.DataFrame(repeated)
                    selection = signals.star.eq(star)
                    validate_repeat(signals[selection], new, pd.DataFrame(cells),
                                    injections[injections.star.eq(star)])
                    for row in new.itertuples(index=False):
                        hit = selection & signals.signal_index.eq(row.signal_index)
                        # Only update signals whose sole unresolved internal check was
                        # the S-index. Earlier and later components stay in the ledger.
                        if signals.loc[hit, "initial_flags"].item() == "missing_sindex":
                            signals.loc[hit, ["activity_json", "flags", "rv_status",
                                              "sindex_source"]] = [
                                row.activity_json, row.flags, row.rv_status, "DACE matched exposures"]
                    receipt["status"] = "repeated with DACE S-index"
            log.writerow(receipt)
            handle.flush()
    pd.concat([pd.DataFrame(columns=COLUMNS), *matches], ignore_index=True).to_csv(
        out / "sindex-epochs.csv.gz", index=False, compression={"method":"gzip","mtime":0})
    updated = out / "signals.next.csv"
    signals.to_csv(updated, index=False)
    updated.replace(out / "signals.csv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.work, args.out)
