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
           "instrument_name", "spectro_smw", "match_method", "bjd_difference_seconds"]


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
    missing = {"ins_mode", "instrument_name", "spectro_smw"} - columns.keys()
    if missing:
        raise ValueError("DACE response lacks " + ", ".join(sorted(missing)))
    if not any(len(columns.get(key, [])) for key in ("rjd", "file_rootname")) and len(columns["ins_mode"]):
        raise ValueError("DACE response lacks both exposure identifiers and barycentric times")
    for key in ("rjd", "file_rootname"):
        if not len(columns.get(key, [])):
            columns[key] = [None] * len(columns["ins_mode"])
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


def timestamps(values):
    # Prefixes and suffixes differ between ESO filenames and RVBank TimeId.
    # Require at least millisecond precision. Normalize only numerically equal
    # fractional seconds (e.g. .1230 and .123), never a time offset.
    parts = values.astype("string").str.extract(
        r"(?<!\d)(\d{4}-\d{2}-\d{2}[T ]\d{2}:\d{2}:\d{2})\.(\d{3,9})(?!\d)")
    return parts[0].str.replace(" ", "T", regex=False) + "." + parts[1].str.pad(
        9, side="right", fillchar="0")


def match_exposures(star, epoch, dace, receipt):
    # DACE calls rjd the barycentric reduced Julian date. This is BJD-2400000,
    # not MJD-OBS. The request is scoped to the target, southern HARPS and
    # public, quality-passing, latest-DRS POSTDRS_A observations.
    left = epoch[["epoch_id", "bjd", "era"]].reset_index(drop=True).copy()
    right = dace.reset_index(drop=True).copy()
    left["rjd_rvbank"] = left.bjd - 2400000.
    left["stamp"] = timestamps(left.epoch_id)
    right["stamp"] = timestamps(right.file_rootname)
    right["rjd"] = pd.to_numeric(right.rjd, errors="coerce")
    right["spectro_smw"] = pd.to_numeric(right.spectro_smw, errors="coerce")
    mode = right.ins_mode.astype("string").str.strip().str.upper()
    instrument = right.instrument_name.astype("string").str.strip().str.upper()
    positive = np.isfinite(right.spectro_smw) & right.spectro_smw.gt(0)
    receipt.update(raw_mode_counts=json.dumps(right.ins_mode.fillna("<missing>").value_counts().to_dict()),
                   raw_instrument_counts=json.dumps(right.instrument_name.fillna("<missing>").value_counts().to_dict()),
                   positive_s_rows=int(positive.sum()),
                   rvbank_timestamp_rows=int(left.stamp.notna().sum()),
                   dace_timestamp_rows=int(right.stamp.notna().sum()),
                   finite_rjd_rows=int(np.isfinite(right.rjd).sum()))
    all_rjd = right.rjd[np.isfinite(right.rjd)]
    if not all_rjd.empty:
        times = np.unique(all_rjd)
        receipt.update(rjd_min=times[0], rjd_max=times[-1],
                       rjd_examples=json.dumps(times[np.linspace(0,len(times)-1,min(3,len(times)),
                                                                 dtype=int)].tolist()),
                       # Float serialization strips trailing zeroes. This is a
                       # precision diagnostic, not an estimate of the clock error.
                       rjd_decimals_shown=json.dumps(all_rjd.map(lambda value:
                           len(np.format_float_positional(value,trim="-").partition(".")[2])
                       ).value_counts().sort_index().to_dict()))
    # Keep explicit instrument/mode validation in addition to the server filters.
    valid_instrument = instrument.isin(["HARPS", "HARPS03", "HARPS15"])
    valid_mode = mode.isin(["HARPS", "HAM"])
    receipt.update(instrument_valid_rows=int(valid_instrument.sum()),
                   mode_valid_rows=int(valid_mode.sum()),
                   usable_s_rows=int((valid_instrument & valid_mode & positive).sum()))
    right["named_era"] = instrument.map({"HARPS03": "PRE", "HARPS15": "POST"})
    right = right[valid_instrument & valid_mode & positive].copy()
    # The match can be proved by instrument timestamp even when reduction
    # programs calculated slightly different barycentric mid-exposure times.
    # Explicit PRE/POST disagreement always rejects the association.
    unique_left = left[left.stamp.notna() & ~left.stamp.duplicated(keep=False)]
    unique_right = right[right.stamp.notna() & ~right.stamp.duplicated(keep=False)]
    exact = unique_left.assign(i=unique_left.index).merge(
        unique_right.assign(j=unique_right.index), on="stamp", validate="one_to_one")
    receipt["exact_timestamp_pairs"] = len(exact)
    era_ok = exact.named_era.isna() | exact.era.eq(exact.named_era)
    selected = [(row.i, row.j, "unique exposure timestamp")
                for row in exact[era_ok].itertuples(index=False)]
    used_left = {i for i, _, _ in selected}
    used_right = {j for _, j, _ in selected}

    # Test BJD uniqueness in the full series, including any rows already
    # paired by timestamp, to avoid resolving an ambiguity through row removal.
    a = left[np.isfinite(left.rjd_rvbank)].sort_values("rjd_rvbank")
    b = right[np.isfinite(right.rjd)].sort_values("rjd")
    at, bt = a.rjd_rvbank.to_numpy(), b.rjd.to_numpy()
    tolerance = 2 / 86400.
    if len(at) and len(bt):
        lo = np.searchsorted(bt, at-tolerance)
        hi = np.searchsorted(bt, at+tolerance, side="right")
        reverse_lo = np.searchsorted(at, bt-tolerance)
        reverse_hi = np.searchsorted(at, bt+tolerance, side="right")
        mutual = (hi-lo == 1) & ((reverse_hi-reverse_lo)[lo.clip(0,len(bt)-1)] == 1)
        nearest_idx = np.searchsorted(bt,at)
        nearest = np.minimum(abs(bt[(nearest_idx-1).clip(0,len(bt)-1)]-at),
                             abs(bt[nearest_idx.clip(0,len(bt)-1)]-at))*86400
        receipt["nearest_rjd_seconds_q0_q25_q50_q90"] = json.dumps(
            np.percentile(nearest,[0,25,50,90]).tolist())
        receipt["mutually_unique_bjd_pairs"] = int(mutual.sum())
        for ia in np.flatnonzero(mutual):
            i, j = a.index[ia], b.index[lo[ia]]
            if i in used_left or j in used_right:
                continue
            lhs, rhs = left.loc[i], right.loc[j]
            if pd.notna(rhs.named_era) and rhs.named_era != lhs.era:
                continue
            if pd.notna(lhs.stamp) and pd.notna(rhs.stamp) and lhs.stamp != rhs.stamp:
                continue  # A conflicting raw exposure identifier overrides a near BJD.
            selected.append((i, j, "mutual barycentric time within 2 s"))
    paired = []
    for i, j, origin in selected:
        lhs, rhs = left.loc[i], right.loc[j]
        paired.append(dict(star=star, epoch_id=lhs.epoch_id, bjd=lhs.bjd,
                           file_rootname=rhs.file_rootname, rjd=rhs.rjd, ins_mode=rhs.ins_mode,
                           instrument_name=rhs.instrument_name, spectro_smw=rhs.spectro_smw,
                           match_method=origin, bjd_difference_seconds=(lhs.rjd_rvbank-rhs.rjd)*86400))
    receipt["accepted_by_method"] = json.dumps(pd.Series([o for _,_,o in selected],
                                                          dtype=str).value_counts().to_dict())
    return pd.DataFrame(paired, columns=COLUMNS)


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


def initial_sindex_sources(work, out):
    """Publish the source for each nonmissing S-index in the actual RV input."""
    used = pd.read_csv(work / "epochs.csv.gz", dtype={"star": str, "epoch_id": str},
                       usecols=["star", "epoch_id", "bjd", "sindex"])
    used = used[used.sindex.notna()].copy()
    used["sindex_source"] = "SPORES-HWO II Table 5"
    supplement = out / "ibanez-sindex-matches.csv"
    if supplement.exists():
        matches = pd.read_csv(supplement, dtype={"star": str, "epoch_id": str})
        if not matches.empty:
            keys = ["star", "epoch_id", "bjd"]
            joined = matches.merge(used, on=keys, how="left", validate="one_to_one",
                                   suffixes=("_published", "_input"))
            np.testing.assert_allclose(joined.sindex_published, joined.sindex_input,
                                       rtol=1e-12, atol=0)
            used.loc[used.star.isin(matches.star), "sindex_source"] = "Ibáñez Bustos et al. Table B.1"
            if used.star.isin(matches.star).sum() != len(matches):
                raise ValueError("Supplemented target contains unaccounted S-index measurements")
    used.sort_values(["star", "bjd"]).to_csv(out / "ingest-sindex-epochs.csv.gz",
        index=False, compression={"method": "gzip", "mtime": 0})
    return used[["star", "sindex_source"]].drop_duplicates().set_index("star", verify_integrity=True).sindex_source


def main(work, out):
    require_survey_workflow()
    cache = work / "catalog-cache"
    cache.mkdir(parents=True, exist_ok=True)
    signals = pd.read_csv(out / "signals.csv", dtype={"star": str}, keep_default_na=False)
    signals["initial_activity_json"] = signals.activity_json
    signals["initial_flags"] = signals["flags"]
    signals["initial_rv_status"] = signals.rv_status
    signals["sindex_source"] = signals.star.map(initial_sindex_sources(work, out)).fillna("unavailable")
    signals["initial_sindex_source"] = signals.sindex_source
    wanted = signals.loc[signals["flags"].eq("missing_sindex"), "star"].unique()
    if len(wanted):
        epochs = pd.read_csv(work / "epochs.csv.gz", dtype={"star": str})
        injections = pd.read_csv(out / "injections.csv", dtype={"star": str})
        n_targets = len(pd.read_csv(out / "targets.csv", usecols=["star"]))
    headers = ["star", "status", "utc", "query_url", "attempts", "http_status",
               "response_bytes", "sha256", "dace_rows", "matched_epochs", "nights", "error",
               "raw_mode_counts", "raw_instrument_counts", "positive_s_rows", "usable_s_rows",
               "rvbank_timestamp_rows", "dace_timestamp_rows", "finite_rjd_rows", "rjd_min", "rjd_max",
               "rjd_examples", "rjd_decimals_shown", "instrument_valid_rows", "mode_valid_rows",
               "exact_timestamp_pairs", "mutually_unique_bjd_pairs", "nearest_rjd_seconds_q0_q25_q50_q90",
               "accepted_by_method"]
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
                paired = match_exposures(star, frame, table, receipt)
                nights = len(np.unique(np.floor(paired.bjd)))
                receipt.update(matched_epochs=len(paired), nights=nights,
                               status="insufficient matched DACE series")
                matches.append(paired)
                if nights >= 10 and paired.spectro_smw.nunique() > 1:
                    # Analyze the homogeneous DACE S time series alone. Mixing with
                    # another published reduction could introduce a systematic offset.
                    frame["sindex"] = frame.bjd.map(paired.set_index("bjd", verify_integrity=True).spectro_smw)
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
