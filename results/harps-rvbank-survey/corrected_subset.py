"""Build a self-contained RV input from the verified April 2026 CDS Table 4 prefix."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from ingest import COLUMNS, PRIMARY, RV_FIELDS, classify, field_projection
from runtime import require_survey_workflow


SOURCE = PRIMARY + " (verified corrected subset)"
PAGES = [
    (1, 25000, "vizier.cfa.harvard.edu.1-25000", "downloaded",
     "9720c636b3aa041fbf5261f673611caa4c2215f01b4acba55965679c41b308b5"),
    (25001, 50000, "vizier.cfa.harvard.edu.25001-50000", "downloaded",
     "6dfa9b599f254ba8820818c61857efe42171eb1743514a275619162f507dd3ae"),
    (50001, 50538, "vizier.idia.ac.za.50001-75000", "used_verified_prefix",
     "5cf7e533e7196c7ea45479a2414c2115fb70de08002bbd2a1b8796c39021558e"),
]
GZIP = {"method": "gzip", "compresslevel": 6, "mtime": 0}
DACE = "DACE matched southern HARPS spectra"
KEY = ["star", "author_epoch_id", "author_bjd_tick"]


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def verified_file(folder, name, receipts):
    hashes = {line.split(maxsplit=1)[1].strip(): line.split()[0]
              for line in (folder/"SHA256SUMS").read_text().splitlines() if line.strip()}
    path = folder/name
    sha = digest(path)
    if sha != hashes[name]:
        raise ValueError(f"Changed evidence file: {path}")
    receipts.append(dict(file=name, url=f"{folder.name}/{name}", bytes=path.stat().st_size,
                         sha256=sha, status="verified_local_evidence"))
    return path


def corrected_pages(snapshot, manifest):
    data, receipts = [], []
    manifest = pd.read_csv(manifest, keep_default_na=False)
    for first, last, suffix, status, sha in PAGES:
        name = f"J-A+A-683-A125-table4.{suffix}.tsv"
        rows = manifest[manifest.file.eq(name) & manifest.status.eq(status) & manifest.sha256.eq(sha)]
        matches = list(snapshot.rglob(name))
        if len(rows) != 1 or len(matches) != 1:
            raise ValueError(f"Expected one original receipt and one byte stream for {name}")
        path = matches[0]
        if digest(path) != sha or path.stat().st_size != int(rows.iloc[0].bytes):
            raise ValueError(f"Original CDS byte count or SHA-256 mismatch in {name}")
        raw = pd.read_csv(path, comment="#", sep="\t", dtype=str, quoting=csv.QUOTE_NONE)
        raw.columns = raw.columns.str.strip().str.lower()
        records = pd.to_numeric(raw.recno, errors="coerce").dropna()
        np.testing.assert_array_equal(np.sort(records.to_numpy()), np.arange(first, last+1))
        part = field_projection(raw.loc[records.index].copy(), RV_FIELDS)
        part["recno"] = records.astype(int)
        part["page"] = name
        data.append(part.sort_values("recno"))
        receipts.append(rows.iloc[0].to_dict())
    return pd.concat(data, ignore_index=True), receipts


def ticks(bjd):
    # A tick is 8.64 milliseconds. The original S ledger and the author ledger
    # contain the same exposure BJD serialized with potentially different CSV readers.
    return (pd.to_numeric(bjd, errors="coerce") * 1e7).round().astype("Int64")


def author_identity(data, evidence, receipts):
    """Recover the author row that the original run independently identified.

    This bridge contains no author velocities or activity measurements. The
    historical crossmatch includes target identity, BJD and a 2-arcsec sky check.
    """
    bridge = pd.read_csv(verified_file(evidence, "correction-crossmatch.csv.gz", receipts),
                         dtype={"epoch_id": str, "star": str}, float_precision="round_trip")
    old = pd.read_csv(verified_file(evidence, "epoch-ledger.csv.gz", receipts),
                      usecols=["star","bjd","epoch_id","disposition"],
                      dtype={"epoch_id": str,"star": str}, float_precision="round_trip")
    bridge = bridge.set_index("recno", verify_integrity=True).loc[data.recno].reset_index()
    np.testing.assert_array_equal(bridge.recno, data.recno)
    np.testing.assert_array_equal(bridge.star.str.strip().fillna(""), data.star.fillna(""))
    if not np.allclose(pd.to_numeric(bridge.bjd), data.bjd, rtol=0, atol=5e-8, equal_nan=True):
        raise ValueError("Frozen CDS crossmatch does not reference the corrected record times")
    np.testing.assert_array_equal(bridge.epoch_id.str.strip().fillna(""), data.epoch_id.fillna(""))
    good = bridge.match.isin(["target_timestamp_and_bjd","unique_target_bjd"]) & bridge.author_index.notna()
    bridge = bridge[good].copy()
    old_pos = bridge.author_index.to_numpy(dtype=int)
    if len(np.unique(old_pos)) != len(old_pos) or np.any(old_pos < 0) or np.any(old_pos >= len(old)):
        raise ValueError("Author-to-CDS exposure map is nonunique or out of range")
    old["star"] = old.star.str.strip()
    old["author_epoch_id"] = old.epoch_id.str.strip().fillna("")
    old["author_bjd_tick"] = ticks(old.bjd)
    old["unique_key"] = ~old.duplicated(KEY, keep=False)
    picked = old.iloc[old_pos].reset_index(drop=True)
    bridge = bridge.reset_index(drop=True)
    norm = lambda v: v.str.strip().str.upper().str.replace(r"\s+", "", regex=True)
    name_ok = norm(bridge.star).eq(norm(picked.star)).fillna(False)
    tolerance = np.where(bridge.match.eq("unique_target_bjd"), 2.5e-7, 2/86400.)
    same_time = (bridge.bjd - picked.bjd).abs().le(tolerance)
    usable = name_ok & same_time & picked.unique_key & picked.disposition.eq("used")
    joined = pd.concat([bridge[["recno","match"]], picked[KEY]], axis=1)
    return joined[usable].copy()


def attach_sindex(data, bindings, evidence, followup, receipts):
    def read(folder, name):
        return pd.read_csv(verified_file(folder, name, receipts), dtype={"star": str, "epoch_id": str},
                           float_precision="round_trip", low_memory=False)

    published = read(evidence, "ingest-sindex-epochs.csv.gz")
    dace_ledgers = [read(evidence, "sindex-epochs.csv.gz")]
    if followup is not None:
        dace_ledgers.append(read(followup, "all-activity-dace-epochs.csv.gz"))
    dace = pd.concat(dace_ledgers, ignore_index=True)
    instrument = dace.instrument_name.astype("string").str.upper().str.strip()
    mode = dace.ins_mode.astype("string").str.upper().str.strip()
    dace = dace[instrument.isin(["HARPS03","HARPS15","HARPS"]) & mode.isin(["HARPS","HAM"])].copy()
    dace["sindex"] = pd.to_numeric(dace.spectro_smw, errors="coerce")
    dace["sindex_source"] = DACE
    listings = pd.concat([published, dace], ignore_index=True)
    listings = listings[listings.sindex.gt(0) & np.isfinite(listings.sindex) & listings.bjd.notna()].copy()
    listings["star"] = listings.star.str.strip()
    listings["author_epoch_id"] = listings.epoch_id.str.strip().fillna("")
    listings["author_bjd_tick"] = ticks(listings.bjd)
    # Reject same-source disagreements on an exposure. The follow-up also
    # contains many exact copies of the previous DACE ledger.
    variants = listings.groupby(KEY + ["sindex_source"], dropna=False).sindex.transform("nunique")
    conflict = variants.gt(1)
    unique = listings.loc[~conflict].drop_duplicates(KEY+["sindex_source"])
    joined = bindings.merge(unique[KEY+["sindex_source","sindex"]], on=KEY, validate="one_to_many")
    meta = data.loc[data.disposition.eq("used"), ["star","recno","bjd","epoch_id","era"]]
    joined = joined.drop(columns="star").merge(meta, on="recno", validate="many_to_one")
    joined["night"] = np.floor(joined.bjd).astype(int)
    coverage = joined.groupby(["star","sindex_source"],as_index=False).agg(
        matched_epochs=("sindex","size"), nights=("night","nunique"),
        unique_values=("sindex","nunique"))
    coverage["testable"] = coverage.nights.ge(10) & coverage.unique_values.gt(1)
    coverage["published_series"] = coverage.sindex_source.ne(DACE)
    # Prefer a testable published series, otherwise a testable DACE series.
    # A single reduction supplies the whole time series for each target.
    coverage["testable_published"] = coverage.testable & coverage.published_series
    coverage = coverage.sort_values(["star","testable","testable_published","nights","published_series","sindex_source"],
                                    ascending=[True,False,False,False,False,True])
    coverage["selected"] = ~coverage.star.duplicated()
    chosen = coverage[coverage.selected][["star","sindex_source"]]
    matches = joined.merge(chosen, on=["star","sindex_source"], validate="many_to_one")
    lookup = matches.set_index("recno",verify_integrity=True)
    data["sindex"] = data.recno.map(lookup.sindex)
    data["sindex_source"] = data.recno.map(lookup.sindex_source).fillna("unavailable")
    return coverage, matches, int(conflict.sum())


def main(args):
    require_survey_workflow()
    evidence, work, out = args.evidence, args.work, args.out
    work.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    evidence_receipts = []
    raw_manifest = verified_file(evidence, "inputs.csv", evidence_receipts)
    data, page_receipts = corrected_pages(args.snapshot, raw_manifest)
    classify(data, PRIMARY)
    bindings = author_identity(data, evidence, evidence_receipts)
    covered, matched, conflicts = attach_sindex(data, bindings, evidence, args.followup, evidence_receipts)
    holdout = json.loads(verified_file(evidence,"holdout.json",evidence_receipts).read_text())
    data["used"] = data.disposition.eq("used").astype(int)
    data["has_s"] = data.used.eq(1) & data.sindex.notna()
    targets = data.groupby("star", sort=True).agg(
        ra_deg=("ra_deg","median"), dec_deg=("dec_deg","median"), n_epochs=("used","sum"),
        n_input=("recno","size"), first_recno=("recno","min"), last_recno=("recno","max"),
        sindex_epochs=("has_s","sum")).reset_index()
    targets["source"] = SOURCE
    targets.to_csv(out/"targets.csv",index=False)
    used = data.loc[data.used.eq(1)].sort_values(["star","bjd"],kind="stable")
    used[COLUMNS].to_csv(work/"epochs.csv.gz",index=False,compression=GZIP)
    data[["recno","page","star","epoch_id","bjd","era","disposition","sindex_source"]].to_csv(
        out/"epoch-ledger.csv.gz",index=False,compression=GZIP)
    covered.to_csv(out/"sindex-source-coverage.csv",index=False)
    matched.to_csv(out/"subset-sindex-epochs.csv.gz",index=False,compression=GZIP)
    pd.DataFrame(page_receipts + evidence_receipts).to_csv(out/"inputs.csv",index=False)
    audit = dict(source=SOURCE, release_version="2026-04-20 corrected Table 4, verified prefix",
                 cds_corrected_table=True, corrected_table_complete=False,
                 correction_status="published corrected values used directly on a verified contiguous prefix",
                 published_input_rows=252615, input_rows=len(data), first_recno=int(data.recno.min()),
                 last_recno=int(data.recno.max()), unretrieved_rows=252615-len(data),
                 source_pages=len(PAGES), targets=len(targets), used_epochs=len(used),
                 first_bjd=float(used.bjd.min()), last_bjd=float(used.bjd.max()),
                 full_release_last_bjd=float(holdout["latest_primary_bjd"]),
                 dispositions=data.disposition.value_counts().to_dict(),
                 sindex_sources=used.sindex_source.value_counts().to_dict(),
                 sindex_epochs=int(used.sindex.notna().sum()), sindex_stars=int(used.loc[used.sindex.notna(),"star"].nunique()),
                 usable_author_identity_links_for_activity=len(bindings),
                 conflicting_same_source_s_records_excluded=conflicts,
                 sindex_channels_from_followup=args.followup is not None,
                 warnings=["The prefix is ordered by CDS record number and is not a random sample of the full archive. "
                           "An author-release exposure link supplies auxiliary S-index identity only."])
    if not len(used):
        raise ValueError("No usable corrected Table 4 velocities")
    (out/"ingest.json").write_text(json.dumps(audit,indent=2)+"\n")
    print(json.dumps(audit),flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for arg in ["snapshot","evidence","work","out"]:
        parser.add_argument("--"+arg, required=True, type=Path)
    parser.add_argument("--followup",type=Path)
    main(parser.parse_args())
