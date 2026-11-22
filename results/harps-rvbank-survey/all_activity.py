"""Evaluate Mount Wilson S on every frozen RV fit using matched public HARPS exposures."""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from threading import local
import time

import numpy as np
import pandas as pd
import requests

from enrich import COLUMNS, fetch, match_exposures
from runtime import require_survey_workflow
from search import INDICATORS, TAU, activity, nightly, nuisance, orbit_basis


S_FLAGS = {"missing_sindex", "activity_sindex", "rotation_proxy_sindex"}
DACE_SOURCE = "DACE matched exposures"
THREAD = local()


def session():
    if not hasattr(THREAD, "client"):
        client = requests.Session()
        client.headers.update({"Accept": "application/json", "User-Agent": "Open-Space-HARPS-survey/1.0"})
        THREAD.client = client
    return THREAD.client


def canonical_names(vetting):
    """Only the single SIMBAD main identifier from a successful identity check qualifies."""
    aliases = {}
    for star, rows in vetting.groupby("star", sort=False):
        unique = rows[["simbad_oid", "simbad_identity_route", "identity"]].drop_duplicates()
        if len(unique) != 1:
            continue
        oid, route, name = unique.iloc[0]
        if oid and route in {"name", "position"} and name.strip():
            # Preserve the actual main identifier, including the stellar component.
            aliases[star] = " ".join(name.split())
    return aliases


def query(star, name, route, frame, cache):
    receipt = {"star": star, "queried_name": name, "query_route": route, "provenance": "follow-up request"}
    pairs = pd.DataFrame(columns=COLUMNS)
    try:
        with TemporaryDirectory(dir=cache) as scratch:
            data = fetch(session(), name, receipt, Path(scratch))
            if data is not None:
                pairs = match_exposures(star, frame, data, receipt)
                receipt["status"] = "matched exposures" if len(pairs) else "no matching exposure"
        receipt.update(matched_epochs=len(pairs), nights=int(np.unique(np.floor(pairs.bjd.astype(float))).size))
    except (requests.RequestException, ValueError, KeyError, TypeError, OSError, IndexError) as exc:
        receipt.update(status="query or exposure matching failed", error=str(exc)[:500])
    time.sleep(.25)
    return receipt, pairs


def coverage(star, frame, original, archived, alias, cache):
    receipts = []
    pairs = archived
    if original is not None:
        record = dict(original, star=star, queried_name=star, query_route="RVBank name",
                      provenance="revised/sindex-receipts.csv")
        receipts.append(record)
    # Earlier valid queries and definitive not-found responses have receipts.
    # Transport failures can be retried once by this follow-up.
    complete = original is not None and (
        str(original.get("matched_epochs", "")).strip() != "" or str(original.get("http_status", "")) == "404")
    if not complete and pairs.empty:
        record, pairs = query(star, star, "RVBank name", frame, cache)
        receipts.append(record)
    if pairs.empty and alias and " ".join(star.split()).casefold() != alias.casefold():
        record, pairs = query(star, alias, "verified SIMBAD main identifier", frame, cache)
        receipts.append(record)
    return receipts, pairs


def flags_of(value):
    return set(filter(None, str(value).split("|")))


def status_of(flags):
    if {"bootstrap_fap", "weak_orbit"} & flags or any(
            f.startswith(("activity_", "rotation_proxy_")) for f in flags):
        return "rejected"
    return "inconclusive" if flags else "pass"


def dace_values(frame, pairs):
    """Attach a previously proven exposure pair to the corresponding normalized input row."""
    if pairs.empty:
        return np.full(len(frame), np.nan)
    if pairs.bjd.duplicated().any() or not pairs.star.eq(frame.star.iloc[0]).all():
        raise ValueError("DACE evidence contains duplicate BJDs or a target mismatch")
    # Old DACE receipts serialize the BJD used by enrich.py. Round-trip input
    # on the receipt side preserves the exact double read from the epochs file.
    table = pairs.set_index("bjd", verify_integrity=True)
    table_ids = table.epoch_id.fillna("").astype(str)
    actual_ids = frame.set_index("bjd", verify_integrity=True).epoch_id.fillna("").astype(str)
    if not table_ids.index.isin(actual_ids.index).all() or not np.array_equal(
            table_ids.to_numpy(), actual_ids.reindex(table_ids.index).to_numpy()):
        raise ValueError(f"Frozen DACE pairs no longer identify the same exposures on {frame.star.iloc[0]}")
    return frame.bjd.map(table.spectro_smw).to_numpy(dtype=float)


def s_only_nights(frame):
    # The original activity() corrects its significance for seven channels.
    # Keep the same channel list while leaving the other six unpopulated.
    data = frame[["bjd", "err", "rv", "era", "sindex"]].copy()
    return nightly(data)


def replay(fits, frame, pairs):
    measured = s_only_nights(frame)
    dace = s_only_nights(frame.assign(sindex=dace_values(frame, pairs))) if len(pairs) else None
    n_dace = int(np.isfinite(dace.sindex).sum()) if dace is not None else 0
    if dace is not None and (n_dace < 10 or dace.sindex.nunique() <= 1):
        dace = None
    n_usable = n_dace if dace is not None else 0
    t, y, err = (measured[k].to_numpy() for k in ("bjd", "rv", "err"))
    span = float(np.ptp(t))
    design = nuisance(t, measured.era.to_numpy())
    n_base = design.shape[1]
    ordered = fits.sort_values("signal_index")
    np.testing.assert_array_equal(ordered.signal_index, np.arange(1, len(ordered)+1))
    result = []
    for fit in ordered.itertuples(index=False):
        # With tperi = tmin - (M0 mod 2pi)*P/2pi this recreates the two
        # original Keplerian quadratures. search() refits both coefficients
        # for every preceding orbit jointly with the current orbit and trend.
        m0 = TAU * (t.min() - fit.tperi_bjd) / fit.period_d
        basis = orbit_basis(t-t.min(), 1/fit.period_d, fit.ecc, m0)
        extended = np.column_stack([design, basis])
        coef = np.linalg.lstsq(extended / err[:, None], y / err, rcond=None)[0]
        # Relative tolerance covers serializing tperi at Julian dates near
        # 2.45 million, including nearly circular fits and very small K.
        np.testing.assert_allclose(np.hypot(*coef[-2:]), fit.k_ms, rtol=1e-4, atol=1e-5)
        current = y - design[:, n_base:] @ coef[n_base:-2]
        design = extended
        original = flags_of(fit.flags)
        prior_s = json.loads(fit.activity_json)["sindex"]
        prior_tested = prior_s["status"] == "tested"
        old_s_flags = original & S_FLAGS
        new_info, new_s_flags = {}, set()
        if dace is not None:
            values, flaglist = activity(dace, current, 1/fit.period_d, span)
            new_info, new_s_flags = values["sindex"], set(flaglist) & S_FLAGS
        if fit.sindex_source == DACE_SOURCE and prior_tested:
            if new_info.get("status") != "tested" or old_s_flags != new_s_flags or prior_s["n"] != new_info["n"]:
                raise ValueError(f"Frozen DACE S diagnostic did not reproduce for {fit.star} #{fit.signal_index}")
        elif prior_tested and fit.rv_status == "pass":
            info, flags = activity(measured, current, 1/fit.period_d, span)
            if info["sindex"]["n"] != prior_s["n"] or old_s_flags != (set(flags) & S_FLAGS):
                raise ValueError(f"Frozen input S diagnostic did not reproduce for {fit.star} #{fit.signal_index}")

        accepted = new_info.get("status") == "tested"
        final_s = old_s_flags.copy()
        if accepted:
            final_s.discard("missing_sindex")
            final_s |= new_s_flags
        final_flags = (original - S_FLAGS) | final_s
        new_status = status_of(final_flags)
        if status_of(original) != fit.rv_status:
            raise ValueError(f"Original disposition no longer reproducible for {fit.star} #{fit.signal_index}")
        if not accepted and new_status != fit.rv_status:
            raise ValueError("Prior disposition changed without new S-index evidence")
        if accepted and fit.sindex_source == DACE_SOURCE and fit.rv_status == "pass" and new_status != "pass":
            raise ValueError("Previously passing identical DACE data changed disposition")
        result.append(dict(star=fit.star, signal_index=fit.signal_index, period_d=fit.period_d,
            k_ms=fit.k_ms, fap=fit.fap, delta_bic=fit.delta_bic,
            original_rv_status=fit.rv_status, new_rv_status=new_status, original_flags=fit.flags,
            new_flags="|".join(sorted(final_flags)), unchanged_non_s_flags="|".join(sorted(original-S_FLAGS)),
            prior_s_source=fit.sindex_source, prior_s_check=prior_s["status"],
            prior_s_flags="|".join(sorted(old_s_flags)), prior_s_json=json.dumps(prior_s),
            dace_matched_epochs=len(pairs), dace_matched_nights=n_dace,
            dace_check=new_info.get("status", "missing_or_insufficient"),
            dace_s_flags="|".join(sorted(new_s_flags)), dace_s_json=json.dumps(new_info),
            final_s_check="tested" if prior_tested or accepted else "missing_or_insufficient",
            final_s_flags="|".join(sorted(final_s))))
    return result, n_usable


def main(primary, evidence, out):
    require_survey_workflow()
    out.mkdir(parents=True, exist_ok=True)
    signals = pd.read_csv(evidence/"signals.csv", dtype={"star": str}, keep_default_na=False)
    vetting = pd.read_csv(evidence/"vetting.csv", dtype=str, keep_default_na=False)
    refs = canonical_names(vetting)
    main_epochs = pd.read_csv(primary, dtype={"star": str, "epoch_id": str}, low_memory=False)
    groups = main_epochs.groupby("star", sort=False)
    fits = signals.groupby("star", sort=False)
    names = sorted(fits.groups)
    previous = pd.read_csv(evidence/"sindex-receipts.csv", dtype=str, keep_default_na=False)
    old_records = previous.set_index("star", verify_integrity=True).to_dict("index")
    old_pairs = pd.read_csv(evidence/"sindex-epochs.csv.gz", dtype={"star":str, "epoch_id":str},
                            float_precision="round_trip", low_memory=False)
    pair_groups = {s:g for s,g in old_pairs.groupby("star", sort=False)}
    empty = pd.DataFrame(columns=COLUMNS)
    checks, ledgers, receipts, targets = [], [], [], []
    with TemporaryDirectory(dir=primary.parent) as cache, ThreadPoolExecutor(max_workers=3) as pool:
        def lookup(star):
            result = coverage(star, groups.get_group(star), old_records.get(star),
                              pair_groups.get(star, empty), refs.get(star), cache)
            return star, result
        for count, (star, (attempts, pairs)) in enumerate(pool.map(lookup, names), 1):
            rows, n_usable = replay(fits.get_group(star), groups.get_group(star), pairs)
            checks.extend(rows)
            receipts.extend(attempts)
            if len(pairs):
                ledgers.append(pairs[COLUMNS].copy())
            targets.append(dict(star=star, fitted_signals=len(rows), queried_names=" | ".join(
                r["queried_name"] for r in attempts), query_statuses=" | ".join(
                r.get("status", "unknown") for r in attempts), matched_epochs=len(pairs),
                accepted_series_nights=n_usable))
            if count % 100 == 0 or count == len(names):
                print(f"S-index review {count}/{len(names)} targets", flush=True)
    audit = pd.DataFrame(checks).sort_values(["star","signal_index"])
    keys = ["star","signal_index"]
    if audit.duplicated(keys).any() or len(audit) != len(signals):
        raise ValueError("Incomplete follow-up fit ledger")
    if not audit[keys].reset_index(drop=True).equals(signals[keys].sort_values(keys).reset_index(drop=True)):
        raise ValueError("Fit keys differ from the recorded search")
    old_pass = audit.original_rv_status.eq("pass")
    new_pass = audit.new_rv_status.eq("pass")
    audit.to_csv(out/"all-activity.csv", index=False)
    audit[new_pass & ~old_pass].to_csv(out/"all-activity-new-passes.csv", index=False)
    pd.DataFrame(receipts).to_csv(out/"all-activity-receipts.csv", index=False)
    pd.DataFrame(targets).to_csv(out/"all-activity-targets.csv", index=False)
    pd.concat([empty,*ledgers], ignore_index=True).to_csv(out/"all-activity-dace-epochs.csv.gz",
            index=False, compression={"method":"gzip","mtime":0})
    summary = dict(fits=len(audit), fitted_targets=len(names), original_rv_passes=int(old_pass.sum()),
        followup_rv_passes=int(new_pass.sum()), newly_passing=int((~old_pass & new_pass).sum()),
        no_longer_passing=int((old_pass & ~new_pass).sum()), initial_tested_fits=int(audit.prior_s_check.eq("tested").sum()),
        dace_tested_fits=int(audit.dace_check.eq("tested").sum()), final_tested_fits=int(audit.final_s_check.eq("tested").sum()),
        fit_dispositions=dict(Counter(audit.new_rv_status)), requests=dict(Counter(r["provenance"] for r in receipts)),
        request_statuses=dict(Counter(str(r.get("status", "unknown")) for r in receipts)),
        targets_with_usable_dace=sum(int(r["accepted_series_nights"] > 0) for r in targets),
        matched_epochs=sum(len(df) for df in ledgers), calibrated_activity_channels=len(INDICATORS))
    (out/"all-activity-summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary", type=Path, required=True)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    main(args.primary, args.evidence, args.out)
