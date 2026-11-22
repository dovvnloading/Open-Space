"""Validate shard coverage before summarizing blind fits, controls and injections."""
import argparse
from collections import defaultdict
import difflib
import hashlib
import json
import os
import re
from pathlib import Path

import pandas as pd
from runtime import require_survey_workflow


# The blind search does not import this module. These approximate published periods
# are used for a retrospective test of recovery, including all misses.
CONTROLS = {
    "HD4308": {"b": 15.56},
    "HD40307": {"b": 4.3123, "c": 9.6184, "d": 20.432},
    "HD69830": {"b": 8.667, "c": 31.56, "d": 197.0},
}


def read(path):
    return pd.read_csv(path, dtype={"star": str}, keep_default_na=False, low_memory=False)


def combine(out, shards, count):
    directories = [shards / f"harps-shard-{i}" for i in range(count)]
    targets = read(out / "targets.csv").set_index("star", verify_integrity=True)
    combined = {}
    for name in ("stars", "signals", "injections"):
        pieces = [read(shard / f"{name}.csv") for shard in directories]
        if any(tuple(p.columns) != tuple(pieces[0].columns) for p in pieces):
            raise ValueError(f"Column mismatch between {name} shards")
        combined[name] = pd.concat(pieces, ignore_index=True)
    stars = combined["stars"].set_index("star", verify_integrity=True)
    if set(stars.index) != set(targets.index):
        raise ValueError("Shard union differs from the input target inventory")
    if not (stars.n_epochs.astype(int) == targets.n_epochs.reindex(stars.index).astype(int)).all():
        raise ValueError("Per-target RV counts differ from ingestion")
    if not stars.eligibility.isin(["searched", "insufficient_nights", "insufficient_span"]).all():
        raise ValueError("Unresolved per-star error or unexpected eligibility")
    searched = stars[stars.eligibility.eq("searched")]
    signals = combined["signals"]
    keys = ["star", "signal_index"]
    if signals.duplicated(keys).any() or set(signals.star) != set(searched.index):
        raise ValueError("An eligible star lacks a fit or a signal appears twice")
    if not (signals.groupby("star").size().reindex(searched.index) == searched.n_scans.astype(int)).all():
        raise ValueError("Fits per target do not match the search ledger")
    injections = combined["injections"]
    if set(injections.star) != set(searched.index):
        raise ValueError("Injection coverage and period-search coverage differ")
    if injections.duplicated(["star", "period_d", "k_ms", "eccentricity"]).any():
        raise ValueError("Duplicate injection cell")
    if (injections.trials.le(0) | injections.recovered.lt(0) |
            injections.recovered.gt(injections.trials)).any():
        raise ValueError("Impossible injection counts")
    for name, frame in combined.items():
        frame.to_csv(out / f"{name}.csv", index=False)


def control_tables(stars, signals, vetting):
    """A prespecified planet can be recovered once across independently fitted archive series."""
    def hd(name):
        m = re.fullmatch(r"HD[\s_-]*0*(\d+)", str(name).strip(), re.I)
        return "HD" + str(int(m[1])) if m else ""

    hosts = defaultdict(set)
    object_ids = defaultdict(set)
    for star in stars.star:
        if hd(star) in CONTROLS:
            hosts[hd(star)].add(star)
    cols = ["star","identity","aliases"]
    if "simbad_oid" in vetting.columns:
        cols.append("simbad_oid")
    for ref in vetting[cols].drop_duplicates().itertuples(index=False):
        oid = str(getattr(ref, "simbad_oid", "")).strip() or str(ref.identity).strip()
        if oid and ref.identity:
            object_ids[ref.star].add(oid)
            for name in [ref.identity,*str(ref.aliases).split("|")]:
                if hd(name) in CONTROLS:
                    hosts[hd(name)].add(ref.star)

    inventory = stars.set_index("star",verify_integrity=True)
    systems, series = [], []
    for host, periods in CONTROLS.items():
        names = sorted(hosts[host])
        identified = set().union(*(object_ids[s] for s in names))
        # A directly named HD series may have no fitted peaks and no SIMBAD row.
        # Other series must be linked to that unqualified HD identifier; divergent
        # OIDs make the entire assignment ambiguous, regardless of measured periods.
        conflict = len(identified)>1
        for companion, period in periods.items():
            base = dict(host=host,companion=companion,reference_period_d=period)
            parts = []
            for star in names:
                r = dict(**base,star=star)
                coverage = inventory.loc[star]
                if conflict:
                    r["outcome"] = "ambiguous target identity"
                elif coverage.eligibility != "searched":
                    r["outcome"] = coverage.eligibility
                else:
                    hits = signals[signals.star.eq(star)].copy()
                    hits["distance"] = (1/hits.period_d - 1/period).abs()
                    hits = hits[hits.distance <= min(.05/period,1/float(coverage.baseline_d))]
                    if hits.empty:
                        r["outcome"] = "period missed"
                    else:
                        hits["significant"] = hits.fap.le(.01) & hits.delta_bic.ge(10)
                        hit = hits.sort_values(["significant","distance","signal_index"],
                                               ascending=[False,True,True]).iloc[0]
                        r.update(signal_index=int(hit.signal_index), fit_period_d=hit.period_d,
                                 fit_k_ms=hit.k_ms,fap=hit.fap,rv_status=hit.rv_status,
                                 flags=hit["flags"],outcome=("period recovered above threshold"
                                 if hit.significant else "period recovered below threshold"))
                parts.append(r)
            series.extend(parts)
            system = dict(**base,matching_targets=" | ".join(names))
            priorities = ["period recovered above threshold","period recovered below threshold",
                          "period missed","ambiguous target identity"]
            unique = {r["outcome"] for r in parts}
            system["outcome"] = next((status for status in priorities if status in unique),
                                     " | ".join(sorted(unique)) or "host not identified")
            # Stable representative values for readers of the original seven-row
            # table. The per-series ledger preserves all other fits and misses.
            representative = next((r for r in parts if r["outcome"]==system["outcome"]),{})
            system.update({k:v for k,v in representative.items() if k not in base})
            systems.append(system)
    fields = ["host","companion","reference_period_d","star","outcome","signal_index",
              "fit_period_d","fit_k_ms","fap","rv_status","flags"]
    return pd.DataFrame(systems), pd.DataFrame(series).reindex(columns=fields)


def review_tables(out, signals, vetting):
    # Use the complete fit key, including earlier peaks on multiplanet targets.
    keep = ["star", "signal_index", "delta_bic", "gaussian_trials_screen", "flags"]
    joined = vetting.merge(signals[keep], on=["star", "signal_index"], validate="one_to_one")
    fields = ["star", "signal_index", "identity", "period_d", "k_ms", "ecc", "fap",
              "minimum_mass_mjup", "dynamical_mass_mjup", "joint_status", "delta_bic",
              "gaussian_trials_screen", "rv_status", "external_status", "flags", "reasons",
              "prior_signal", "alias_matches", "rotation_check", "binary_check", "nasa_check",
              "literature_check", "bibliography_matches"]
    ranked = joined[fields].sort_values(["fap", "delta_bic", "star", "signal_index"],
                                        ascending=[True, False, True, True])
    survivors = ranked[ranked.rv_status.eq("pass")]
    survivors.to_csv(out / "rv-survivors.csv", index=False)
    inconclusive = ranked[(ranked.fap <= .01) & (ranked.delta_bic >= 10) &
                          ranked.rv_status.eq("inconclusive") &
                          ranked.external_status.eq("inconclusive")]
    inconclusive.to_csv(out / "inconclusive-periods.csv", index=False)
    papers = read(out / "literature.csv.gz")
    papers[papers.star.isin(survivors.star)].sort_values(["star", "bibcode"]).to_csv(
        out / "survivor-literature.csv", index=False)
    return len(survivors), len(inconclusive)


def query_coverage(out):
    receipt = read(out / "catalog-inputs.csv")
    rows = []
    for source, data in receipt.groupby("source", sort=True):
        ok = data.status.eq("ok")
        failures = data.status.eq("failed")
        rows.append(dict(source=source, queries_planned=len(data), successful=int(ok.sum()),
                         failed=int(failures.sum()), skipped=int((~ok & ~failures).sum()),
                         successful_rows=int(pd.to_numeric(data.loc[ok, "rows"]).sum())))
    pd.DataFrame(rows, columns=["source", "queries_planned", "successful", "failed",
                                "skipped", "successful_rows"]).to_csv(out / "source-coverage.csv", index=False)
    return int(receipt.status.eq("ok").sum()), len(receipt)


def activity_coverage(out, signals):
    from search import INDICATORS, REQUIRED
    detail = signals.activity_json.map(json.loads)
    flags = signals["flags"].map(lambda value: set(value.split("|")))
    rows = []
    all_required = pd.Series(True, index=signals.index)
    for name in INDICATORS:
        tested = detail.map(lambda row: row[name]["status"] == "tested")
        flagged = flags.map(lambda row: bool({"activity_" + name, "rotation_proxy_" + name} & row))
        if name in REQUIRED:
            all_required &= tested
        rows.append(dict(indicator=name, fits_tested=int(tested.sum()),
                         target_series_tested=signals.loc[tested, "star"].nunique(),
                         fits_missing_or_insufficient=int((~tested).sum()),
                         fits_flagged=int(flagged.sum())))
    pd.DataFrame(rows).to_csv(out / "activity-coverage.csv", index=False)
    return int(all_required.sum()), signals.loc[all_required, "star"].nunique()


def catalog_patch(out, efficiency):
    """Prepare a reviewable catalog edit only after measuring injection recovery."""
    source = Path("README.md").read_text()
    heading = "### Planets around other stars\n\n"
    target = "results/harps-rvbank-survey/README.md"
    if target in source:
        return  # A later run must not duplicate an existing catalog entry.
    lines = source.splitlines(keepends=True)
    if len(lines) < 9 or source.count(heading) != 1:
        raise ValueError("Review the catalog structure before preparing its update")
    line = lines[8]
    for pattern in [r"(?<=has )\d+(?= entries)", r"\d+(?= new measurements)"]:
        line, replacements = re.subn(pattern, lambda m: str(int(m.group()) + 1), line)
        if replacements != 1:
            raise ValueError("Cannot identify catalog counts on line 9")
    lines[8] = line
    cell = efficiency[(efficiency.period_d == 30) & (efficiency.k_ms == 10) &
                      (efficiency.eccentricity == 0)].iloc[0]
    bullet = (f"- **HARPS sensitivity to 30-day companions.** The blind period search "
              f"recovered {int(cell.recoveries):,} of {int(cell.trials):,} circular "
              f"10 m/s injections across {int(cell.stars):,} stars. "
              f"*Measured* · [Details]({target})\n")
    revised = "".join(lines).replace(heading, heading + bullet, 1)
    diff = difflib.unified_diff(source.splitlines(keepends=True),
                                revised.splitlines(keepends=True),
                                fromfile="a/README.md", tofile="b/README.md")
    (out / "catalog.patch").write_text("".join(diff))


def summarize(out):
    stars, signals, injections, vetting = [read(out / f"{n}.csv") for n in
                                            ("stars", "signals", "injections", "vetting")]
    keys = ["star", "signal_index"]
    if (vetting.duplicated(keys).any() or signals.duplicated(keys).any() or
            set(map(tuple, vetting[keys].to_numpy())) != set(map(tuple, signals[keys].to_numpy()))):
        raise ValueError("External vetting must have one disposition per fitted peak")
    searched = stars[stars.eligibility.eq("searched")]
    if searched.empty:
        raise ValueError("No eligible targets. Inspect coverage before interpreting the run.")
    controls, control_series = control_tables(stars, signals, vetting)
    control_series.to_csv(out / "published-control-series.csv", index=False)
    rv_survivors, unresolved = review_tables(out, signals, vetting)
    source_ok, source_total = query_coverage(out)
    controls.to_csv(out / "published-controls.csv", index=False)
    # Fractions describe the initial detection stage at tested P, K and e on eligible
    # stars whose baseline admits P. They do not estimate population occurrence.
    efficiency = injections.groupby(["period_d", "k_ms", "eccentricity"], as_index=False).agg(
        stars=("star", "nunique"), trials=("trials", "sum"), recoveries=("recovered", "sum"))
    efficiency["recovery_fraction"] = efficiency.recoveries / efficiency.trials
    efficiency.to_csv(out / "completeness-grid.csv", index=False)
    inputs = read(out / "inputs.csv")
    downloaded = inputs[inputs.status.eq("downloaded")]
    if downloaded.empty or not downloaded.sha256.map(lambda s: bool(re.fullmatch(r"[0-9a-f]{64}", s))).all():
        raise ValueError("Missing spectroscopic input receipt")
    run = os.environ["GITHUB_RUN_ID"]
    repo = os.environ["GITHUB_REPOSITORY"]
    revision = (out / "revision.txt").read_text().strip()
    ingest = json.loads((out / "ingest.json").read_text())
    dispositions = {str(k): int(v) for k, v in vetting.external_status.value_counts().items()}
    counts = dict(
        source=ingest["source"], code_revision=revision, actions_run=run,
        catalog_targets=len(stars), spectra=int(stars.n_epochs.sum()),
        searched_targets=len(searched), searched_spectra=int(searched.n_epochs.sum()),
        fits=len(signals), injection_trials=int(injections.trials.sum()),
        rv_survivors=rv_survivors, significant_inconclusive=unresolved,
        queued_for_literature_review=int(vetting.external_status.eq("requires literature review").sum()),
        external_queries_successful=source_ok, external_queries_planned=source_total,
        status_counts=dispositions, published_controls=len(controls),
        controls_above_threshold=int(controls.outcome.eq("period recovered above threshold").sum()),
    )
    accelerations = read(out / "accelerations.csv")
    fitted_mass = pd.to_numeric(vetting.dynamical_mass_mjup, errors="coerce").notna()
    vetting.loc[fitted_mass].to_csv(out / "dynamical-masses.csv", index=False)
    counts.update(
        trend_targets=len(accelerations),
        trends_with_hgca=int(accelerations.astrometry_check.eq("matched").sum()),
        trends_with_pma_p_below_001=int(pd.to_numeric(accelerations.pma_p, errors="coerce").lt(.001).sum()),
        adopted_conditional_dynamical_masses=int(pd.to_numeric(vetting.dynamical_mass_mjup, errors="coerce").notna().sum()),
        rv_pass_by_external_status={k:int(v) for k,v in vetting.loc[vetting.rv_status.eq("pass"),
                                                              "external_status"].value_counts().items()},
        target_eligibility={k:int(v) for k,v in stars.eligibility.value_counts().items()},
    )
    counts["all_six_activity_fits"], counts["all_six_activity_stars"] = activity_coverage(out, signals)
    (out / "survey.json").write_text(json.dumps(counts, indent=2) + "\n")
    body = f"""# HARPS radial-velocity companion survey

**The search fitted {len(signals):,} periodic signals in {len(searched):,} of {len(stars):,} catalog targets.** The targets that met the cadence cuts contributed {counts['searched_spectra']:,} usable spectra. A fitted peak can arise from noise, activity, binarity, sampling or a known companion. Individual novelty claims require review of the source papers and object identities.

The input was {ingest['source']}. See [input checksums](inputs.csv), [exposure dispositions](epoch-ledger.csv.gz) and [target eligibility](stars.csv). The [RV fit table](signals.csv) records periods, semiamplitudes, false-alarm estimates and diagnostic checks. The [external vetting table](vetting.csv) records companion-mass constraints and object-level checks, including unresolved evidence. The [acceleration table](accelerations.csv) records long-term trends.

{rv_survivors:,} fits passed all internal RV tests. The [RV survivors](rv-survivors.csv) include their catalog and literature dispositions. The [associated bibliography](survivor-literature.csv) contains indexed references. Another {unresolved:,} fits met the statistical thresholds and retained unresolved checks in both the RV and external ledgers. See [inconclusive periods](inconclusive-periods.csv). Of {source_total:,} planned external queries, {source_ok:,} returned usable responses. [Coverage by source](source-coverage.csv) · [Query receipts](catalog-inputs.csv)

The initial period search recovered {counts['controls_above_threshold']} of {len(controls)} prespecified published control periods above the adopted FAP and model-improvement thresholds. The denominator includes every prespecified control, including misses and targets outside coverage. Each planet counts once across its [independently fitted archive series](published-control-series.csv). See [control recovery](published-controls.csv).

The [injection grid](completeness-grid.csv) summarizes {counts['injection_trials']:,} injected orbits in observed residuals. Its denominator varies with orbital period and available baseline. It measures recovery of the blind periodogram peak before external vetting. The [individual injection cells](injections.csv) retain both recoveries and misses.

## Signal dispositions

"""
    body += "".join(f"- {status}: {count:,}\n" for status, count in sorted(dispositions.items()))
    body += f"""
## Reproduce the run

Use the source at commit `{revision}`, the [survey workflow](../../../.github/workflows/harps-rvbank-survey.yml) and [installed software versions](packages.txt). The workflow records its commands, including shard count and selection criteria. Compare the exact archive and catalog input checksums between runs. [Actions run](https://github.com/{repo}/actions/runs/{run})
"""
    (out / "README.md").write_text(body)
    catalog_patch(out, efficiency)
    files = sorted(p for p in out.iterdir() if p.is_file() and p.name != "SHA256SUMS")
    with (out / "SHA256SUMS").open("w") as output:
        for path in files:
            output.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.name}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="step", required=True)
    join = sub.add_parser("combine")
    join.add_argument("--shards", type=Path, required=True)
    join.add_argument("--count", type=int, required=True)
    for cmd in (join, sub.add_parser("summarize")):
        cmd.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    require_survey_workflow()
    if args.step == "combine":
        combine(args.out, args.shards, args.count)
    else:
        summarize(args.out)


if __name__ == "__main__":
    main()
