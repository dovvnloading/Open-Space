"""Validate coverage and write the result and its catalog entry from workflow evidence."""

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
from collections import Counter
from statistics import median
from adjudicate import audit

FOLDER = "results/nanograv15-transients"
SECTION = "### Neutron stars, black holes and cosmic explosions"
LABELS = ("Measured", "Candidate", "Prediction", "Correction", "Null result")


def read(path):
    return json.loads(path.read_text())


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n")


def table(path, rows, fields):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value
                             for key, value in row.items()})


def index(rows, key):
    values = {row[key]: row for row in rows}
    if len(values) != len(rows):
        raise ValueError(f"duplicate {key} in evidence")
    return values


def inventory(evidence):
    manifest = read(evidence / "inputs.json")
    sources = manifest["sources"]
    targets = {r["pulsar"] for r in sources if r["dataset"] == "nb"}
    if len(targets) != 68:
        raise ValueError("incomplete narrowband input inventory")
    if targets != {r["pulsar"] for r in sources if r["dataset"] == "wb"}:
        raise ValueError("wideband pulsar inventory differs from narrowband")
    coverage = {}
    for ds in ("nb", "wb", "epta"):
        expected = {row["pulsar"] for row in sources if row["dataset"] == ds} & targets
        rows = [r for part in sorted(evidence.glob(f"coverage-{ds}-*.json")) for r in read(part)]
        by_name = index(rows, "pulsar")
        if any(r["dataset"] != ds or r["status"] not in ("ok", "inconclusive") for r in rows):
            raise ValueError(f"invalid {ds} coverage entry")
        if by_name.keys() != expected:
            raise ValueError(f"{ds} coverage ledger does not match the selected archive targets")
        coverage[ds] = by_name
    return manifest, coverage


def archive(evidence, dest):
    # Deposit only generated evidence. Raw arrival times and dense derived matrices stay in Actions.
    dest.mkdir(parents=True, exist_ok=False)
    for file in sorted(evidence.rglob("*")):
        if file.is_file():
            target = dest / file.relative_to(evidence)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(file, target)


def catalog(root, label, bullet):
    path = root / "README.md"
    lines = path.read_text().splitlines(keepends=True)
    link = FOLDER + "/README.md"
    matches = [i for i, line in enumerate(lines) if link in line]
    if len(matches) > 1:
        raise ValueError("multiple existing catalog links to this study")
    entry = f"- **NANOGrav timing transients.** {bullet} *{label}* · [Details]({link})\n"
    if matches:
        lines[matches[0]] = entry
    else:
        section = [i for i, line in enumerate(lines) if line.rstrip() == SECTION]
        if len(section) != 1:
            raise ValueError("cannot identify the neutron-star section")
        lines.insert(section[0] + 2, entry)
    pattern = re.compile(r"(?<!\*)\*(" + "|".join(map(re.escape, LABELS)) + r")\*(?!\*)")
    counts = dict.fromkeys(LABELS, 0)
    for line in lines:
        found = pattern.findall(line) if re.match(r"\s*[-*+]\s+", line) else []
        if found:
            counts[found[-1]] += 1
    phrases = {"entries": sum(counts.values()), "new measurements": counts["Measured"],
               "candidate discoveries": counts["Candidate"], "forecasts": counts["Prediction"],
               "corrections": counts["Correction"], "carefully recorded null": counts["Null result"]}
    if "The catalog below has" not in lines[8]:
        raise ValueError("catalog summary moved from line 9")
    for phrase, count in phrases.items():
        lines[8], n = re.subn(r"\b[\d,]+ (?=" + re.escape(phrase) + r"\b)", f"{count} ", lines[8])
        if n != 1:
            raise ValueError(f"cannot reconcile summary count for {phrase}")
    path.write_text("".join(lines))


def render(root, evidence):
    folder = root / FOLDER
    manifest, cov = inventory(evidence)
    singles = index([read(p) for p in sorted((evidence / "single").glob("*.json"))], "pulsar")
    if singles.keys() != {s for s, row in cov["nb"].items() if row["status"] == "ok"}:
        raise ValueError("single-pulsar outputs differ from the successful narrowband reductions")
    ok = [r for r in singles.values() if r["status"] == "ok"]
    array = read(evidence / "array" / "array_summary.json")
    if not ok or not array.get("n_trials") or array.get("status") != "analyzed":
        raise ValueError("no usable individual or common-mode search; inspect the workflow evidence")
    survivor_rows = [dict(row, pulsar=p["pulsar"], search="individual") for p in ok for row in p["survivors"]]
    survivor_rows.extend(dict(row, search="array") for row in array["survivors"])
    primary_count = sum(len(r["primary_crossings"]) for r in ok)
    limits = [dict(row, pulsar=p["pulsar"]) for p in ok for row in p["sensitivity"]]
    reported_limits = [r for r in limits if r.get("reported_upper_limit_us") is not None]
    coverage = [dict(pulsar=name, nb=cov["nb"][name]["status"], wb=cov["wb"][name]["status"],
                     epta=cov["epta"].get(name, {}).get("status", "outside_release"),
                     epta_selection=cov["epta"].get(name, {}).get("release_selection", ""),
                     individual=singles.get(name, {}).get("status", "no_reduction"))
                for name in sorted(cov["nb"])]
    n_wb = sum(r["wb"] == "ok" for r in coverage)
    n_epta = sum(r["epta"] == "ok" for r in coverage)
    coherent_limit = any(row.get("first_tested_amplitude_90pct") is not None
                         for row in array["sensitivity"])
    complete_null = (not primary_count and not array["threshold_clusters"]
                     and len(ok) == n_wb == 68 and reported_limits and coherent_limit
                     and n_epta == len(cov["epta"]) and n_epta)
    label = "Null result" if complete_null else "Measured"
    dest = folder / "evidence"
    if dest.exists():
        # Replace only the study's own previously generated evidence on a replay.
        if not (dest / "run.json").is_file():
            raise ValueError("existing evidence directory has no prior run receipt")
        shutil.rmtree(dest)
    archive(evidence, dest)
    run = {key: os.environ.get(key) for key in
           ("GITHUB_SHA", "GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT", "GITHUB_REPOSITORY", "GITHUB_REF",
            "PRIMARY_RUN", "PRIMARY_SHA", "NG_RUN")}
    write_json(dest / "run.json", run)
    table(dest / "coverage.csv", coverage, ("pulsar", "nb", "wb", "epta", "epta_selection", "individual"))
    write_json(dest / "surviving-fits.json", survivor_rows)
    disposition = []
    for p in sorted(ok, key=lambda item: item["pulsar"]):
        hits = p["primary_crossings"]
        reasons = Counter(reason for hit in hits for reason in hit["unpassed"])
        best = min(p["peaks"], key=lambda row: row["p_global"], default={})
        disposition.append(dict(pulsar=p["pulsar"], trials=p["trials"],
                                primary=len(hits), neighborhoods=len(p["peak_neighborhoods"]),
                                surviving_templates=len(p["survivors"]), failed_checks=dict(reasons),
                                best_shape=best.get("template"), best_epoch_mjd=best.get("epoch_mjd"),
                                best_amplitude_us=best.get("amplitude_us"), best_global_p=best.get("p_global"),
                                noise_diagnostics_pass=p["noise_diagnostics"]["passes"]))
    table(dest / "individual-dispositions.csv", disposition,
          ("pulsar", "trials", "primary", "neighborhoods", "surviving_templates", "failed_checks",
           "best_shape", "best_epoch_mjd", "best_amplitude_us", "best_global_p", "noise_diagnostics_pass"))
    table(dest / "individual-sensitivity.csv", limits,
          ("pulsar", "template", "width_d", "epoch_min_mjd", "epoch_max_mjd", "injections",
           "status", "reported_upper_limit_us", "exclusion_withheld", "efficiency_nearest_matching_template"))
    review = audit(dest, ok, array)
    coherence = read(dest / "coherence/summary.json")
    sign = coherence["observed_full_max_tail"]
    selections = Counter(row["epta_selection"] for row in coverage if row["epta"] == "ok")
    if review["individual"] == review["all_support"] == 0:
        summary = (f"No transient passed every check in {len(ok)} NANOGrav pulsars. "
                   f"The blind search found {review['array_internal']} internally significant common-pattern fits. "
                   f"{review['full_sign']} passed the survey-wide pulsar-sign test.")
    else:
        summary = (f"In {len(ok)} NANOGrav pulsars, {review['individual']} individual fits and "
                   f"{review['all_support']} of {review['array_internal']} common-pattern fits passed every check.")
    bullet = summary + f" Signal injections give conditional amplitude limits for {len(reported_limits)} pulsar, waveform and width combinations."
    sums = {path.relative_to(dest).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(dest.rglob("*")) if path.is_file()}
    (dest / "SHA256SUMS").write_text("".join(f"{sha}  {name}\n" for name, sha in sums.items()))
    body = [
        "# Blind pulsar timing transient search", "", f"**{summary}**", "",
        f"The strongest array alignment had a survey-wide sign-randomization p-value of "
        f"**{sign['p_sign_max_bank']:.4g}**. In {sign['exceedances']:,} of {coherence['draws']:,} trials, "
        f"randomly reversing entire pulsar time series produced a peak at least as strong somewhere among "
        f"the {array['n_trials']:,} array templates. This test preserves transient structure within each pulsar "
        f"and assumes sign symmetry across pulsars. After restricting the bank to the "
        f"{len(coherence['diagnostic_passing_pulsars'])} pulsars passing noise diagnostics, the strongest "
        f"alignment had p = {coherence['observed_subset_max_tail']['p_sign_max_bank']:.3g}. "
        "[Randomization results](evidence/coherence/summary.json)", "",
        f"The individual search produced {primary_count} primary template crossings. "
        "Nearby epochs, widths and physical templates can describe the same excursion. "
        "The number of fits is not the number of distinct astrophysical events.", "",
        f"Wideband reductions succeeded for {n_wb} pulsars. Independent EPTA data were usable for {n_epta} of "
        f"{len(cov['epta'])} overlapping pulsars ({selections['DR2full']} DR2full solutions and "
        f"{selections['DR2new']} modern-backend DR2new solutions). For a failed full-release fit the same "
        "phase and residual-quality checks were applied to the published modern-backend subset. "
        "The coverage files retain both attempts. Catalog membership alone does not guarantee observations at a given epoch. "
        "[Coverage](evidence/coverage.csv) · [Per-pulsar dispositions](evidence/individual-dispositions.csv) · "
        "[Individual records](evidence/single) · [Surviving fits](evidence/surviving-fits.json)", "",
        "## Timing signatures and interpretation", "",
        f"Of the {review['array_internal']} internally surviving clock, ephemeris or memory fits, "
        f"{review['full_sign']} passed the whole-bank pulsar-sign control, "
        f"{review['noise_subset_and_sign']} also passed the noise-screened repeat, and "
        f"{review['external_holm']} had multiplicity-adjusted support in EPTA "
        f"({review['external_testable']} had enough external sky coverage). "
        f"The independent-per-visit dispersion fit supported the original direction in {review['dm_support']} "
        f"of {review['dm_testable']} testable fits. On the same visits, {review['matched_raw_support']} "
        f"had directional support before fitting an independent DM per visit. "
        f"Only {review['dm_expected_snr_ge3']} dispersion-controlled fits had enough precision for an expected "
        "signal-to-noise ratio of at least three at the discovery amplitude. "
        "[Every surviving fit, with MJD, amplitude and checks](evidence/EVENTS.md) · "
        "[Physical review by waveform and epoch](COHERENT_REVIEW.md) · [Post-discovery controls](ADJUDICATION.md)", "",
        "The surviving-fits table includes the fitted epoch, amplitude, significance and checks for every surviving template. "
        "Individual survivors require independent EPTA support. Array survivors pass the internal leave-one-pulsar-out check. "
        "Each array record separately reports wideband and EPTA support, inconclusive checks and missing epochs. "
        "A wideband recurrence tests the reduction of the same telescope observations. "
        "An EPTA comparison tests independent telescope data at the discovery epoch. "
        "[Common-mode results and alternatives](evidence/array) · [Published event epochs](PRIOR_EVENTS.md)", "",
        "As a post-search comparison, fixed 60-day timing templates were evaluated at the two published "
        "J1713+0747 dips and the J1643−1224 profile event. The table distinguishes a fixed-epoch response "
        "from a blind trial that crossed the survey threshold. [Known-event controls](evidence/known-event-controls.csv)", "",
        "## Sensitivity and scope", "",
        "![Survey-wide sign test and injected memory recovery](evidence/search-calibration.png)", "",
        f"The individual search reports {len(reported_limits)} pulsar, morphology and width combinations with a conditional amplitude limit. "
        "Blank limits mark combinations without a calibrated exclusion. The individual records give the injected amplitudes, units, "
        "recovery rule and assumptions. The coherent-search tables give injection counts and recovery fractions by signal family. "
        "Use the epoch range and sampling distribution accompanying each limit. "
        "[Individual sensitivity](evidence/individual-sensitivity.csv) · [Array sensitivity](evidence/array/array_limits.csv)", "",
        "## Inputs and reproduction", "",
        "| Release | DOI | Downloaded archive SHA-256 |", "| --- | --- | --- |",
    ]
    external_rows = [r for r in read(dest / "event-review.json") if r["external_support_after_holm"]]
    if external_rows:
        lines = ["### Patterns with independent directional support", "",
                 "EPTA tests the discovery direction at a frozen epoch. The sign test compares that discovery "
                 "with the strongest coincidence anywhere in the NANOGrav bank. Each row links to the full "
                 "record of participants and other checks. [Event ledger](evidence/EVENTS.md)", "",
                 "| Pattern | MJD | NANOGrav amplitude ± SE | EPTA amplitude ± SE | Units | EPTA Holm p | NANOGrav sign-test p |",
                 "| --- | ---: | ---: | ---: | --- | ---: | ---: |"]
        for row in external_rows:
            r = row["coherence"]
            lines.append(f"| {row['family']} {row['shape']} | {row['epoch_mjd']:.0f} | "
                         f"{float(r['discovery_direction_amplitude']):.3g} ± {float(r['discovery_direction_se']):.3g} | "
                         f"{float(r['epta_direction_amplitude']):.3g} ± {float(r['epta_direction_se']):.3g} | "
                         f"{row['units']} | {float(r['epta_p_holm']):.3g} | {float(r['p_sign_max_bank']):.3g} |")
        lines.append("")
        position = body.index("## Sensitivity and scope")
        body[position:position] = lines
    thresholds = ["### Individual amplitude limits", "",
                  "| Timing signature | Width in days | Pulsars with limits | Median limit in µs | Range in µs |",
                  "| --- | ---: | ---: | ---: | ---: |"]
    for name, width in (("spin_step_or_memory", None), ("shapiro_log", 60.), ("dm_tail", 60.)):
        values = [r["reported_upper_limit_us"] for r in reported_limits
                  if r["template"] == name and r["width_d"] == width]
        med = f"{median(values):.4g}" if values else "Unresolved"
        span = f"{min(values):.4g}–{max(values):.4g}" if values else "Unresolved"
        thresholds.append(f"| {name} | {width or '—'} | {len(values)} | {med} | {span} |")
    thresholds.extend(["", "For a spin step the amplitude is the timing drift accumulated in one Julian year. "
                       "The Shapiro coefficient multiplies the logarithmic profile. The DM tail amplitude is the initial "
                       "delay at 1,400 MHz. The linked data give every tested width and morphology. [Method](SINGLE_SEARCH.md)", "",
                       "### Coherent amplitude limits", "",
                  "| Spatial pattern | Time dependence | First amplitude reaching the recovery criterion | "
                  "Largest upper bound over eligible grid cells | Units |",
                  "| --- | --- | ---: | ---: | --- |"])
    def number_or_gap(value):
        return "Unresolved" if value is None else f"{value:.4g}"
    for row in array["sensitivity"]:
        thresholds.append(f"| {row['family']} | {row['shape']} | "
                          f"{number_or_gap(row.get('first_tested_amplitude_90pct'))} | "
                          f"{number_or_gap(row.get('worst_eligible_grid_upper95'))} | {row['units']} |")
    thresholds.extend(["", "The recovery column uses off-grid injections into the measured residuals. "
                       "The upper-bound column gives simultaneous 95% coverage under the noise model for an exact grid waveform. "
                       "[Method](ARRAY_SEARCH.md)", ""])
    insert = body.index("## Inputs and reproduction")
    body[insert:insert] = thresholds
    for release, item in manifest["archives"].items():
        name = {"ng": "NANOGrav 15-year", "epta": "EPTA DR2"}[release]
        body.append(f"| {name} | {item['doi']} | `{item['sha256']}` |")
    body += ["", "Individual timing-file hashes, clock corrections and ephemerides are recorded with the archive receipts. "
             "[Input receipts](evidence/inputs.json) · [Output hashes](evidence/SHA256SUMS) · [Run revision](evidence/run.json)", "",
             "Run the manual workflow at the recorded analysis revision with the same software and inputs. "
             "The workflow contains the exact fetch, reduction, search, injection and verification commands. "
             "[Workflow](../../.github/workflows/nanograv15-transients.yml) · [Reproduction](RUN.md)", ""]
    (folder / "README.md").write_text("\n".join(body))
    catalog(root, label, bullet)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    render(args.root.resolve(), args.evidence.resolve())
