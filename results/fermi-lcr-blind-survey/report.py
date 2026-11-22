"""Write the public findings, limits and quantitative sensitivity from frozen outputs."""
import csv
import hashlib
import json
import os
import sys
from collections import Counter
from datetime import date, timedelta
from pathlib import Path

import periods


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def yes(value):
    return value in ("True", "true", True)


def day(mjd):
    return str(date(1858, 11, 17) + timedelta(days=int(float(mjd))))


def describe(output):
    def obj(path):
        return json.loads((output / path).read_text())

    coverage, stats = obj("coverage.json"), obj("period-summary.json")
    period_table = read_csv(output / "periods.csv")
    state_table = read_csv(output / "states.csv")
    transitions = read_csv(output / "transitions.csv")
    review = obj("review/index.json")
    meta = obj("index/catalog_origin.json")
    receipts = [json.loads(line) for line in (output / "inputs.jsonl").read_text().splitlines()]
    valid = [r for r in receipts if r["status"] == "ok"]
    eligible = [r for r in period_table if r["status"] == "tested"]
    eligible_states = sum(int(r["n_joint_splits"]) > 0 for r in state_table)
    recovered = stats["aggregate_recovery"]
    global_hits = [r for r in eligible if yes(r["global_ok"])]
    passed = [r for r in global_hits if yes(r["passes"])]
    if not valid or not eligible:
        raise RuntimeError("No validated periodic search. A scientific null cannot be written.")
    roster = stats["roster"]
    ntested = len(eligible)
    cut = stats["global_threshold"]
    state_count = coverage["state_sources"]
    if passed:
        result = f"{len(passed)} source(s) passed the catalog-wide periodicity and stability tests; literature and counterpart dispositions are recorded below."
    else:
        result = f"No periodic signal passed all checks among {ntested:,} Fermi gamma-ray sources."
    state_claim = (f"{len(transitions)} sustained factor-ten transition(s) met the two-cadence flux criteria in {state_count} source(s)"
                   if state_count else "No sustained factor-ten transition met the two-cadence flux criteria")
    if passed or state_count:
        label = "Measured"
        title = "A search for repeating gamma rays and tenfold changes in brightness."
        if state_count:
            title = "Some gamma-ray sources stayed at least ten times brighter or fainter for a year."
        detail = f"{len(passed)} of {ntested:,} sources passed the timing tests. {state_claim} across {eligible_states:,} usable sources."
    else:
        label = "Null result"
        title = "No persistent gamma-ray clock passed a survey-wide test."
        detail = f"In {ntested:,} usable sources no rhythm passed the full set of checks. The search also found no sustained tenfold switch among {eligible_states:,} sources with adequate coverage."
    (output / "catalog-entry.json").write_text(json.dumps(dict(label=label, title=title, detail=detail), indent=2) + "\n")
    source = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat/LightCurveRepository/"
    run_url = ("https://github.com/" + os.environ["GITHUB_REPOSITORY"] + "/actions/runs/" + os.environ["GITHUB_RUN_ID"]
               if os.environ.get("GITHUB_RUN_ID") else None)
    text = [
        "# Repeating gamma rays and sustained changes in Fermi-LAT sources", "",
        f"**{result}** {len(global_hits)} cleared the simulated catalog-wide threshold before stability and quality checks. {state_claim} among {eligible_states:,} sources with adequate state coverage. The census measures gamma-ray flux changes. A binary interpretation requires independent evidence.", "",
        "## Coverage and measured sensitivity", "",
        f"The roster contained {roster:,} sources selected by {meta['selection']} from {meta['n_dr2']:,} DR2 entries. The archive returned usable data for {len(valid):,} of {len(receipts):,} requested source–cadence pairs. Observed bin centers span {day(min(r['first_mjd'] for r in valid))} through {day(max(r['last_mjd'] for r in valid))} across available products. Each row records its own endpoints. [Archive](%s) · [Input receipts](outputs/inputs.jsonl) · [Frozen roster](outputs/index/catalog.json)" % source, "",
        "| Product | Valid curves | Detections | Upper-limit bins |", "|---|---:|---:|---:|"]
    for cadence in ("3day", "weekly", "monthly"):
        rows = [r for r in valid if r["cadence"] == cadence]
        text.append(f"| {cadence} | {len(rows):,} | {sum(r['n_flux'] for r in rows):,} | {sum(r['n_ul'] for r in rows):,} |")
    text += ["", f"The shortest and longest grid periods tested anywhere in the catalog were {min(float(r['pmin_days']) for r in eligible):.1f} and {max(float(r['pmax_days']) for r in eligible):.1f} days. Each source has its own range. The trial correction used the maximum over all {sum(int(r['n_frequencies']) for r in eligible):,} source–frequency combinations in each of {stats['replicates']:,} simulated catalogs. The critical normalized power was {cut:.4g}. [Source scores, p-values and exclusions](outputs/periods.csv) · [Catalog maxima](outputs/period-summary.json)", "",
             "The table gives the fraction of injected signals recovered at the survey threshold with the correct frequency and stable halves. Amplitude is the fractional sinusoidal semi-amplitude multiplying the simulated flux. Counterpart, short-cadence and censored-flux tests can reduce final acceptance further. [Per-source injection counts](outputs/injections.csv)", "",
             "| Period (days) | Semi-amplitude | Recovered / injected | Recovery |", "|---:|---:|---:|---:|"]
    for r in recovered:
        text.append(f"| {r['period_days']:.0f} | {r['semi_amplitude']:.0%} | {r['recovered']:,} / {r['n']:,} | {r['fraction']:.1%} |")
    text += ["", "## Signals, state changes and prior reports", ""]
    pg = next((r for r in period_table if r["source"] == "4FGL J1555.7+1111"), None)
    if pg and pg["status"] == "tested":
        text += [f"For the published PG 1553+113 comparison, the strongest blind peak was {float(pg['period_days']):.1f} days, with catalog-wide p = {float(pg['p_global']):.4g}. Its independent-half check {'passed' if yes(pg['halves_ok']) else 'failed'}. This source received the same search grid and threshold as every other source.", ""]
    if global_hits:
        text += ["| Source | Best period (days) | Catalog p | Stable halves | Spectral fit | Short cadences | Censor bounds | Annual veto clear |",
                 "|---|---:|---:|---|---|---|---|---|"]
        for r in sorted(global_hits,key=lambda r:float(r["p_global"])):
            cols = ["pass" if yes(r[k]) else "fail" for k in ("halves_ok", "spectrum_ok", "cadence_ok", "censor_ok")]
            text.append(f"| {r['source']} | {float(r['period_days']):.1f} | {float(r['p_global']):.4g} | " + " | ".join(cols) + (" | fail |" if yes(r["annual_alias"]) else " | pass |"))
        text += [""]
    if not passed and not state_count:
        text += ["No source qualified for literature or counterpart follow-up under the prespecified rules.", ""]
    if transitions:
        text += ["| Source | Direction | Split MJD | Minimum flux ratio |", "|---|---|---:|---:|"]
        for r in sorted(transitions,key=lambda r:-float(r["ratio_floor"])):
            text.append(f"| {r['source']} | {r['direction']} | {float(r['split_mjd']):.1f} | {float(r['ratio_floor']):.2f} |")
        text += [""]
    if review:
        text += ["Each dossier records alias searches, comparison papers, catalog and counterpart matches, service failures and ambiguous identifications. A source with an incomplete review has no novelty designation.", ""]
        for rev in review:
            text.append(f"- {rev['name']} — {rev['status'].replace('_',' ')}. [Evidence](outputs/review/{rev['file']})")
        text += [""]
    text += ["The state census tested adjacent year-long states across a monthly grid of split dates. Each half-year on the bright side had to exceed both half-years on the faint side by tenfold in monthly and weekly products. It used per-bin uncertainty and upper-limit bounds, with missing bins padded conservatively. At least 5 of 6 monthly and 21 of 26 weekly bins per half-year were required. A high or low state can reverse outside this two-year interval. These descriptive bounds do not assign a stochastic false-alarm probability to a flux change. [Per-source census](outputs/states.csv) · [Passing intervals](outputs/transitions.csv)", "",
        "## Timing model and limits", "", periods.METHOD.strip(), "",
        "## Inputs and reproduction", "",
        f"The SHA-256 of the original DR2 source-table response is `{meta['sha256']}`. The per-response receipt gives the URL, SHA-256 of the uncompressed response bytes, bin counts and any archive error. The normalized time, flux, error, upper-limit, test-statistic and fit arrays used for the search are deposited in 24 compressed shards. [File checksums](outputs/SHA256SUMS) · [Selection evidence](outputs/index/selection-evidence.json) · [Coverage](outputs/coverage.json)", "",
        "Run from the repository root with Python 3.12 and the dependencies below. The replay uses only the deposited curves and repeats the same deterministic simulations. It writes into a new directory.", "",
        "```bash", "study=results/fermi-lcr-blind-survey", 'python -m pip install -r "$study/requirements.txt"', 'OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python "$study/replay.py" "$study/outputs" replay-fermi', "(cd \"$study/outputs\" && sha256sum -c SHA256SUMS)", "python -m open_space.verify", "python -m unittest discover -s tests", "```", "",
        "To repeat a live retrieval use the workflow or the following commands. The archive is mutable, so a later retrieval can produce different inputs. The normalized snapshot fixes the exact search inputs. " + (f"[Original workflow run]({run_url})" if run_url else ""), "",
        "```bash", 'python "$study/acquire.py" index --out fresh-index', 'for shard in $(seq 0 23); do', '  python "$study/acquire.py" fetch --index fresh-index/catalog.json --out "data-$shard" --shard "$shard" --shards 24 --workers 1', '  OPENBLAS_NUM_THREADS=1 python "$study/survey.py" shard --data "data-$shard" --index fresh-index/catalog.json --out "fresh/part-$shard" --shard "$shard" --shards 24', '  cp "data-$shard/manifest.jsonl" "fresh/part-$shard/"', 'done', 'python "$study/survey.py" reduce --index fresh-index/catalog.json --work fresh --out fresh-results', 'python "$study/review.py" fresh-results', "```", ""]
    (output / "README.md").write_text("\n".join(text))
    files = sorted(p for p in output.rglob("*") if p.is_file() and p.name not in ("README.md", "SHA256SUMS"))
    (output / "SHA256SUMS").write_text("".join(f"{digest(p)}  {p.relative_to(output)}\n" for p in files))


if __name__ == "__main__":
    describe(Path(sys.argv[1]))
