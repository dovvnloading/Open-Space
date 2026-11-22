"""Reduce frozen observation ledgers without contacting any archive."""
import argparse
import csv
import math
import os
import re
from collections import Counter, defaultdict

from scipy.stats import beta

from archive import dump
from campaign import EVIDENCE, HERE, N_CONTROL, read, read_obs
from search import DURATIONS, _time_support
from vet import separation

ROOT = HERE.parents[1]
TAG = "results/chandra_acis_fxt/README.md"
FINAL_VETO = {"foreground_veto", "known_transient"}


def bounds(k, n):
    if not n:
        return None, None
    return (float(beta.ppf(.025, k, n-k+1)) if k else 0.,
            float(beta.ppf(.975, k+1, n-k)) if k < n else 1.)


def rows_to_csv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def distinct(rows):
    """Collapse different data streams of the same ObsID for coincident bursts."""
    retained = defaultdict(list)
    for row in sorted(rows, key=lambda r: -r["seed_score"]):
        for other in retained[row["obsid"]]:
            if min(row["stop_met_s"], other["stop_met_s"]) <= max(row["start_met_s"], other["start_met_s"]):
                continue
            if separation(row["ra_deg"], row["dec_deg"], other["ra_deg"], other["dec_deg"]) <= math.hypot(
                row["position_error_arcsec"], other["position_error_arcsec"]
            ) * math.sqrt(2):
                break
        else:
            retained[row["obsid"]].append(row)
    return [row for group in retained.values() for row in group]


def kept(rows):
    return [r for r in rows if r["instrument_pass"] and r["vet"]["status"] not in FINAL_VETO]


def collect():
    population = read(EVIDENCE / "inventory/inventory.json")
    exceptions = read(EVIDENCE / "inventory/inventory_exceptions.json")
    observations = []
    trials = defaultdict(list)
    tally = defaultdict(Counter)
    exposure = Counter()
    receipt_rows, injections = [], []
    complete = searched = 0
    for entry in population:
        oid = entry["obsid"]
        path = EVIDENCE / "observations" / f"{oid}.json.gz"
        row = read_obs(path) if path.exists() else {"status": "unattempted", "products": [], "errors": []}
        acquired = {(r["filename"], r["sha256"]) for r in row.get("acquired", [])}
        scanned = {(r["input"]["filename"], r["input"]["sha256"]) for r in row["products"] if r.get("complete")}
        finished = bool(acquired) and acquired == scanned and row["status"] in ("complete", "followup_incomplete")
        observations.append({"obsid": oid, "status": row["status"], "scanned_all_products": finished,
                             "n_products": len(row["products"]), "errors": " | ".join(row["errors"])})
        complete += finished and row["status"] == "complete"
        for src in row.get("acquired", []):
            receipt_rows.append({**src, "included_in_inference": finished})
        if not finished:
            continue
        searched += 1
        intervals = defaultdict(list)
        local = defaultdict(list)
        for product in row["products"]:
            result = product["search"]
            for length, windows in result["start_support_by_duration"].items():
                intervals[length].extend(windows)
            streams = [(-1, result["triggers"])] + [(c["replica"], c["triggers"]) for c in result["controls"]]
            if sorted(i for i, _ in streams) != list(range(-1, N_CONTROL)):
                raise ValueError(f"{oid} control ensemble incomplete")
            for replica, events in streams:
                for event in events:
                    event.update(obsid=oid, input_name=product["input"]["filename"], input_url=product["input"]["url"],
                                 input_sha256=product["input"]["sha256"], cache=event.get("catalog_cache", row["catalog_cache"]),
                                 proposal=entry.get("proposal_type"), target=entry.get("target_name"))
                    # Full diagnostics and catalog bodies remain in the ObsID ledger.
                    compact = {k: v for k, v in event.items() if k not in ("vet", "flux")}
                    compact["vet"] = {k: event["vet"].get(k) for k in ("status", "host")}
                    compact["flux"] = {k: event.get("flux", {}).get(k) for k in
                                       ("peak_flux_erg_cm2_s", "flux_counting_sigma_erg_cm2_s")}
                    local[replica].append(compact)
            for j in result["injections"]:
                key = (j["duration_s"], j["counts"], j["sigma_sky_px"])
                tally[key]["planned"] += 1
                if j["attempted"]:
                    tally[key]["attempted"] += 1
                    tally[key]["recovered"] += bool(j["recovered"])
                injections.append(dict(j, obsid=oid, input_sha256=product["input"]["sha256"]))
        for length, spans in intervals.items():
            exposure[length] += sum((b-a)*f for a, b, f in _time_support(spans))
        for replica, events in local.items():
            # Keep every raw trigger. Deduplicate only within the selection being counted.
            trials[replica].extend(events)
    return population, exceptions, observations, trials, tally, exposure, receipt_rows, injections, complete, searched


def run():
    (population, exceptions, observations, trials, tally, exposure, receipts,
     injections, complete, searched) = collect()
    out = EVIDENCE
    eligible = {i: distinct(kept(trials[i])) for i in range(-1, N_CONTROL)}
    all_pass = {i: distinct([t for t in trials[i] if t["instrument_pass"]]) for i in range(-1, N_CONTROL)}
    controls = []
    for i in range(N_CONTROL):
        controls.append({"replica": i, "scanned_obsids": searched, "triggers": len(trials[i]),
                         "instrument_pass": len(all_pass[i]), "surviving_or_unresolved": len(eligible[i]),
                         "incomplete_queries": sum(t["vet"]["status"] == "catalog_incomplete" for t in trials[i]),
                         "instrument_maximum": max((t["seed_score"] for t in all_pass[i]), default=None),
                         "maximum_score": max((t["seed_score"] for t in eligible[i]), default=None)})
    rows_to_csv(out / "control_trials.csv", list(controls[0]), controls)
    real = eligible[-1]
    leads = []
    for t in real:
        reference = "maximum_score" if complete == searched else "instrument_maximum"
        n_extreme = sum(c[reference] is not None and c[reference] >= t["seed_score"] for c in controls)
        host = t["vet"].get("host") or {}
        near = [other for other in real if other["obsid"] != t["obsid"] and
                separation(t["ra_deg"], t["dec_deg"], other["ra_deg"], other["dec_deg"]) <= math.hypot(
                    t["position_error_arcsec"], other["position_error_arcsec"]) * math.sqrt(2)]
        fap = (1 + n_extreme) / (1 + N_CONTROL)
        fap_lo, fap_hi = bounds(n_extreme, N_CONTROL)
        leads.append({**{k:t[k] for k in ("obsid", "ra_deg", "dec_deg", "position_error_arcsec", "event_mjd",
                      "duration_s", "trial_duration_s", "peak_rate_ct_s", "peak_live_s", "trigger_key",
                      "input_name","input_url","input_sha256", "proposal","target")},
                      "host": host.get("name"), "host_z": host.get("redshift"), "status": t["vet"]["status"],
                      "host_offset_arcsec": host.get("separation_arcsec"), "host_chance_upper95":host.get("chance_upper95"),
                      "peak_flux_erg_cm2_s": t.get("flux", {}).get("peak_flux_erg_cm2_s"),
                      "counting_sigma": t.get("flux", {}).get("flux_counting_sigma_erg_cm2_s"),
                      "other_trigger_obsids": sorted({r["obsid"] for r in near}),
                      "null_reference": reference, "survey_fap": fap, "control_exceedances": n_extreme, "fap_lo95": fap_lo, "fap_hi95": fap_hi})
    dump(out / "leads.json", leads)
    rows_to_csv(out / "obsid_coverage.csv", ["obsid","status","scanned_all_products","n_products","errors"], observations)
    rows_to_csv(out / "event_receipts.csv", ["obsid","filename","url","resolved_url","bytes","sha256","fits_sha256",
                "revision","cycle","included_in_inference"], receipts)
    grid = []
    for (duration, counts, sigma), n in sorted(tally.items()):
        lo, hi = bounds(n["recovered"], n["attempted"])
        grid.append({"duration_s":duration, "input_counts":counts,"psf_sigma_sky_px":sigma,
                     "scheduled":n["planned"],"attempted":n["attempted"], "recovered":n["recovered"],
                     "fraction": n["recovered"]/n["attempted"] if n["attempted"] else None, "lo95":lo, "hi95":hi})
    rows_to_csv(out / "completeness.csv", ["duration_s","input_counts","psf_sigma_sky_px","scheduled","attempted",
                "recovered","fraction","lo95","hi95"], grid)
    attempted = sum(i["attempted"] for i in injections)
    recovered = sum(bool(i["recovered"]) for i in injections)
    nreal = len(trials[-1])
    weighted = exposure["60"]
    has_measurement = weighted > 0 and attempted > 0
    # A discovery disposition requires documented external review. Unmatched flares
    # remain searchable leads even if all the positional services returned no rows.
    has_null = has_measurement and recovered > 0 and not real and complete == searched
    label = "Null result" if has_null else "Measured" if has_measurement else None
    exp_text = f"{weighted/1e6:.6f} Ms of eligible 60-second flare-start live time"
    if has_null:
        headline = f"No FXT survived in {searched:,} scanned ObsIDs covering {exp_text}."
    elif has_measurement:
        headline = f"The scan covered {exp_text} in {searched:,} ObsIDs and retained {len(real):,} unresolved leads."
    else:
        headline = "The run did not establish a positive exposure with measured injection completeness."
    avg = sum(len(eligible[i]) for i in range(N_CONTROL))/N_CONTROL
    summary = {"label":label, "inventory_obsids":len(population), "inventory_exceptions":len(exceptions),
               "completed_obsids":complete, "scanned_obsids":searched, "release_window_complete":
               bool(population) and complete == len(population) and not any(e["reason"] != "check_mode_from_FITS" for e in exceptions),
               "input_products":len(receipts), "exposure_by_duration_s":dict(exposure),
               "real_raw_triggers":nreal,"real_instrument_pass":len(all_pass[-1]),"scientific_leads":len(real),
               "control_replicas":N_CONTROL, "mean_control_survivors":avg,
               "control_events_per_Ms":avg*1e6/weighted if weighted else None,
               "conditional_null_upper95_events_per_Ms": -math.log(.05)*1e6/weighted if has_null else None,
               "detected_event_upper95_per_Ms_by_duration":
                   {str(d): -math.log(.05)*1e6/exposure[str(d)] if exposure[str(d)] else None for d in DURATIONS}
                   if has_null else None,
               "scheduled_injections":len(injections),"attempted_injections":attempted,"recovered_injections":recovered,
               "run_url":f"https://github.com/{os.environ['GITHUB_REPOSITORY']}/actions/runs/{os.environ['GITHUB_RUN_ID']}"}
    dump(out / "summary.json", summary)
    lines = ["# Chandra ACIS fast X-ray transient search", "", f"**{headline}**", "",
             f"The frozen release inventory contains {len(population):,} ObsIDs. All data products were scanned in "
             f"{searched:,}. Catalog calls and needed photometry completed in {complete:,}. "
             f"The release-date or mode ledger contains {len(exceptions):,} exceptions.",
             "", f"Across {N_CONTROL} full-survey scrambles the mean number passing the same instrumental and "
             f"counterpart screen was {avg:.3f}. The real fields produced {nreal:,} raw triggers. "
             f"Of {attempted:,} executed injections, {recovered:,} passed instrumental selection. "
             f"{len(injections)-attempted:,} additional draws lacked usable observing time or geometry.", "",
             "| Trial duration (s) | Eligible field-start live seconds |", "|---:|---:|"]
    lines.extend(f"| {s} | {exposure[str(s)]:,.1f} |" for s in DURATIONS)
    if has_null:
        lines.extend(["", f"The one-sided 95% Poisson upper bound is {summary['conditional_null_upper95_events_per_Ms']:.3g} "
                      "detected events per Ms for this field mixture, using 60-second start support. "
                      "Intrinsic sky rates require effective survey area and a source population model."])
    lines.extend(["","[Input checksums](evidence/event_receipts.csv) · [Coverage](evidence/obsid_coverage.csv) · "
                  "[Scrambles](evidence/control_trials.csv) · [Injection grid](evidence/completeness.csv) · "
                  "[All unresolved leads with positions, durations, fluxes, hosts and input hashes](evidence/leads.json)",
                  "", "Each compressed ObsID ledger records all raw triggers, vetoes, controls, catalog receipts "
                  "and injection outcomes. Leads require examination of the target and proposal, independent visits, "
                  "the detector image, and contemporary notices before a claim of novelty.",
                  "", "[Methods and reproduction](METHODS.md) · [Actions run](" + summary["run_url"] + ")"])
    (HERE / "README.md").write_text("\n".join(lines)+"\n")
    update_catalog(headline, label)


def update_catalog(headline, label):
    path = ROOT / "README.md"
    lines = path.read_text().splitlines(keepends=True)
    indices = [i for i,s in enumerate(lines) if s.startswith("- ") and TAG in s]
    if len(indices) != 1:
        raise ValueError("Expected exactly one Chandra survey catalog entry")
    idx = indices[0]
    suffix = f" *{label}* ·" if label else ""
    lines[idx] = f"- **Chandra ACIS search for minute-to-hour X-ray flares.** {headline}{suffix} [Details]({TAG})\n"
    tally = Counter(m.group(1) for s in lines if s.startswith("- ") for m in
                    re.finditer(r"(?<!\*)\*(Candidate|Measured|Null result|Correction|Prediction)\*(?!\*)", s))
    # Rewrite only the existing counts in the catalog summary.
    edits = [(r"\d+(?= entries)",sum(tally.values())),(r"\d+(?= new measurements)",tally["Measured"]),
             (r"\d+(?= candidate discoveries)",tally["Candidate"]),(r"\d+(?= forecasts)",tally["Prediction"]),
             (r"\d+(?= corrections)",tally["Correction"]),(r"\d+(?= carefully recorded null)",tally["Null result"])]
    for pattern, value in edits:
        lines[8], n = re.subn(pattern, str(value), lines[8])
        if n != 1:
            raise ValueError("Catalog line 9 changed; inspect the summary before committing")
    path.write_text("".join(lines))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--check-coverage", action="store_true")
    if p.parse_args().check_coverage:
        if not read(EVIDENCE / "summary.json")["release_window_complete"]:
            raise SystemExit("Partial coverage or source exceptions; see deposited ledger before retrying")
    else:
        run()
