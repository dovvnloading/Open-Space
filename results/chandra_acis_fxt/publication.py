"""Publish measured survey outcomes from exact, manually adjudicated lead identities."""
import argparse
from collections import Counter, defaultdict
import csv
from functools import lru_cache
import hashlib
import json
import math
import os
from pathlib import Path
import re
from urllib.parse import quote, urlsplit
import unicodedata

from scipy.stats import chi2

from archive import check_archive_url, dump
from campaign import HERE, identity, read, read_obs
import report
from review import fmt
from search import DURATIONS, SEED_P

PRIMARY = "dated"
GROUPS = ("date_bracketed_sky", "unresolved_date_sky", "non_sky_or_calibration_review")
NAMES = {PRIMARY:"Dated census", GROUPS[0]:"Release bracket established", GROUPS[1]:"Release date unresolved",
         GROUPS[2]:"Detector or calibration review"}
HEX = re.compile(r"[a-f0-9]{64}")
ALLOWED = {"candidate", "unresolved", "rejected"}


def need(condition, why):
    if not condition:
        raise ValueError(why)


@lru_cache(maxsize=8192)
def file_sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def csv_rows(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def key(row):
    oid, digest, trigger = int(row["obsid"]), row["input_sha256"], row["trigger_key"]
    need(oid > 0 and HEX.fullmatch(digest) and re.fullmatch(r"real/(0|[1-9]\d*)", trigger),
         "Invalid retained-trigger identity")
    return oid, digest, trigger


def local(value):
    path = (HERE / value).resolve()
    need(path.is_relative_to(HERE) and path.is_file(), f"Missing local evidence {value}")
    return path


def href(path):
    return quote(path.relative_to(HERE).as_posix(), safe="/-._")


def link(path, label):
    return f"[{label}]({href(path)})"


def slug(text):
    return "".join(c for c in text.strip().lower() if c in "-_ " or
                   not unicodedata.category(c).startswith(("P", "S"))).replace(" ", "-")


def references(items, ident, row, root):
    """At least one local reference must address this exact trigger."""
    need(isinstance(items, list) and items, "A manual decision requires source evidence")
    oid, h, trigger = ident
    direct = root / "photometry" / str(oid) / row["input_name"] / trigger
    specific = False
    links, hashes = [], {}
    for n, value in enumerate(items, 1):
        u = urlsplit(value)
        if u.scheme:
            need(u.scheme == "https" and u.hostname and not u.username and not u.password
                 and not any(c in value for c in " ()<>\r\n"), "Invalid external evidence reference")
            links.append(f"[source {n}]({value})")
            continue
        need(not u.netloc and not u.query and not Path(u.path).is_absolute(), "Expected a survey-relative reference")
        path = local(u.path)
        need(path not in {HERE/"README.md",HERE/"evidence/publication.json"}, "Circular publication evidence")
        hashes[href(path)] = file_sha(path)
        content = path.read_text() if path.suffix == ".md" else ""
        sections = re.split(r"(?m)^(?=#{1,6} )", content)
        if u.fragment:
            sections = [s for s in sections if s and s.startswith("#") and
                        slug(s.splitlines()[0].lstrip("# ")) == u.fragment]
            need(len(sections) == 1, f"Invalid or ambiguous review-note heading in {value}")
        if path.parent == direct and path.name in ("diagnostics.png", "diagnostics.json", "event_rows.csv.gz",
                                                   "lightcurves.csv", "flux.json", "peak_weights.fits.gz"):
            specific = True
        for section in sections:
            lines = section.splitlines()
            chunks = [line for line in lines if line.startswith("|")]
            chunks.append("\n".join(line for line in lines if not line.startswith("|")))
            if any(h in chunk and re.search(rf"\b{oid}\b", chunk) and
                   re.search(rf"(?<!\w){re.escape(trigger)}(?!\d)", chunk) for chunk in chunks):
                specific = True
        links.append(f"[evidence {n}]({href(path)}" + ("#" + quote(u.fragment, safe="._-") if u.fragment else "") + ")")
    need(specific, f"No local review section, table row or diagnostic identifies {ident}")
    return " · ".join(links), hashes


def take(files, path, csv_file=False):
    files[href(path)] = file_sha(path)
    return csv_rows(path) if csv_file else read(path)


def read_group(name, root, reports, meta, files):
    primary = name == PRIMARY
    exposure = meta["exposure_by_duration_s" if primary else "eligible_start_support_s"]
    exposure = {str(d):float(exposure.get(str(d), 0)) for d in DURATIONS}
    need(all(math.isfinite(v) and v >= 0 for v in exposure.values()), "Invalid start support")
    coverage = take(files, reports / ("obsid_coverage.csv" if primary else "coverage.csv"), True)
    receipts = take(files, reports / ("event_receipts.csv" if primary else "event_inputs.csv"), True)
    controls = take(files, reports / ("control_trials.csv" if primary else "control_counts.json"), primary)
    grid = take(files, reports / ("completeness.csv" if primary else "injections.json"), primary)
    leads = take(files, reports / "leads.json")
    ids = {int(r["obsid"]) for r in coverage}
    valid = {int(r["obsid"]) for r in coverage if r["scanned_all_products"].lower() == "true"}
    complete = sum(int(r["obsid"]) in valid and r["status"] == "complete" for r in coverage)
    need(len(ids) == len(coverage) and len(valid) == meta["scanned_obsids" if primary else "fully_scanned"]
         and complete == meta["completed_obsids" if primary else "followup_complete"], "Inconsistent observing coverage")
    need(len(ids) == meta["inventory_obsids"] if primary else ids == set(meta["obsids"]), "Inconsistent cohort membership")
    eligible = {(int(r["obsid"]),r["sha256"]):r for r in receipts if r["included_in_inference"].lower() == "true"}
    need({k[0] for k in eligible} == valid, "Receipt membership differs from scanned ObsIDs")
    by_duration = defaultdict(Counter)
    for r in grid:
        planned, attempted, recovered = (int(r[k]) for k in ("scheduled" if primary else "planned","attempted","recovered"))
        need(0 <= recovered <= attempted <= planned and int(r["duration_s"]) in DURATIONS, "Invalid injection grid")
        by_duration[int(r["duration_s"])].update(planned=planned, attempted=attempted, recovered=recovered)
    total = sum(by_duration.values(),Counter())
    need(total["attempted"] == meta["attempted_injections"] and total["recovered"] == meta["recovered_injections"],
         "Injection totals disagree with summary")
    if primary:
        need(total["planned"] == meta["scheduled_injections"] and len(receipts) == meta["input_products"],
             "Primary acquisition or injection-count mismatch")
    need(sorted(int(r["replica"]) for r in controls) == list(range(24)), "Incomplete or repeated control replica")
    ctrl = [int(c["surviving_or_unresolved" if primary else "retained_or_unresolved"]) for c in controls]
    mean = meta["mean_control_survivors" if primary else "mean_control_retained_or_unresolved"]
    need(all(i >= 0 for i in ctrl) and (not valid and mean is None or math.isclose(sum(ctrl)/24,mean)),
         "Scrambled-control count mismatch")
    if primary:
        need(meta["control_replicas"] == 24 and all(int(c["scanned_obsids"]) == len(valid) for c in controls),
             "Control exposure mismatch")
    need(len(leads) == meta["scientific_leads" if primary else "surviving_or_unresolved"] and
         len({key(r) for r in leads}) == len(leads), "Lead count or identity mismatch")
    return dict(name=name, root=root, reports=reports, exposure=exposure, ids=ids, scanned=len(valid),
                complete=complete, valid=valid, inputs=eligible, grid=by_duration, total=total, control=ctrl,
                control_rows=controls, leads=leads, raw=meta["real_raw_triggers"], instrument=meta.get("real_instrument_pass"))


@lru_cache(maxsize=1)
def obs(path):
    return read_obs(path)


def measurements(g, files):
    primary = g["name"] == PRIMARY
    column = "maximum_score" if g["scanned"] == g["complete"] else "instrument_maximum"
    if primary:
        maxima = [float(r[column]) if r[column] else None for r in g["control_rows"]]
    elif g["leads"]:
        maxima = [None] * 24
        for oid in sorted(g["valid"]):
            p = g["root"] / "observations" / f"{oid}.json.gz"
            files[href(p)] = file_sha(p)
            for prod in obs(p)["products"]:
                for c in prod["search"]["controls"]:
                    scores = [t["seed_score"] for t in c["triggers"] if t["instrument_pass"] and
                              (column == "instrument_maximum" or t["vet"]["status"] not in report.FINAL_VETO)]
                    i = c["replica"]
                    if scores:
                        maxima[i] = max([*scores] + ([] if maxima[i] is None else [maxima[i]]))
    out = {}
    for lead in sorted(g["leads"], key=key):
        ident = key(lead)
        oid, h, trigger = ident
        p = g["root"] / "observations" / f"{oid}.json.gz"
        files[href(p)] = file_sha(p)
        ledger = obs(p)
        receipt = g["inputs"][(oid,h)]
        check_archive_url(receipt["url"])
        prod, = [r for r in ledger["products"] if r["input"]["sha256"] == h]
        event = prod["search"]["triggers"][int(trigger.split("/")[1])]
        need(prod["complete"] and prod["input"]["url"] == lead["input_url"] == receipt["url"] and
             receipt["filename"] == lead["input_name"] == prod["input"]["filename"] and int(receipt["bytes"]) == prod["input"]["bytes"]
             and HEX.fullmatch(receipt["fits_sha256"]) and receipt["fits_sha256"] == prod["input"]["fits_sha256"]
             and event["trigger_key"] == trigger and event["instrument_pass"]
             and event["vet"]["status"] not in report.FINAL_VETO, "Lead and exact acquired trigger disagree")
        for field in ("ra_deg","dec_deg","position_error_arcsec","event_mjd","duration_s","trial_duration_s"):
            need(event[field] == lead[field] and math.isfinite(float(event[field])), "Lead measurement mismatch")
        n = sum(m is not None and m >= event["seed_score"] for m in maxima)
        fap = (1+n)/25
        need(math.isclose(fap,lead["survey_fap" if primary else "supplement_fap"]), "FAP differs from control maxima")
        if primary:
            need(lead["control_exceedances"] == n and lead["null_reference"] == column and
                 lead["status"] == event["vet"]["status"] and
                 lead["peak_flux_erg_cm2_s"] == event.get("flux",{}).get("peak_flux_erg_cm2_s"),
                 "Primary follow-up differs from deposited ledger")
        lo, hi = report.bounds(n,24)
        if primary:
            need(math.isclose(lo,lead["fap_lo95"]) and math.isclose(hi,lead["fap_hi95"]), "FAP precision mismatch")
        out[ident] = dict(lead=lead, event=event, receipt=receipt, group=g, fap=fap, n=n, lo=lo, hi=hi,
                          control_reference=column, search_seed=prod["search"]["seed"],
                          revision=prod["input"].get("revision"), ledger=href(p), ledger_sha256=files[href(p)])
    return out


def load(supplement):
    root = HERE / "evidence"
    files = {}
    meta = take(files, root/"summary.json")
    need(supplement.is_relative_to(root/"supplemental_search"), "Select one archived supplemental search")
    other = take(files, supplement/"summary.json")
    released = take(files, supplement/"release_bounds.json")["obsids"]
    freeze = take(files, root/"inventory/campaign.json")
    need(freeze["implementation"] == identity() and freeze["implementation"]["vet.py"] ==
         "9f3615c635e42d60425be4738dd3ffd62dfe8aee74bc717736df951af2a996e8", "Use the corrected, deposited primary campaign")
    revision = (root/freeze["followup_revision"]).resolve()
    need(revision.is_relative_to(root/"history"), "Invalid corrective revision")
    corrected = take(files,revision/"result.json")
    groups = [read_group(PRIMARY, root, root, meta, files)]
    for ident, item in corrected["changed_or_attempted"].items():
        p = root/"observations"/f"{int(ident)}.json.gz"
        need((file_sha(p) if p.is_file() else None) == item["after_sha256"], "Revised observation ledger has changed")
    need(set(other["groups"]) == set(GROUPS), "Unexpected supplemental cohorts")
    assigned = defaultdict(set)
    for row in released:
        cohort = GROUPS[2] if row["sky_status"] != "nominal_imaging" else GROUPS[0] if row["proven_in_window"] else GROUPS[1]
        assigned[cohort].add(int(row["obsid"]))
    need(sum(map(len,assigned.values())) == len(released), "Duplicate supplemental release bound")
    for name in GROUPS:
        g = read_group(name,supplement,supplement/name,other["groups"][name],files)
        need(g["ids"] == assigned[name], "Supplemental date or calibration membership differs from its receipts")
        groups.append(g)
    seen = set()
    leads = {}
    for g in groups:
        need(not seen & g["ids"], "Overlapping observation cohorts")
        seen |= g["ids"]
        leads.update(measurements(g,files))
    sig = hashlib.sha256(json.dumps(files,sort_keys=True,separators=(",",":")).encode()).hexdigest()
    return groups, leads, sig, files


def flux_components(flux):
    """Reconcile the signed estimator with both measured nonnegative components."""
    net, sigma, gross, background = (float(flux.get(k, math.nan)) for k in
        ("peak_flux_erg_cm2_s", "flux_counting_sigma_erg_cm2_s", "gross_flux_erg_cm2_s",
         "scaled_background_flux_erg_cm2_s"))
    need(all(math.isfinite(x) for x in (net, sigma, gross, background)) and sigma > 0 and
         gross > 0 and background >= 0 and math.isclose(net, gross-background, rel_tol=1e-9, abs_tol=1e-30),
         "Invalid signed peak photometry or inconsistent aperture and background components")
    return gross, background


REVIEW_FIELDS = (
    "cohort obsid trigger_key disposition ra_deg dec_deg position_error_arcsec event_mjd duration_s trial_duration_s "
    "ccd_id start_met_s stop_met_s peak_start_met_s peak_window_s peak_live_s peak_counts peak_annulus_counts "
    "on_counts annulus_counts annulus_area_ratio before_counts before_s after_counts after_s "
    "photometry_status photometry_error peak_flux_erg_cm2_s counting_sigma_erg_cm2_s gross_peak_flux_erg_cm2_s "
    "scaled_peak_background_erg_cm2_s catalog_status host host_redshift host_offset_arcsec nearest_catalog "
    "nearest_name nearest_kind nearest_offset_arcsec search_random_seed seed_local_p seed_score control_reference "
    "control_exceedances_of_24 survey_fap exceedance_lo95 exceedance_hi95 event_filename event_url event_bytes "
    "event_sha256 decompressed_fits_sha256 archive_revision source_ledger source_ledger_sha256 calibrated_weights_sha256 "
    "catalog_cache review human_evidence"
).split()


def reviewed_window(ident, result, decision):
    oid, h, tag = ident
    event, source = result["event"], result["receipt"]
    vet, flux = event["vet"], event.get("flux") or {}
    host = vet.get("host") or {}
    nearest = next(iter(vet.get("matches") or []), {})
    net = flux.get("peak_flux_erg_cm2_s")
    sigma = flux.get("flux_counting_sigma_erg_cm2_s")
    gross = flux.get("gross_flux_erg_cm2_s")
    background = flux.get("scaled_background_flux_erg_cm2_s")
    quantities = (net, sigma, gross, background)
    need(all(v is None or math.isfinite(v) for v in quantities) and (sigma is None or sigma >= 0) and
         (gross is None or gross >= 0) and (background is None or background >= 0) and
         (flux.get("event_sha256") is None or flux["event_sha256"] == h), "Invalid stored flux value or input identity")
    if all(v is not None for v in quantities):
        flux_components(flux)
        status = "calibrated"
    else:
        status = "partial" if any(v is not None for v in quantities) else "unavailable"
    event_fields = ("ra_deg dec_deg position_error_arcsec event_mjd duration_s trial_duration_s ccd_id "
                    "start_met_s stop_met_s peak_start_met_s peak_window_s peak_live_s peak_counts peak_annulus_counts "
                    "on_counts annulus_counts annulus_area_ratio before_counts before_s after_counts after_s seed_local_p seed_score")
    return {
        "cohort":result["group"]["name"], "obsid":oid, "trigger_key":tag, "disposition":decision["decision"],
        **{k:event.get(k) for k in event_fields.split()}, "photometry_status":status, "photometry_error":flux.get("error"),
        "peak_flux_erg_cm2_s":net, "counting_sigma_erg_cm2_s":sigma, "gross_peak_flux_erg_cm2_s":gross,
        "scaled_peak_background_erg_cm2_s":background, "catalog_status":vet["status"], "host":host.get("name"),
        "host_redshift":host.get("redshift"), "host_offset_arcsec":host.get("separation_arcsec"),
        "nearest_catalog":nearest.get("catalog"), "nearest_name":nearest.get("name"), "nearest_kind":nearest.get("kind"),
        "nearest_offset_arcsec":nearest.get("separation_arcsec"), "search_random_seed":result["search_seed"],
        "control_reference":result["control_reference"], "control_exceedances_of_24":result["n"], "survey_fap":result["fap"],
        "exceedance_lo95":result["lo"], "exceedance_hi95":result["hi"], "event_filename":source["filename"],
        "event_url":source["url"], "event_bytes":source["bytes"], "event_sha256":h, "decompressed_fits_sha256":source["fits_sha256"],
        "archive_revision":result["revision"], "source_ledger":result["ledger"], "source_ledger_sha256":result["ledger_sha256"],
        "calibrated_weights_sha256":flux.get("calibrated_fits_sha256"), "catalog_cache":event.get("catalog_cache"),
        "review":decision["reason"], "human_evidence":json.dumps(decision["evidence"],ensure_ascii=False,separators=(",",":")),
    }


def adjudicate(manifest, leads, sig):
    need(manifest.get("reviewed_inputs_sha256") == sig, "Adjudications refer to a different measured release")
    decisions, hashes = {}, {}
    for d in manifest["decisions"]:
        k = key(d)
        need(k in leads and k not in decisions and d["decision"] in ALLOWED, "Unknown, duplicate or unsupported disposition")
        r = leads[k]
        need(isinstance(d["reason"],str) and d["reason"].strip(), "Missing written decision")
        links, sources = references(d["evidence"],k,r["lead"],r["group"]["root"])
        hashes.update(sources)
        if d["decision"] == "candidate":
            ev = r["event"]
            flux = ev.get("flux",{})
            need(r["group"]["name"] != GROUPS[2] and ev["vet"].get("errors") == {} and ev["vet"]["status"] != "catalog_incomplete"
                 and ev.get("diagnostics",{}).get("files"), "Candidate lacks completed physical or catalog characterization")
            flux_components(flux)
            directory = r["group"]["root"] / "photometry" / str(k[0]) / r["lead"]["input_name"] / k[2]
            need(flux.get("event_sha256") == k[1] and file_sha(directory/"peak_weights.fits.gz") == flux["calibrated_fits_sha256"]
                 and ev["diagnostics"].get("input_sha256") == k[1], "Candidate characterization byte identity changed")
            for name,h in ev["diagnostics"]["files"].items():
                need(Path(name).name == name and file_sha(directory/name) == h,"Candidate diagnostic bytes changed")
        decisions[k] = (d["decision"],fmt(d["reason"]),links)
    need(decisions.keys() == leads.keys(), f"{len(leads.keys()-decisions.keys())} retained real leads still require adjudication")
    return decisions, hashes


def instrument_rate_bounds(total_passes, exposure):
    """Conservative upper rates for distinct detectable events in the primary fields."""
    need(isinstance(total_passes, int) and total_passes >= 0, "Invalid instrumental pass count")
    mean_upper = float(chi2.ppf(0.95, 2 * (total_passes + 1)) / 2)
    return {str(d): mean_upper * 1e6 / exposure[str(d)] if exposure[str(d)] > 0 else None
            for d in DURATIONS}


def publish(groups, leads, decisions, raw_inputs, manual_hashes, sig, manifest, supplement):
    dated = groups[0]
    counts = {g["name"]:Counter(decisions[key(r)][0] for r in g["leads"]) for g in groups}
    sky = groups[:3]
    combined = sum((counts[g["name"]] for g in sky),Counter())
    need(dated["scanned"] > 0 and dated["exposure"]["60"] > 0 and dated["total"]["attempted"] > 0,
         "No positive primary exposure with measured completeness")
    def zero(g):
        c = counts[g["name"]]
        return (g["name"] != GROUPS[2] and g["scanned"] > 0 and g["scanned"] == g["complete"]
                and not c["candidate"] and not c["unresolved"] and g["total"]["recovered"] > 0)
    if combined["candidate"]:
        label = "Candidate"
    elif zero(dated) and not combined["unresolved"] and all(g["scanned"] == g["complete"] for g in sky):
        label = "Null result"
    else:
        label = "Measured"
    n, x = dated["scanned"], dated["exposure"]["60"]/1e6
    c = counts[PRIMARY]
    if c["candidate"]:
        headline = f"Written review retained {c['candidate']} FXT candidate{'s' if c['candidate']!=1 else ''} in {n:,} dated ObsIDs covering {x:.6f} Ms of 60-second start support."
    elif label == "Null result":
        headline = f"No FXT candidate survived in {n:,} dated ObsIDs covering {x:.6f} Ms of eligible 60-second flare-start time."
    else:
        headline = (f"No new extragalactic FXT was established in {n:,} dated ObsIDs and {x:.6f} Ms of 60-second start support. "
                    f"{c['unresolved']} lead{'s remain' if c['unresolved']!=1 else ' remains'} unresolved.")
    if combined["candidate"] > c["candidate"]:
        extra = combined["candidate"]-c["candidate"]
        headline += f" The separate release-date supplement contains {extra} candidate{'s' if extra!=1 else ''}."
    lines = ["# Chandra ACIS search for minute-to-hour X-ray flares","",f"**{headline}**","",
             f"The CXC census selected {len(dated['ids']):,} ACIS imaging ObsIDs by first public release from "
             f"2021-01-01 through 2026-11-21 inclusive. All selected products completed the blind search in {n:,} ObsIDs. "
             f"Catalog and required physical follow-up completed in {dated['complete']:,}. The other "
             f"{len(dated['ids'])-n:,} census ObsIDs contribute zero exposure to this measurement. Simultaneous CCDs "
             "and interleaved event streams count once per ObsID.", "",
             "## Search and controls", "",
             f"The seed threshold was a local binomial tail of {SEED_P:g}. The dated data produced {dated['raw']:,} raw triggers "
             f"and {dated['instrument']:,} distinct instrumental passes. The automated counterpart screen retained "
             f"{len(dated['leads']):,} real leads. Across 24 scrambled surveys the corresponding counts ranged from "
             f"{min(dated['control'])} to {max(dated['control'])}, with a mean of {sum(dated['control'])/24:.3f} per survey."]
    if dated["leads"]:
        best = max((leads[key(r)] for r in dated["leads"]),key=lambda v:v["event"]["seed_score"])
        lines.append(f"The highest retained real seed score was {best['event']['seed_score']:.6g} "
                     f"(−log₁₀ of the seed tail). {best['n']} of 24 control maxima equaled or exceeded it.")
    lines.extend(["", "FAP uses the automated control maximum per survey and the plus-one estimator, with resolution 1/25. "
                  "With incomplete follow-up it uses every instrumental pass in each control. The quoted exact 95% intervals "
                  "describe the control exceedance probability under the time-scramble model. Human review was applied to "
                  "the real leads. The measured injection efficiency covers the instrumental selection.", "",
                  "## Exposure, recovery and conditional limits"])
    from completeness_figure import render
    figure = HERE/"evidence/injection_recovery.png"
    render(dated["reports"]/"completeness.csv", figure, dated["total"]["attempted"], dated["total"]["recovered"])
    lines.extend(["", "![Measured instrumental recovery by injected photon count and flare duration]("
                  + href(figure) + ")", "",
                  "Cells show the recovered fraction and the number recovered out of those attempted. Each cell pools "
                  "the three injected PSF widths in their measured proportions. The full grid below retains each width "
                  "and its binomial interval. These efficiencies apply to injected event counts and sampled locations."])
    bounds = {}
    for g in groups:
        name, ctr = g["name"], counts[g["name"]]
        cohort = NAMES[name]
        lines.extend(["",f"### {cohort}","",
                      f"{len(g['ids']):,} assigned ObsIDs, {g['scanned']:,} fully scanned and {g['complete']:,} with "
                      f"follow-up complete. The real data produced {g['raw']:,} raw triggers and {len(g['leads'])} retained leads "
                      f"before review. Manual dispositions were {ctr['candidate']} candidate, {ctr['unresolved']} unresolved and "
                      f"{ctr['rejected']} rejected."])
        bounds[name] = {}
        if not g["scanned"]:
            continue
        tot = g["total"]
        lines.append(f"The 24 scrambled catalogs retained {sum(g['control'])} events in total, "
                     f"or {sum(g['control'])/24:.3f} per catalog. Injections recovered {tot['recovered']:,} of "
                     f"{tot['attempted']:,} attempted ({100*tot['recovered']/tot['attempted']:.1f}%)."
                     if tot["attempted"] else f"The 24 scrambled catalogs retained {sum(g['control'])} events in total. No injection was attempted.")
        if name != GROUPS[2] and g["exposure"]["60"] > 0:
            lines.append(f"The automated retained-control rate was "
                         f"{sum(g['control'])*1e6/(24*g['exposure']['60']):.4g} per Ms of 60-second start support.")
        if name == GROUPS[2]:
            lines.append("This detector-control sample contributes zero celestial exposure.")
            continue
        lines.extend(["","| Trial (s) | Eligible start time (s) | Injection recovery | 95% detected events/Ms, upper |",
                      "|---:|---:|---:|---:|"])
        for duration in DURATIONS:
            t, grid = g["exposure"][str(duration)], g["grid"][duration]
            upper = -math.log(.05)*1e6/t if zero(g) and t>0 and grid["recovered"]>0 else None
            bounds[name][str(duration)] = upper
            lines.append(f"| {duration} | {t:,.3f} | {grid['recovered']}/{grid['attempted']} | "
                         + (f"{upper:.5g}" if upper is not None else "—") + " |")
        p = g["reports"]
        lines.extend(["", " · ".join((link(p/("event_receipts.csv" if name==PRIMARY else "event_inputs.csv"),"Event URLs and SHA-256"),
                    link(p/("control_trials.csv" if name==PRIMARY else "control_counts.json"),"Control outcomes"),
                    link(p/("completeness.csv" if name==PRIMARY else "injections.json"),"Counts and PSF recovery grid")))])
    lines.extend(["","Each numeric limit is the zero-count one-sided 95% Poisson upper bound divided by that duration’s "
                  "eligible flare-start exposure. It applies to detectable events in those fields under the stated classification "
                  "and sensitivity. No solid-angle or intrinsic-population correction is applied. A dash marks an unsupported bound. "
                  "The supplemental cohorts contribute no exposure to the dated result. " + link(supplement/"release_bounds.json","Release-date evidence") + ".",
                  ""])
    instrument_bounds = instrument_rate_bounds(dated["instrument"], dated["exposure"])
    adjudication_bounds = {}
    if dated["complete"] == dated["scanned"]:
        for scope, n_events in (("accepted_candidates",counts[PRIMARY]["candidate"]),
                                ("candidates_and_unresolved",counts[PRIMARY]["candidate"]+counts[PRIMARY]["unresolved"])):
            adjudication_bounds[scope] = instrument_rate_bounds(n_events, dated["exposure"])
    lines.extend(["## Upper rates in the dated fields", "",
                  f"Counting all {dated['instrument']:,} distinct instrumental passes as possible FXTs gives a conservative "
                  "upper rate for distinct, instrument-detectable events in the dated fields. Each duration uses the full pass "
                  "count from all trial lengths, including known foreground events and ambiguous triggers. For a Poisson "
                  "process the 95% upper mean for N detected events is half the 95th percentile of a chi-square distribution "
                  "with 2(N+1) degrees of freedom. Dividing it by that duration's eligible start time bounds the mean number "
                  "of detectable events per Ms of field observation. The calculation assumes distinct events are resolved "
                  "and makes no correction for missed flares, sky area, or source-population selection."])
    if adjudication_bounds:
        lines.append(f"The review accepted {counts[PRIMARY]['candidate']} candidates and left "
                     f"{counts[PRIMARY]['unresolved']} windows unresolved. The first column of rates bounds the yield "
                     "of the completed search and manual candidate selection. The next counts every unresolved window "
                     "as a possible new FXT. Both are conditional on the counterpart and discovery criteria. The "
                     "injections measured instrumental recovery only.")
    else:
        lines.append("Incomplete counterpart or photometric follow-up prevents an upper bound for the fully reviewed selection.")
    lines.extend(["", "| Assumed duration (s) | Accepted candidates | Candidates + unresolved | All instrumental passes |",
                  "|---:|---:|---:|---:|",])
    for d in DURATIONS:
        ds = str(d)
        values = [adjudication_bounds.get(scope,{}).get(ds) for scope in ("accepted_candidates","candidates_and_unresolved")]
        values.append(instrument_bounds[ds])
        lines.append(f"| {d} | " + " | ".join(f"{v:.5g}" if v is not None else "—" for v in values) + " |")
    lines.append("All rates in this table are one-sided 95% Poisson upper limits per million seconds of eligible flare-start time.")
    unknown = groups[2]
    if unknown["scanned"] and unknown["exposure"]["60"] > 0:
        certain = dated["exposure"]["60"] + groups[1]["exposure"]["60"]
        possible = certain + unknown["exposure"]["60"]
        lines.extend(["", f"Across the measured sky observations, assigning none or all of the {unknown['scanned']} "
                      "scanned observations with unresolved release dates to the window gives "
                      f"{certain/1e6:.6f}–{possible/1e6:.6f} Ms of eligible 60-second start time. "
                      "This interval describes release-date uncertainty in measured coverage."])
    candidate_file = HERE/"evidence/fxt_candidates.csv"
    reviewed_file = HERE/"evidence/reviewed_fxt_triggers.csv"
    reviewed_rows = []
    written = {key(row):row for row in read(manifest)["decisions"]}
    need(written.keys() == decisions.keys() == leads.keys(), "Written decisions and publication identity sets differ")
    lines.extend(["", "## Manual review", "", link(reviewed_file,"Every manually reviewed trigger, disposition and measurement")
                  + " · " + link(candidate_file,"New FXT candidates with exact inputs and measurements"), "",
                  "Both tables use the same columns. Blank flux fields indicate unavailable measurements. "
                  "The human evidence column contains a JSON list of source links or survey-relative paths and section anchors."])
    signed = [v["event"]["flux"]["peak_flux_erg_cm2_s"] for v in leads.values()
              if v["event"].get("flux",{}).get("peak_flux_erg_cm2_s") is not None]
    if signed:
        lines.extend(["",f"Of {len(signed):,} calibrated retained windows, {sum(f <= 0 for f in signed):,} had a "
                      "nonpositive background-subtracted 60-second peak flux estimate. An annulus photon with a large "
                      "response weight can dominate this estimate. The blind search tests photon counts. The quoted "
                      "counting errors condition on the selected minute and omit response uncertainty."])
    if not leads:
        lines.extend(["","The automated selection retained no real trigger for manual review."])
    for ident in sorted(leads,key=lambda k:(("candidate","unresolved","rejected").index(decisions[k][0]),k)):
        r, (state,reason,refs) = leads[ident], decisions[ident]
        oid, h, trig = ident
        ev, source = r["event"], r["receipt"]
        note = written[ident]
        need(state == note["decision"] and fmt(note["reason"]) == reason, "Written disposition changed during publication")
        reviewed_rows.append(reviewed_window(ident,r,note))
        if state != "candidate":
            lines.extend(["",f"- **{state.capitalize()} · ObsID {oid}, {trig}.** {reason} {refs}"])
            continue
        host, matches = ev["vet"]["host"], ev["vet"].get("matches",[])
        counterpart = (f"Proposed host {fmt(host['name'])} at {fmt(host['separation_arcsec'])}″, catalog redshift {fmt(host.get('redshift'))}."
                       if host else (f"No identified host. Nearest matched counterpart {fmt(matches[0]['name'])} at "
                                     f"{fmt(matches[0]['separation_arcsec'])}″." if matches else
                                     "No catalog counterpart within the screening region."))
        f = ev["flux"]
        gross, background = flux_components(f)
        flux_note = (" The background subtraction leaves the source energy flux unconstrained in this minute."
                     if f["peak_flux_erg_cm2_s"] <= 0 else "")
        lines.extend(["",f"### Candidate · ObsID {oid}, {trig}","",
                 f"{NAMES[r['group']['name']]}. ICRS position ({ev['ra_deg']}, {ev['dec_deg']})°, "
                 f"90% radius {ev['position_error_arcsec']}″. Peak MJD {ev['event_mjd']} UTC. Aperture-photon "
                 f"t95−t05 duration {ev['duration_s']} s in a {ev['trial_duration_s']} s search window. "
                 f"Signed peak 0.5–7 keV aperture flux estimate {f['peak_flux_erg_cm2_s']} ± "
                 f"{f['flux_counting_sigma_erg_cm2_s']} erg cm⁻² s⁻¹ (counting uncertainty). Gross aperture and scaled "
                 f"background estimates are {gross} and {background} erg cm⁻² s⁻¹.{flux_note} {counterpart}", "",
                 f"The aperture contained {ev['before_counts']} photons in {ev['before_s']:.1f} s of preceding good time, "
                 f"{ev['on_counts']} in the search window and {ev['after_counts']} in {ev['after_s']:.1f} s afterward. "
                 f"The surrounding annulus contained {ev['annulus_counts']} photons during the window.", "",
                 f"The automated survey FAP was {r['fap']:.3g} ({r['n']}/24 control maxima). The exact 95% interval "
                 f"for that exceedance probability is {r['lo']:.4g}–{r['hi']:.4g}. {reason} {refs}", "",
                 f"Input {source['bytes']} bytes, SHA-256 `{h}`. [{fmt(source['filename'])}]({source['url']})"])
    candidate_rows = [r for r in reviewed_rows if r["disposition"] == "candidate"]
    expected = {(g,d):n for g,totals in counts.items() for d,n in totals.items()}
    columns=set(REVIEW_FIELDS)
    need(all(r.keys()==columns for r in reviewed_rows) and
         {(r["obsid"],r["event_sha256"],r["trigger_key"]) for r in reviewed_rows}==leads.keys() and
         Counter((r["cohort"],r["disposition"]) for r in reviewed_rows)==expected and
         len(candidate_rows)==combined["candidate"], "Published table counts and cohort dispositions disagree")
    report.rows_to_csv(reviewed_file, REVIEW_FIELDS, reviewed_rows)
    report.rows_to_csv(candidate_file, REVIEW_FIELDS, candidate_rows)
    lines.extend(["","## Reproduction","",
                  "The archived EVT2 manifests bind the compressed and decompressed bytes to each observation. The full ledgers retain "
                  "all detector vetoes, catalog query receipts, eligible observing intervals and controls. "
                  + link(HERE/"evidence/obsid_coverage.csv","Primary coverage") + " · "
                  + link(HERE/"evidence/leads.json","All automated primary leads") + " · " + link(manifest,"Written adjudications")
                  + " · " + link(HERE/"METHODS.md","Frozen search, calibration and replay commands") + ".", "",
                  "Rebuild the primary report, then render this publication in the Chandra workflow.", "", "```bash",
                  "python results/chandra_acis_fxt/report.py",
                  "python results/chandra_acis_fxt/publication.py --supplement " + href(supplement) +
                  " --decisions " + href(manifest), "```",""])
    supporting = [(title,p) for title,p in (("Per-lead review",HERE/"evidence/LEAD_REVIEW.md"),
                         ("Independent visits",HERE/"evidence/independent_visits/independent_visits.json"),
                         ("Contemporary notices",HERE/"evidence/alerts/review.json"),
                         ("Matched Gaia stellar, galaxy and quasar classifications",HERE/"evidence/gaia_dr3_properties/README.md"),
                         ("S4 frame occupancy in real data and all 24 scrambles",HERE/"evidence/frame_forensics/population/final/README.md"),
                         ("Field and detector measurements for ObsIDs 23324 and 32380",HERE/"PRIORITY.md")) if p.is_file()]
    field = HERE/"evidence/field_followup/27089/45008126305-1"
    for name,title in (("photon_check.json","GW170817 field · ObsID 27089 photon comparison"),
                       ("history/history.json","GW170817 field · all-years visits and exclusions")):
        if (field/name).is_file():
            supporting.append((title,field/name))
    for doc in sorted(field.glob("source_replay/*/README.md")):
        if (doc.parent/"manifest.json").is_file():
            supporting.append(("GW170817 field · time-averaged rates and original GCN circulars",doc))
    for history in sorted((HERE/"evidence/pointed_visits").glob("**/history/history.json")):
        base=history.parent.parent
        if (base/"execution.json").is_file() and (base/"leads.json").is_file():
            ids=", ".join(str(x["obsid"]) for x in read(base/"leads.json"))
            supporting.append((f"All-mission independent aperture measurements for ObsIDs {ids}",history))
            for energy in sorted(base.glob("*/**/ENERGY.md")):
                supporting.append((f"Source and annular energy counts for ObsID {energy.relative_to(base).parts[0]}",energy))
    for result in sorted((HERE/"evidence/frame_forensics").glob("*/*/real/*/*/result.json")):
        doc = result.parent/"README.md"
        item = read(result)
        if doc.is_file() and item.get("status")=="complete" and (item.get("measurement") or {}).get("instrument_ccd")==8:
            target = item["target"]
            supporting.append((f"ACIS-S4 frame measurement · ObsID {int(target['obsid'])}, {target['trigger']}",doc))
    for doc in sorted((HERE/"evidence/csc_deblend").glob("*/*/*/README.md")):
        if doc.parent.name.isdecimal() and (doc.parent/"measurement.json").is_file():
            supporting.append((f"Trial and CSC neighbor photon maps for ObsID {doc.parent.name}",doc))
    if supporting:
        lines.extend(["## Supporting measurements","",*(f"- {link(path,title)}" for title,path in supporting),""])
    outcome = {"label":label, "headline":headline,"reviewed_inputs_sha256":sig,"input_sha256":raw_inputs,
               "decisions_sha256":file_sha(manifest),"candidate_table_sha256":hashlib.sha256(candidate_file.read_bytes()).hexdigest(),
               "reviewed_table_sha256":hashlib.sha256(reviewed_file.read_bytes()).hexdigest(),
               "reviewed_rows":len(reviewed_rows),"candidate_rows":len(candidate_rows),
               "supporting_measurements_sha256":{href(path):file_sha(path) for _,path in supporting},
               "recovery_figure_sha256":hashlib.sha256(figure.read_bytes()).hexdigest(),"cited_local_evidence_sha256":manual_hashes,"dispositions":counts,
               "conditional_upper95_detected_events_per_Ms":bounds,
               "review_selection_upper95_events_per_Ms":adjudication_bounds,
               "instrument_pass_conservative_upper95_detected_events_per_Ms":instrument_bounds, "implementation_sha256":file_sha(Path(__file__))}
    # update_catalog checks the existing bullet and rewrites only line 9 and this entry.
    report.update_catalog(headline,label)
    original = (report.ROOT/"README.md").read_text().splitlines()
    bullet, = [s for s in original if s.startswith("- ") and report.TAG in s]
    kinds = "Candidate|Measured|Null result|Correction|Prediction"
    labels = re.compile(rf"(?<!\*)\*({kinds})\*(?!\*)")
    need(labels.findall(bullet) == [label] and len(original)>8, "Catalog update did not preserve a single counted entry")
    tally = Counter(k for line in original if line.startswith("- ") for k in labels.findall(line))
    for phrase, value in (("entries",sum(tally.values())),("new measurements",tally["Measured"]),
                          ("candidate discoveries",tally["Candidate"]),("forecasts",tally["Prediction"]),
                          ("corrections",tally["Correction"]),("carefully recorded null",tally["Null result"])):
        need(re.search(rf"\b{value} {phrase}\b", original[8]), "Catalog summary count mismatch")
    dump(HERE/"evidence/publication.json",outcome)
    path = HERE/"README.md"
    temp = path.with_suffix(".md.tmp")
    temp.write_text("\n".join(lines))
    temp.replace(path)
    return outcome


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--supplement",required=True,help="Survey-relative directory of the deposited supplemental run")
    choice = p.add_mutually_exclusive_group(required=True)
    choice.add_argument("--signature",action="store_true")
    choice.add_argument("--decisions",help="Survey-relative JSON with reviewed_inputs_sha256 and decisions")
    opt = p.parse_args()
    need(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")=="dovvnloading/Open-Space" and
         os.getenv("GITHUB_REF")=="refs/heads/research/chandra-acis-fast-transients","Use the Chandra research workflow")
    supplement = (HERE/opt.supplement).resolve()
    groups, leads, sig, files = load(supplement)
    if opt.signature:
        print(sig)
        return
    manifest = local(opt.decisions)
    decisions, sources = adjudicate(read(manifest),leads,sig)
    result = publish(groups,leads,decisions,files,sources,sig,manifest,supplement)
    print(result["headline"])


if __name__ == "__main__":
    main()
