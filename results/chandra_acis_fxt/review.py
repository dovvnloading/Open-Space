"""Build an offline human-review index from a deposited Chandra survey report."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
from urllib.parse import quote, urlsplit


def fmt(value, digits=5):
    if value is None or value == "":
        return "—"
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    if isinstance(value, bool):
        return "yes" if value else "no"
    value = html.escape(" ".join(str(value).split()), quote=False)
    return re.sub(r"([\\|`*\[\]<>])", r"\\\1", value)


def counts(items):
    return ", ".join(f"{fmt(k)} {v}" for k, v in sorted(Counter(items).items())) or "none"


def table(lines, headers, rows):
    lines.extend(["", "| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)])
    for row in rows:
        lines.append("| " + " | ".join(row) + " |")


class Deposit:
    def __init__(self, root):
        self.root = root.resolve()
        self.gaps = []

    def path(self, relative):
        path = (self.root / relative).resolve()
        if not path.is_relative_to(self.root):
            raise ValueError("path escapes the evidence directory")
        return path

    def link(self, relative, label):
        try:
            path = self.path(relative)
            if path.is_file():
                target = quote(path.relative_to(self.root).as_posix(), safe="/")
                return f"[{fmt(label)}]({target})"
        except (OSError, ValueError):
            pass
        return f"{fmt(label)} (file unavailable)"

    def read(self, relative, kind):
        try:
            path = self.path(relative)
            if kind == "csv":
                with path.open(encoding="utf-8", newline="") as source:
                    return list(csv.DictReader(source))
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rt", encoding="utf-8") as source:
                result = json.load(source)
            if not isinstance(result, kind):
                raise ValueError("unexpected top-level schema")
            return result
        except (OSError, EOFError, UnicodeError, ValueError, csv.Error) as exc:
            self.gaps.append(f"{relative} ({type(exc).__name__})")
            return None


def key(row):
    return "/".join(str(row.get(k, "")) for k in ("obsid", "input_sha256", "trigger_key"))


def raw_events(depot, coverage, leads):
    """Read each required ledger once, keeping only retained real trigger details."""
    wanted = {key(row) for row in leads}
    included = {str(row["obsid"]) for row in coverage if row.get("scanned_all_products", "").lower() == "true"}
    ids = included | {str(row.get("obsid", "")) for row in leads}
    statuses, vetoes, found = defaultdict(Counter), Counter(), {}
    read_count = 0
    for oid in sorted(ids):
        if not oid.isdigit():
            depot.gaps.append(f"Invalid ObsID {oid}")
            continue
        path = f"observations/{oid}.json.gz"
        ledger = depot.read(path, dict)
        if ledger is None:
            continue
        if str(ledger.get("obsid")) != oid or not isinstance(ledger.get("products"), list):
            depot.gaps.append(f"{path} (invalid identity or products)")
            continue
        try:
            local_statuses, local_vetoes = defaultdict(Counter), Counter()
            for product in ledger["products"]:
                src, search = product["input"], product["search"]
                streams = [("real", search["triggers"])]
                if oid in included:
                    streams.extend(("scrambled", replica["triggers"]) for replica in search["controls"])
                for stream, events in streams:
                    for event in events:
                        stage = "pass" if event["instrument_pass"] else "veto"
                        status = event["vet"]["status"]
                        local_statuses[status][(stream, stage)] += 1
                        if stage == "veto":
                            local_vetoes.update((stream, reason) for reason in event.get("vetoes", []))
                        identity = f"{oid}/{src['sha256']}/{event.get('trigger_key', '')}"
                        if stream == "real" and identity in wanted:
                            found[identity] = (event, src, event.get("catalog_cache", ledger.get("catalog_cache")))
            if oid in included:
                read_count += 1
                for status, tally in local_statuses.items():
                    statuses[status].update(tally)
                vetoes.update(local_vetoes)
        except (KeyError, TypeError, AttributeError) as exc:
            depot.gaps.append(f"{path} (incomplete schema, {type(exc).__name__})")
    return statuses, vetoes, found, read_count, len(included)


def survey(lines, depot, summary, coverage, controls, recovery, raw):
    lines.extend(["## Survey and controls", "",
                  depot.link("summary.json", "Survey summary") + " · " + depot.link("obsid_coverage.csv", "Observation coverage")
                  + " · " + depot.link("event_receipts.csv", "Input receipts") + ".", "",
                  f"Inventory ObsIDs {fmt(summary.get('inventory_obsids'))}. Fully scanned {fmt(summary.get('scanned_obsids'))}. "
                  f"Follow-up complete {fmt(summary.get('completed_obsids'))}. "
                  f"Release window complete {fmt(summary.get('release_window_complete'))}. "
                  f"Inventory exceptions {fmt(summary.get('inventory_exceptions'))}.",
                  "", "Recorded coverage by ObsID status. " + (counts(r.get("status", "unknown") for r in coverage)
                                                               if coverage is not None else "unavailable") + "."])
    exposure = summary.get("exposure_by_duration_s") or {}
    if exposure:
        table(lines, ["Trial duration (s)", "Eligible field-start live time (s)"],
              ((fmt(duration), fmt(seconds, 12)) for duration, seconds in sorted(exposure.items(), key=lambda p: float(p[0]))))
    lines.extend(["", f"Real raw triggers {fmt(summary.get('real_raw_triggers'))}. Distinct instrument passes "
                  f"{fmt(summary.get('real_instrument_pass'))}. Distinct retained real leads "
                  f"{fmt(summary.get('scientific_leads'))}."])
    statuses, vetoes, _, nread, included = raw
    lines.extend(["", f"Raw trigger statuses across {nread} of {included} ObsIDs marked fully scanned in the coverage table. "
                  "Scrambled columns pool all replicas. Multiple products can trigger on the same burst."])
    if statuses:
        columns = [("real", "pass"), ("real", "veto"), ("scrambled", "pass"), ("scrambled", "veto")]
        table(lines, ["Catalog status", "Real pass", "Real veto", "Scrambled pass", "Scrambled veto"],
              ([fmt(status)] + [str(count[c]) for c in columns] for status, count in sorted(statuses.items())))
    if vetoes:
        reasons = sorted({reason for _, reason in vetoes})
        table(lines, ["Detector veto (reasons can overlap)", "Real raw", "Scrambled raw"],
              ((fmt(reason), str(vetoes["real", reason]), str(vetoes["scrambled", reason])) for reason in reasons))
    lines.extend(["", depot.link("control_trials.csv", "Per-replica control counts and maxima") + ". "
                  f"Replicas {fmt(summary.get('control_replicas'))}. Mean retained per replica "
                  f"{fmt(summary.get('mean_control_survivors'))}. Retained control events per Ms "
                  f"{fmt(summary.get('control_events_per_Ms'))}."])
    if controls:
        stages = [("triggers", "Raw triggers"), ("instrument_pass", "Distinct instrument passes"),
                  ("surviving_or_unresolved", "Distinct retained"), ("incomplete_queries", "Raw triggers with incomplete catalogs")]
        values = []
        for field, label in stages:
            try:
                v = [int(c[field]) for c in controls]
                values.append([label, str(min(v)), str(max(v)), fmt(sum(v) / len(v))])
            except (ValueError, KeyError, TypeError):
                depot.gaps.append(f"control_trials.csv (invalid {field})")
        table(lines, ["Control selection", "Minimum per replica", "Maximum", "Mean"], values)
    lines.extend(["", depot.link("completeness.csv", "Injection recovery by duration, counts and PSF width") + ". "
                  f"Scheduled {fmt(summary.get('scheduled_injections'))}. Attempted {fmt(summary.get('attempted_injections'))}. "
                  f"Recovered {fmt(summary.get('recovered_injections'))}. Recovery measures the instrumental selection "
                  "on the sampled photon-count and width grid."])
    if recovery:
        grid = defaultdict(Counter)
        for cell in recovery:
            try:
                row = {k: int(cell[k]) for k in ("scheduled", "attempted", "recovered")}
                grid[int(cell["duration_s"])].update(row)
            except (KeyError, ValueError, TypeError):
                depot.gaps.append("completeness.csv (invalid cell)")
        table(lines, ["Injected duration (s)", "Scheduled", "Attempted", "Recovered"],
              ([str(d)] + [str(c[k]) for k in ("scheduled", "attempted", "recovered")] for d,c in sorted(grid.items())))


def check_followup(lines, depot, identity, visits, alert, alerts):
    vpath = "independent_visits/independent_visits.json"
    rows = (visits.get("visits") or {}).get(identity)
    if rows is None:
        lines.append("- Independent visits. No result indexed for this lead.")
    elif not rows:
        lines.append("- Independent visits. No eligible nearby product listed. See the " +
                     depot.link(vpath, "receipt selection and its exceptions") + ".")
    else:
        lines.append("- Independent visits. " + counts(v.get("status", "unknown") for v in rows) + ". "
                     + depot.link(vpath, "All visit records") + ".")
        measured = [v for v in rows if isinstance(v.get("per_lead_bonferroni_p"), (float, int))]
        if measured:
            v = min(measured, key=lambda r: r["per_lead_bonferroni_p"])
            lines.append(f"  - Smallest per-lead corrected excess probability {fmt(v.get('per_lead_bonferroni_p'))} "
                         f"in {fmt(v.get('temporal_relation'))} ObsID {fmt(v.get('visit_obsid'))}. "
                         f"Conditional 95% upper rate {fmt(v.get('conditional_aperture_rate_upper95_ct_s'))} ct/s. "
                         + depot.link("independent_visits/" + str(v.get("record", "")), "Visit measurement and input receipt") + ".")
        failed = next((v for v in rows if v.get("status") == "error"), None)
        if failed:
            lines.append("  - " + depot.link("independent_visits/" + str(failed.get("record", "")), "First failed visit receipt") + ".")
    if alert is None:
        lines.append(f"- GCN circulars. No current result indexed for this lead. Archive status {fmt(alerts.get('status'))}.")
    else:
        matches = alert.get("matches") or []
        lines.append(f"- GCN circulars. {fmt(alert.get('status'))}. "
                     + counts(m.get("relation", "unknown") for m in matches) + ". "
                     + depot.link("alerts/review.json", "Full GCN evidence and parser coverage") + ".")
        relevant = [m for m in matches if str(m.get("relation", "")).startswith("position")]
        for m in relevant[:4]:
            positions = [p["offset_arcsec"] for p in m.get("positions", []) if p.get("offset_arcsec") is not None]
            lines.append(f"  - {fmt(m.get('relation'))}, closest offset {fmt(min(positions) if positions else None)}″. "
                         + depot.link("alerts/" + str(m.get("file", "")), f"GCN {m.get('id')}") + ".")
        if len(relevant) > 4:
            lines.append(f"  - {len(relevant)-4} additional positional matches appear in the full GCN evidence.")


def lead_section(lines, depot, lead, raw, control_n, visits, alert, alerts):
    oid, name, trigger = (str(lead.get(k, "")) for k in ("obsid", "input_name", "trigger_key"))
    event, src, cache = raw or ({}, {}, None)
    vet = event.get("vet") or {}
    digest = lead.get("input_sha256")
    lines.extend(["", f"### ObsID {fmt(oid)} · {fmt(trigger)} · {fmt(name)}", "",
                  f"- Deposited disposition **{fmt(lead.get('status'))}**. {fmt(vet.get('reason'))}.",
                  f"- ICRS RA {fmt(lead.get('ra_deg'),12)}°, Dec {fmt(lead.get('dec_deg'),12)}°, "
                  f"position error {fmt(lead.get('position_error_arcsec'))}″. Peak MJD (UTC) {fmt(lead.get('event_mjd'),15)}.",
                  f"- Duration from aperture photon t95−t05 {fmt(lead.get('duration_s'))} s, "
                  f"trial {fmt(lead.get('trial_duration_s'))} s. Peak live time {fmt(lead.get('peak_live_s'))} s, "
                  f"rate {fmt(lead.get('peak_rate_ct_s'))} ct/s. Peak aperture flux (0.5–7 keV) "
                  f"{fmt(lead.get('peak_flux_erg_cm2_s'))} erg cm⁻² s⁻¹, counting σ {fmt(lead.get('counting_sigma'))} erg cm⁻² s⁻¹.",
                  f"- Catalog host {fmt(lead.get('host'))}, offset {fmt(lead.get('host_offset_arcsec'))}″, "
                  f"redshift {fmt(lead.get('host_z'))}, chance-projection upper bound {fmt(lead.get('host_chance_upper95'))}.",
                  f"- Survey false-alarm estimate {fmt(lead.get('survey_fap'))}, from "
                  f"{fmt(lead.get('control_exceedances'))}/{fmt(control_n)} maxima at or above the real seed score "
                  f"{fmt(event.get('seed_score'))}. Reference column {fmt(lead.get('null_reference'))}. "
                  f"95% interval for exceedance probability [{fmt(lead.get('fap_lo95'))}, {fmt(lead.get('fap_hi95'))}]."])
    url = src.get("resolved_url") or lead.get("input_url") or ""
    online = f"[Archive bytes]({quote(url, safe=':/?&=%+')})" if url.startswith(("https://","http://")) else "Archive URL unavailable"
    lines.append(f"- {online}, {fmt(src.get('bytes'))} bytes. Compressed input SHA-256 `{fmt(digest)}`. "
                 f"Uncompressed FITS SHA-256 `{fmt(src.get('fits_sha256'))}`. "
                 + depot.link(f"observations/{oid}.json.gz", "ObsID ledger") + f", raw key `{fmt(trigger)}`.")
    if raw is None:
        lines.append("- Matching raw ledger entry unavailable. Aperture counts and catalog decisions cannot be inspected here.")
    else:
        lines.extend([f"- Source/annulus counts {fmt(event.get('on_counts'))}/{fmt(event.get('annulus_counts'))}, "
                      f"area ratio {fmt(event.get('annulus_area_ratio'))}. Before/after source counts "
                      f"{fmt(event.get('before_counts'))}/{fmt(event.get('after_counts'))} over "
                      f"{fmt(event.get('before_s'))}/{fmt(event.get('after_s'))} s. Immediately adjacent flank counts "
                      f"{fmt(event.get('adjacent_before_counts'))}/{fmt(event.get('adjacent_after_counts'))}. "
                      f"Peak source/annulus counts {fmt(event.get('peak_counts'))}/{fmt(event.get('peak_annulus_counts'))}.",
                      f"- CCD {fmt(event.get('ccd_id'))}, detector pixels {fmt(event.get('distinct_detector_pixels'))}, "
                      f"distinct frames {fmt(event.get('distinct_frames'))}. Detector veto flags "
                      f"{', '.join(fmt(x) for x in event.get('vetoes', [])) or 'none recorded'}. "
                      f"Galactic latitude {fmt((vet.get('position') or {}).get('galactic_latitude_deg'))}°."])
        matches = vet.get("matches") or []
        lines.append("- Catalog matches by type. " + (counts(m.get("kind", "unknown") for m in matches)
                                                      if "matches" in vet else "unavailable") + ".")
        for m in matches[:5]:
            lines.append(f"  - {fmt(m.get('catalog'))} {fmt(m.get('name'))}, {fmt(m.get('kind'))}. "
                         f"Offset {fmt(m.get('separation_arcsec'))}″, chance upper bound {fmt(m.get('chance_upper95'))}, "
                         f"redshift {fmt(m.get('redshift'))}.")
        if len(matches) > 5:
            lines.append(f"  - The ledger lists {len(matches)-5} additional matches.")
        errors = vet.get("errors") or {}
        if errors:
            lines.append("- Catalog failures. " + ". ".join(f"{fmt(k)} — {fmt(v)}" for k, v in errors.items()) + ".")
        if cache:
            bodies = [depot.link(Path(cache) / r["body_file"], label)
                      for label, r in (vet.get("queries") or {}).items() if isinstance(r,dict) and r.get("body_file")]
            if bodies:
                lines.append("- Catalog response bytes. " + " · ".join(bodies) + ". The ledger records individual SHA-256 values.")
    folder = Path("photometry") / oid / name / trigger
    assets = [depot.link(folder / p, label) for p,label in (("diagnostics.png","Detector plot"),("event_rows.csv.gz","Event rows"),
              ("lightcurves.csv","Light curves"),("flux.json","Flux receipt"))]
    lines.append("- " + " · ".join(assets) + ".")
    for field in ("diagnostics","flux"):
        if (event.get(field) or {}).get("error"):
            lines.append(f"- {field.title()} failure. {fmt(event[field]['error'])}.")
    nearby = lead.get("other_trigger_obsids")
    if nearby:
        lines.append("- Other retained bursts with overlapping positional uncertainties. ObsIDs " + ", ".join(map(fmt, nearby)) + ".")
    check_followup(lines, depot, key(lead), visits, alert, alerts)


def render(root):
    depot = Deposit(root)
    summary = depot.read("summary.json", dict) or {}
    lead_doc = depot.read("leads.json", list)
    leads = [r for r in (lead_doc or []) if isinstance(r, dict)]
    if lead_doc is not None and len(lead_doc) != len(leads):
        depot.gaps.append("leads.json (non-object lead entries)")
    coverage = depot.read("obsid_coverage.csv", "csv")
    controls = depot.read("control_trials.csv", "csv")
    recovery = depot.read("completeness.csv", "csv")
    visits = depot.read("independent_visits/independent_visits.json", dict) or {}
    alerts = depot.read("alerts/review.json", dict) or {}
    current_alerts = {}
    if lead_doc is not None and alerts:
        digest = hashlib.sha256(depot.path("leads.json").read_bytes()).hexdigest()
        if digest != alerts.get("leads_sha256"):
            depot.gaps.append("alerts/review.json (lead-file SHA-256 mismatch)")
        else:
            current_alerts = {key(row):row for row in alerts.get("leads", [])}
    raw = raw_events(depot, coverage or [], leads)
    lines = ["# Chandra FXT lead review", "",
             "Each entry retains the survey disposition and links to the observation, catalog and follow-up evidence. "
             "A dash marks an unavailable value. The survey false-alarm estimate is (exceedances + 1)/(replicas + 1)."]
    survey(lines, depot, summary, coverage, controls, recovery, raw)
    if visits:
        lines.extend(["", f"Independent visits submitted {fmt(visits.get('leads_submitted'))} leads. "
                      f"Receipt-selection exceptions {len(visits.get('unassessable_receipts') or [])}. "
                      + depot.link("independent_visits/independent_visits.json", "Coverage conditions and exceptions") + "."])
    if alerts:
        lines.extend(["", f"GCN archive status {fmt(alerts.get('status'))}. Decoded {fmt(alerts.get('decoded_circulars'))} circulars. "
                      f"Einstein Probe documents with parsed coordinates {fmt(alerts.get('ep_with_parsed_coordinates'))} of "
                      f"{fmt(alerts.get('ep_documents'))}. " + depot.link("alerts/review.json", "Alert archive receipts and limitations") + "."])
        if alerts.get("error"):
            lines.append("Archive error. " + fmt(alerts["error"]) + ".")
    lines.extend(["", "## Retained leads", "", depot.link("leads.json", "Deposited distinct-lead list") + "."])
    if lead_doc is None:
        lines.extend(["", "The distinct-lead list is unavailable."])
    elif not lead_doc:
        lines.extend(["", "The deposited list contains zero retained leads."])
    if summary.get("scientific_leads") is not None and lead_doc is not None and len(leads) != summary["scientific_leads"]:
        depot.gaps.append("leads.json (length differs from summary.json)")
    if len({key(r) for r in leads}) != len(leads):
        depot.gaps.append("leads.json (duplicate lead keys)")
    for lead in leads:
        lead_section(lines, depot, lead, raw[2].get(key(lead)), summary.get("control_replicas"),
                     visits, current_alerts.get(key(lead)), alerts)
    if depot.gaps:
        lines.extend(["", "## Evidence gaps", ""])
        for gap in sorted(set(depot.gaps)):
            lines.append("- " + fmt(gap))
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=Path(__file__).resolve().parent / "evidence")
    root = parser.parse_args().evidence
    root.mkdir(parents=True, exist_ok=True)
    result = root / "LEAD_REVIEW.md"
    tmp = result.with_suffix(".md.tmp")
    tmp.write_text(render(root), encoding="utf-8")
    tmp.replace(result)
    print(result)
