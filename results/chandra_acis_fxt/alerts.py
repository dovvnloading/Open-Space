"""Shortlist existing reports for independent reading after the blind search.

python results/chandra_acis_fxt/alerts.py [--replay]
Raw inputs under evidence/alerts/source are needed for replay. A report can
quote a field center, unrelated object, or several times. No match from this
program by itself identifies or excludes a counterpart.
"""
import argparse
from datetime import date, datetime, timedelta, timezone
import hashlib
from html.parser import HTMLParser
import io
import json
import os
from pathlib import Path
import re
import tarfile
from urllib.parse import urljoin, urlparse

from astropy.time import Time
from archive import dump
from vet import CAMPAIGN_REF, read_input, separation

HERE = Path(__file__).resolve().parent
INDEX = "https://gcn.nasa.gov/circulars"
UTC = timezone.utc
MAX_OFFSET = 360.0  # Parsed compact localizations within six arcminutes.
EP = re.compile(r"\b(?:Einstein[ -]Probe|WXT|EP[- /]?(?:WXT|FXT|Alert|trigger|\d{6}[A-Za-z]?))\b", re.I)
EVENT_DATE = re.compile(r"\b(?:EP(?:\s+Alert)?|GRB)[ -]?(\d{6})(?:[A-Za-z]|\.\d+)?\b", re.I)
DECIMAL = r"[+-]?\d{1,3}(?:\.\d+)?"
ANGLE = rf"(?:[+-]?\d{{1,3}}:\d{{2}}:\d{{2}}(?:\.\d+)?|{DECIMAL})"
RA = r"\bR\.?\s*A\.?"
DEC = r"D\.?\s*E\.?\s*C\.?"
OPT_FRAME = r"(?:\s*\(J2000(?:\.0)?\))?"
DEG = r"(?:\s*(?:deg(?:rees?)?\.?|°))?"
LABELLED = re.compile(rf"{RA}{OPT_FRAME}\s*[:=]\s*({ANGLE}){DEG}\s*[,;]?\s*{DEC}{OPT_FRAME}"
                      rf"\s*[:=]\s*({ANGLE})(?![\d.:])", re.I)
PAIR = re.compile(rf"\(?{RA}\s*,\s*{DEC}\)?{OPT_FRAME}\s*[:=]\s*\(?\s*"
                  rf"({ANGLE}){DEG}\s*,\s*({ANGLE})(?![\d.:])", re.I)
RADIUS = re.compile(r"(?:uncertainty|error(?:\s+radius)?|radius)\s*(?:(?:of|is)\s+)?"
                    r"([0-9]+(?:\.[0-9]+)?)\s*(arcsec(?:onds?)?|arcmin(?:utes?)?|deg(?:rees?)?)\b", re.I)
ONSET = re.compile(r"\b(?:T[_ ]?0|triggered|trigger\s+time|onset|started)\b"
                   r"[^.;\n]{0,65}?(20\d\d-\d\d-\d\d[T ]\d\d:\d\d:\d\d(?:\.\d+)?)"
                   r"(?:\s*\(?UTC\)?|Z)", re.I)


class ArchiveLink(HTMLParser):
    def __init__(self):
        super().__init__()
        self.href = None
        self.label = ""
        self.links = []

    def handle_starttag(self, tag, attrs):
        if tag == "a":
            self.href = dict(attrs).get("href")
            self.label = ""

    def handle_data(self, data):
        if self.href:
            self.label += data

    def handle_endtag(self, tag):
        if tag == "a" and self.href:
            if re.match(r"^JSON\s*\(", self.label.strip(), re.I):
                self.links.append(urljoin(INDEX, self.href))
            self.href = None


def archive_link(raw):
    page = ArchiveLink()
    page.feed(raw.decode("utf-8"))
    choices = set(page.links)
    if len(choices) != 1:
        raise ValueError("NASA index lacks a unique JSON archive link")
    url = choices.pop()
    loc = urlparse(url)
    host = loc.hostname or ""
    if loc.scheme != "https" or not any(host == d or host.endswith("." + d)
                                         for d in ("nasa.gov", "amazonaws.com", "cloudfront.net")):
        raise ValueError("Unexpected NASA archive host")
    return url


def angle(text, ra):
    if ":" not in text:
        return float(text)
    degrees, minutes, seconds = map(float, text.split(":"))
    if not (0 <= minutes < 60 and 0 <= seconds < 60):
        raise ValueError("invalid sexagesimal components")
    return (-1 if text.startswith("-") else 1) * (abs(degrees) + minutes/60 + seconds/3600) * (15 if ra else 1)


def positions(body):
    body = " ".join(body.replace("−", "-").split())
    for pattern in (LABELLED, PAIR):
        for m in pattern.finditer(body):
            try:
                ra, dec = angle(m[1], True), angle(m[2], False)
            except ValueError:
                continue
            if not (0 <= ra < 360 and -90 <= dec <= 90):
                continue
            radius = RADIUS.search(body[m.end(): m.end()+145])
            error = None
            if radius:
                unit = radius[2].lower()
                error = float(radius[1]) * (3600 if unit.startswith("deg") else 60 if unit.startswith("arcmin") else 1)
            yield {"ra_deg": ra, "dec_deg": dec, "nearby_quoted_radius_arcsec": error,
                   "coordinate_text": body[m.start():m.end()+145]}


def dates_in_name(text):
    days = set()
    for year_month_day in EVENT_DATE.findall(text):
        try:
            days.add(date.fromisoformat(f"20{year_month_day[:2]}-{year_month_day[2:4]}-{year_month_day[4:]}"))
        except ValueError:
            continue
    return days


def onset_times(body):
    times = set()
    for match in ONSET.finditer(body):
        try:
            times.add(datetime.fromisoformat(match[1].replace(" ", "T")).replace(tzinfo=UTC))
        except ValueError:
            continue
    return sorted(times)


def self_test():
    assert archive_link(b'<a href="/fixture.tar.gz">JSON (31 MB)</a>') == "https://gcn.nasa.gov/fixture.tar.gz"
    a, = positions("R.A. = 13.5 deg, DEC = −2.25 deg with an uncertainty of 3 arcmin")
    assert (a["ra_deg"], a["dec_deg"], a["nearby_quoted_radius_arcsec"]) == (13.5, -2.25, 180)
    b, = positions("(RA, Dec) = (00:54:00, -02:15:00) J2000")
    assert (b["ra_deg"], b["dec_deg"]) == (13.5, -2.25)
    assert not list(positions("the slope is 13.5, -2.25"))
    assert dates_in_name("EP240201a / GRB 240201B") == {date(2024, 2, 1)}
    assert onset_times("T0=2024-02-01T10:15:03 UTC")[0].hour == 10
    assert not onset_times("follow-up exposure at 2024-02-01T10:15:03 UTC")


def review(lead, peak, name_days, published, ep_text, positions_found, onsets):
    names_near = any(abs((d - peak.date()).days) <= 1 for d in name_days)
    publication_near = abs((published.date() - peak.date()).days) <= 1
    if not (ep_text or names_near or publication_near):
        return None
    near = []
    for coords in positions_found:
        dist = separation(lead["ra_deg"], lead["dec_deg"], coords["ra_deg"], coords["dec_deg"])
        if dist <= MAX_OFFSET:
            near.append({**coords, "offset_arcsec": dist})
    # The report publishes peak time and seed width, not both seed endpoints.
    offsets = [(t-peak).total_seconds() for t in onsets]
    close_time = any(abs(s) <= lead["trial_duration_s"] for s in offsets)
    if not near and not ((names_near or publication_near) and not positions_found):
        return None
    if near:
        relation = ("position_and_explicit_time" if close_time else
                    "position_and_named_event_day" if names_near else "position_only")
    else:
        relation = "date_neighbour_without_parsed_coordinates"
    return {"relation": relation, "event_name_dates": sorted(d.isoformat() for d in name_days),
            "near_publication_date": publication_near, "positions": near,
            "explicit_onset_utc": [t.isoformat() for t in onsets], "onset_minus_peak_s": offsets}


def main(args):
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != CAMPAIGN_REF
            or os.getenv("GITHUB_EVENT_NAME") not in {"push", "workflow_dispatch"}):
        raise RuntimeError("Alert review requires the research branch's Chandra Actions campaign")
    self_test()
    if args.self_test:
        return
    path = args.evidence / "leads.json"
    body = path.read_bytes()
    leads = json.loads(body)
    if not isinstance(leads, list):
        raise ValueError("leads.json must contain an array")
    out = args.evidence / "alerts"
    out.mkdir(parents=True, exist_ok=True)
    rows = [{"obsid": l["obsid"], "input_sha256": l["input_sha256"], "trigger_key": l["trigger_key"],
             "ra_deg": l["ra_deg"], "dec_deg": l["dec_deg"], "event_mjd": l["event_mjd"],
             "trial_duration_s": l["trial_duration_s"], "status": "unverified", "matches": []} for l in leads]
    summary = {"status": "no_leads" if not rows else "incomplete", "leads_sha256": hashlib.sha256(body).hexdigest(),
               "receipts": {}, "json_members": 0, "decoded_circulars": 0, "ep_documents": 0,
               "ep_with_parsed_coordinates": 0, "unreadable_members": [], "leads": rows,
               "interpretation": "Human-review shortlist. Examines EP/WXT text at all epochs and other circulars "
               "with publication or EP/GRB event-name date within one UTC calendar day of peak. Only labelled "
               "decimal-degree pairs or colon-separated sexagesimal pairs within six arcminutes are selected. "
               "Explicit UTC onset near peak means within one trial width. Reported coordinates and times may "
               "describe different objects. Unparsed dates, tables, wide error regions, GCN Notices, non-GCN "
               "publications and unreported transients remain outside this check. No automatic novelty verdict."}
    dump(out / "review.json", summary)
    if not rows:
        return
    peaks = [Time(l["event_mjd"], format="mjd", scale="utc").to_datetime(timezone=UTC) for l in leads]
    relevant_days = {peak.date() + timedelta(days=offset) for peak in peaks for offset in (-1, 0, 1)}
    src = out / "source"
    src.mkdir(exist_ok=True)
    try:
        page = read_input("index", INDEX, src, args.replay, summary["receipts"])
        raw = read_input("archive", archive_link(page), src, args.replay, summary["receipts"])
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r|gz") as archive:
            for member in archive:
                if not member.isfile() or not member.name.lower().endswith(".json"):
                    continue
                summary["json_members"] += 1
                try:
                    if member.size > 2_000_000:
                        raise ValueError("unexpectedly large circular")
                    data = archive.extractfile(member).read()
                    doc = json.loads(data)
                    circular_id = int(doc["circularId"])
                    published = datetime.fromtimestamp(float(doc["createdOn"])/1000, UTC)
                    title, text = doc["subject"], doc["body"]
                    if not isinstance(text, str) or not isinstance(title, str):
                        raise ValueError("unexpected text fields")
                    event = str(doc.get("eventId") or "")
                    days = dates_in_name(event + " " + title)
                    ep = bool(EP.search(event + " " + title) or EP.search(text))
                except (ValueError, UnicodeError, TypeError, KeyError, OSError, OverflowError) as exc:
                    summary["unreadable_members"].append({"member": member.name, "error": str(exc)[:160]})
                    continue
                summary["decoded_circulars"] += 1
                summary["ep_documents"] += int(ep)
                stamp = published.isoformat()
                for field, func in (("min_published_utc", min), ("max_published_utc", max)):
                    summary[field] = func(summary.get(field, stamp), stamp)
                if not (ep or days & relevant_days or published.date() in relevant_days):
                    continue
                coords = list(positions(text))
                onsets = onset_times(text)
                summary["ep_with_parsed_coordinates"] += int(ep and bool(coords))
                for row, lead, peak in zip(rows, leads, peaks, strict=True):
                    match = review(lead, peak, days, published, ep, coords, onsets)
                    if match:
                        sha = hashlib.sha256(data).hexdigest()
                        member_path = out / "members" / f"{circular_id}-{sha[:16]}.json"
                        member_path.parent.mkdir(exist_ok=True)
                        if not member_path.exists():
                            member_path.write_bytes(data)
                        row["matches"].append({**match, "id": circular_id, "subject": title, "event": event,
                                               "version": doc.get("version") or doc.get("revision"), "published_utc": stamp,
                                               "archive_member": member.name, "sha256": sha,
                                               "file": str(member_path.relative_to(out))})
        if not summary["ep_with_parsed_coordinates"]:
            raise ValueError("no EP document with parseable coordinates; archive coverage unverified")
        summary["status"] = "partial_archive" if summary["unreadable_members"] else "parsed_archive"
        for row, peak in zip(rows, peaks, strict=True):
            covered = summary["min_published_utc"] <= peak.isoformat() <= summary["max_published_utc"]
            row["peak_within_snapshot_publication_range"] = covered
            if summary["unreadable_members"] or not covered:
                row["status"] = "incomplete_coverage"
            else:
                row["status"] = ("review_shortlist" if row["matches"] else "no_match_in_parsed_subset")
    except Exception as exc:
        summary["error"] = f"{type(exc).__name__}: {exc}"
        raise
    finally:
        dump(out / "review.json", summary)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=HERE / "evidence")
    parser.add_argument("--replay", action="store_true", help="Require archived source bytes; no HTTP")
    parser.add_argument("--self-test", action="store_true", help="Check text and time formats without fetching")
    main(parser.parse_args())
