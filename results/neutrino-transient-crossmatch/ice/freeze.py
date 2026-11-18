"""Freeze primary IceCube tables, latest alert revisions and native 90% NUNIQ masks.

Usage: python freeze.py --cache .ice-cache [--out /path/to/ice]
Requires numpy and astropy. No catalogue positions are consulted.
"""
import argparse
import csv
import hashlib
import io
import json
import re
import tarfile
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from astropy.io import fits

CUTOFF = "2026-11-17T12:09:02+00:00"
CUT = datetime.fromisoformat(CUTOFF)
GCN = "https://gcn.gsfc.nasa.gov/amon_icecube_gold_bronze_events.html"
CIRC = "https://gcn.nasa.gov/circulars/archive.json.tar.gz"
ROC = "https://roc-2.icecube.wisc.edu/public/alerts/"
CDS = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJS/269/25"
DATACITE = "https://api.datacite.org/dois/"
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--cache", type=Path, required=True)
parser.add_argument("--out", type=Path, default=Path(__file__).resolve().parent)
args = parser.parse_args()
cache, out = args.cache, args.out
out.mkdir(parents=True, exist_ok=True)
cache.mkdir(parents=True, exist_ok=True)
receipts = []


def sha(data):
    return hashlib.sha256(data).hexdigest()


def get(name, url):
    path = cache / name
    if not path.exists():
        with urllib.request.urlopen(url, timeout=120) as response:
            data = response.read()
        path.write_bytes(data)
    data = path.read_bytes()
    receipts.append(dict(file=name, url=url, sha256=sha(data), bytes=len(data)))
    return data


def table(name, rows):
    with (out / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def freeze_icecat():
    raw = get("ipac", CDS)
    (out / "icecat1_table2.tsv").write_bytes(raw)
    assert sum(bool(re.match(rb"^IC\d", line)) for line in raw.splitlines()) == 275
    descriptions = []
    for filename, doi in [("icecat1dc", "10.7910/DVN/SCRUCD"), ("dc2", "10.7910/dvn/rx28yt")]:
        d = json.loads(get(filename, DATACITE + doi))["data"]["attributes"]
        descriptions.append(dict(doi=doi, title=d["titles"][0]["title"], metadata_version=d["metadataVersion"],
                                 version=d["version"], publication_year=d["publicationYear"], url=d["url"]))
    table("icecat_releases.csv", descriptions)


def notice_rows():
    raw = get("legacy.html", GCN)
    records = defaultdict(list)
    history = []
    for tr in re.findall(r"<tr\b.*?</tr>", raw.decode(), re.I | re.S):
        v = [c.strip() for c in re.findall(r"<td[^>]*>(.*?)</td>", tr, re.I | re.S)]
        if len(v) != 13:
            continue
        event = re.search(r">(\d+_\d+)<", v[0]).group(1)
        date = datetime.strptime(v[2] + " " + v[3], "%y/%m/%d %H:%M:%S.%f").replace(tzinfo=timezone.utc)
        if not (2024 <= date.year and date <= CUT):
            continue
        row = dict(run_event=event, rev=int(v[1]), utc=date.isoformat(timespec="milliseconds").replace("+00:00","Z"),
                   stream=v[4], ra_deg=float(v[5]), dec_deg=float(v[6]), error90_arcmin=float(v[7]),
                   error50_arcmin=float(v[8]), energy_notice=v[9], signalness_notice=v[10], far_per_year=v[11],
                   notice_url=f"https://gcn.gsfc.nasa.gov/notices_amon_g_b/{event}.amon")
        records[event].append(row)
        history.append(row)
    table("gcn_notice_revisions.csv", sorted(history, key=lambda r:(r["utc"],r["rev"])))
    return {event:max(revs,key=lambda r:r["rev"]) for event,revs in records.items()}


def circular_updates(latest):
    payload = get("circulars.tar.gz", CIRC)
    pubs = defaultdict(list)
    audit = []
    with tarfile.open(fileobj=io.BytesIO(payload), mode="r|gz") as archive:
        for member in archive:
            if not member.isfile():
                continue
            data = archive.extractfile(member).read()
            x = json.loads(data)
            timestamp = datetime.fromtimestamp(x["createdOn"]/1000, timezone.utc)
            if not (2024 <= timestamp.year and timestamp <= CUT):
                continue
            sub, body = x["subject"], x["body"]
            is_retraction = "icecube" in sub.lower() and "retraction" in sub.lower()
            is_position = "icecube observation of a high-energy neutrino candidate track-like event" in sub.lower()
            is_correction = int(x["circularId"]) == 35498
            if not (is_retraction or is_position or is_correction):
                continue
            date_id = re.search(r"IceCube[- ](\d{6}[A-Z])",sub,re.I)
            date_id = "IC"+date_id[1].upper() if date_id else ""
            ids = set(re.findall(r"\b(\d{6}_\d+)\b",sub+" "+body)) & set(latest)
            rows = {}
            for axis in ["RA","Dec"]:
                m = re.search(r"(?im)^\s*" + axis + r"\s*:\s*([+\-]?[\d.]+)\s*(?:deg)?\s*\(\s*\+?([\d.]+)\s*[,/]?\s*-([\d.]+)",body)
                if m: rows[axis] = [float(v) for v in m.groups()]
            audit.append(dict(circular_id=x["circularId"],utc=timestamp.isoformat(),alert_id=date_id,run_event="|".join(sorted(ids)),
                              kind="retraction" if is_retraction else "correction" if is_correction else "discovery",
                              sha256_archive_member=sha(data),url=f"https://gcn.nasa.gov/circulars/{x['circularId']}.json"))
            if not ids and date_id:
                candidates = [k for k,r in latest.items() if r["utc"][2:10].replace("-","")==date_id[2:8]]
                if len(candidates)==1: ids=set(candidates)
            for key in ids:
                pubs[key].append((timestamp,x["circularId"],date_id,is_retraction,rows))
    table("gcn_circular_audit.csv",sorted(audit,key=lambda r:r["circular_id"]))
    return pubs


def convert_map(name):
    data = get(name, ROC+name)
    with fits.open(io.BytesIO(__import__("gzip").decompress(data))) as hdus:
        h, a = hdus[1].header, hdus[1].data
        assert h["ORDERING"]=="NUNIQ" and h["TUNIT2"]=="sr-1" and "highest posterior density" in h["COMMENTS"].lower()
        uniq, density = a["UNIQ"], a["PROBDENSITY"]
        assert len(np.unique(uniq))==len(uniq)
        valid = np.isfinite(density) & (density>=0)
        order = np.floor(np.log2(uniq[valid]/4)/2).astype(int)
        area = np.pi / (3 * np.exp2(2*order))
        uniq, density = uniq[valid], density[valid]
        sort = np.argsort(-density,kind="stable")
        cum = np.cumsum((area*density)[sort])
        assert np.isfinite(cum[-1]) and cum[-1]>0
        index = np.searchsorted(cum,0.9*cum[-1])
        # Include the crossing pixel and all equal-density neighbours.
        mask = density >= density[sort[index]]
        header_area = h["CONTOUR_AREA_90"]
        measured_area = area[mask].sum()*(180/np.pi)**2
        key = f'{h["RUNID"]}_{h["EVENTID"]}'
        return dict(run_event=key,alert_id="IC"+name.split("_")[0].removeprefix("IceCube-"),utc=h["DATE-OBS"]+"Z",
                    ra_deg=float(h["RA"]),dec_deg=float(h["DEC"]),uniq90=sorted(map(int,uniq[mask])),
                    input_mass=float(cum[-1]),enclosed_mass=float(np.sum((area*density)[mask])/cum[-1]),area_deg2=float(measured_area),
                    published_area_deg2=float(header_area),pixel_deg2_at_boundary=float(area[sort[index]]*(180/np.pi)**2),
                    sha256_fits_gz=sha(data),map_url=ROC+name)


def freeze_modern():
    records = notice_rows()
    pubs = circular_updates(records)
    index = get("maps",ROC)
    names = sorted(set(re.findall(r'IceCube-\d{6}[A-Z]_skymap_probdensity_multiorder\.fits\.gz',index.decode())))
    selected = [n for n in names if n[8:14] <= CUT.strftime("%y%m%d")]
    with ThreadPoolExecutor(max_workers=5) as pool:
        maps = {m["run_event"]:m for m in pool.map(convert_map,selected)}
    rows = []
    for key,r in sorted(records.items(),key=lambda v:v[1]["utc"]):
        updates = sorted(pubs[key])
        circulars = [str(u[1]) for u in updates]
        is_retracted = any(u[3] for u in updates)
        refined = [u for u in updates if "RA" in u[4] and "Dec" in u[4]]
        provenance = refined[-1] if refined else None
        bounds = provenance[4] if provenance else {}
        alert_id = next((u[2] for u in reversed(updates) if u[2]),"")
        m = maps.get(key)
        if m: assert m["alert_id"] == alert_id and abs((datetime.fromisoformat(m["utc"])-datetime.fromisoformat(r["utc"])).total_seconds())<1
        rows.append(dict(run_event=key,alert_id=alert_id,utc=m["utc"] if m else r["utc"],stream=r["stream"],last_notice_revision=r["rev"],
                         status="retracted" if is_retracted else "active" if refined else "unreconciled",
                         ra_deg=m["ra_deg"] if m else bounds.get("RA",[r["ra_deg"]])[0],dec_deg=m["dec_deg"] if m else bounds.get("Dec",[r["dec_deg"]])[0],
                         ra_plus90=bounds.get("RA",["","",""])[1],ra_minus90=bounds.get("RA",["","",""])[2],
                         dec_plus90=bounds.get("Dec",["","",""])[1],dec_minus90=bounds.get("Dec",["","",""])[2],
                         circular_ids="|".join(circulars),refined_circular=provenance[1] if provenance else "",
                         region="native_hpd90" if m else "circular_bounds_only" if refined else "notice_radius_only"))
    table("gcn_alerts.csv",rows)
    with (out/"gcn_hpd90.jsonl").open("w") as f:
        for m in sorted(maps.values(),key=lambda m:m['utc']):
            f.write(json.dumps(m,sort_keys=True,separators=(",",":"))+"\n")


if __name__=="__main__":
    freeze_icecat()
    freeze_modern()
    table("source_receipts.csv",sorted(receipts,key=lambda r:r["file"]))
