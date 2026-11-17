"""Reconcile Vincken et al. Table 5, Zenodo 17425610 and the two pinned Q1 catalogues."""
import argparse
import csv
import hashlib
import html
import json
import re
from pathlib import Path
from urllib.request import urlopen

from screen import output_csv, separation

ARCHIVE = "https://zenodo.org/records/17425610/files/q1_missed_targets_AstroVink_catalog.csv?download=1"
PAPER = "https://arxiv.org/html/2604.21977v1"
MD5 = "50f33adc6d2a0ccec1608fa786f03f62"


def fetch(directory, name, url):
    path = directory / name
    if not path.exists():
        with urlopen(url, timeout=120) as src:
            path.write_bytes(src.read())
    return path


def plain(s):
    s = re.sub(r'<math\b[^>]*alttext="([^"]+)"[^>]*>.*?</math>', r'\1', s, flags=re.S)
    return html.unescape(re.sub(r"<[^>]+>", "", s)).strip().replace("\u2009", "").replace("−", "-")


def run(directory, output):
    archive = fetch(directory, "astrovink.csv", ARCHIVE)
    paper = fetch(directory, "astrovink_arxiv.html", PAPER)
    if hashlib.md5(archive.read_bytes()).hexdigest() != MD5:
        raise ValueError("Unexpected AstroVink catalogue")
    with archive.open() as stream:
        released = {r["name"]: r for r in csv.DictReader(stream)}
    slde = [(r["id_str"], float(r["right_ascension"]), float(r["declination"]))
            for r in csv.DictReader((output/"slde_screen.csv").open())]
    agile = [(r["name"], float(r["ra_deg"]), float(r["dec_deg"]))
             for r in csv.DictReader((output/"agilelens_screen.csv").open())]
    pattern = r'<tr id="A5\.T5\.2\.\d+"[^>]*>(.*?)</tr>'
    rows = []
    for tr in re.findall(pattern, paper.read_text(), re.S):
        cells = [plain(td) for td in re.findall(r"<td\b[^>]*>(.*?)</td>", tr, re.S)]
        if len(cells) != 6 or not cells[0].startswith("EUCLJ"):
            continue
        name, ra, dec, score, grade, discovery = cells
        entry = released.get(name, {})
        ra, dec = float(ra), float(dec)
        if entry and (entry["grade"] != grade or separation(ra, dec, float(entry["ra"]), float(entry["dec"])) > 0.08):
            raise ValueError("Table and Zenodo coordinates or grade disagree on " + name)
        row = {"name": name, "ra_deg": ra, "dec_deg": dec, "grade": grade, "vi_score": score,
               "paper_discovery": discovery, "zenodo_id": entry.get("id_str", "")}
        near_archive = min((separation(ra, dec, float(r["ra"]), float(r["dec"])), r["id_str"], r["name"])
                           for r in released.values() if r["grade"] in {"A","B"})
        row.update({"archive_nearest_ab_arcsec": near_archive[0], "archive_nearest_ab_id": near_archive[1],
                    "archive_nearest_ab_name": near_archive[2]})
        for label, others in [("slde", slde), ("agile", agile)]:
            dist, other = min((separation(ra, dec, ra2, de2), key) for key, ra2, de2 in others)
            row["nearest_" + label] = other
            row[label + "_arcsec"] = dist
        rows.append(row)
    if len(rows) != 36 or len({r["name"] for r in rows}) != len(rows):
        raise ValueError("Unexpected Table 5 row count or duplicate names")
    output_csv(output, "astrovink_screen.csv", rows)
    summary = {"inputs": [{"url": u, "sha256": hashlib.sha256(p.read_bytes()).hexdigest()}
                          for u, p in [(ARCHIVE, archive), (PAPER, paper)]],
               "paper_rows": len(rows), "release_ab_rows": sum(r["grade"] in {"A","B"} for r in released.values()),
               "paper_rows_not_in_release": [r["name"] for r in rows if not r["zenodo_id"]],
               "within_2arcsec": {k: sum(r[k + "_arcsec"] < 2 for r in rows) for k in ["slde","agile"]},
               "within_10arcsec": {k: sum(r[k + "_arcsec"] < 10 for r in rows) for k in ["slde","agile"]}}
    (output/"astrovink_summary.json").write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_dir", type=Path)
    p.add_argument("--output", type=Path, default=Path(__file__).resolve().parent)
    args = p.parse_args()
    run(args.input_dir, args.output)
