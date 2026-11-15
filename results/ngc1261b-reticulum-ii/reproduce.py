"""Rebuild the Reticulum II cross-identifications (stdlib; Poppler for the PDF)."""
import argparse
import csv
import hashlib
import math
import re
import subprocess
from collections import Counter
from pathlib import Path
from urllib.request import urlopen

SOURCES = {
    "ibata.txt": (
        "https://content.cld.iop.org/journals/0004-637X/967/2/89/revision1/apjad382dt1_mrt.txt",
        "76742daaf44d69e1afe71b56f49ec93c88631a6ae10261cc4adcb6f8485644ba",
    ),
    "ji.txt": (
        "https://content.cld.iop.org/journals/1538-3881/165/3/100/revision1/ajacad84t1_mrt.txt",
        "4a0c0576b3e93c38beb1334e5a0df690d7509b4a6611c479b965ce0e4068a4d6",
    ),
    "simon.pdf": (
        "https://arxiv.org/pdf/1504.02889",
        "2330033fe384198c6d3e52d586836c24531b23fb25cbab8e7e73cc62a66d9e8d",
    ),
}


def get_sources(cache):
    cache.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in SOURCES.items():
        path = cache / name
        if path.exists():
            data = path.read_bytes()
        else:
            with urlopen(url, timeout=90) as response:
                data = response.read()
        if hashlib.sha256(data).hexdigest() != expected:
            raise ValueError(f"Input changed: {url}. Inspect before updating its digest.")
        if not path.exists():
            path.write_bytes(data)


def separation(ra1, dec1, ra2, dec2):
    # Haversine on the unit sphere. Native coordinates as tabulated in each archive.
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    h = math.sin((d1-d2)/2)**2 + math.cos(d1)*math.cos(d2)*math.sin((r1-r2)/2)**2
    return 7200 * math.degrees(math.asin(math.sqrt(min(1.0, h))))


def write_csv(path, rows):
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(cache, output):
    get_sources(cache)
    # Many Gaia IDs exceed the exact range of a binary64 float. Keep them as strings.
    ibata = [line.split() for line in (cache/"ibata.txt").read_text().splitlines()
             if line[:19].strip().isdigit() and len(line.split()) == 12]
    ji = [line for line in (cache/"ji.txt").read_text().splitlines()
          if line[:3].strip().isdigit() and line[15:16] == "-"]
    by_id = {line[93:112].strip(): line for line in ji}
    if len(ibata) != 24540 or len(ji) != 129 or len(by_id) != len(ji):
        raise ValueError("Unexpected primary catalog schema or population")
    pdf_text = subprocess.run(["pdftotext", "-layout", str(cache/"simon.pdf"), "-"],
                              check=True, capture_output=True).stdout.decode()
    table_row = re.compile(r"^\s*DES (J\d+\.\d+[-−]\d+\.\d+)\s+(\d+\.\d+)"
                           r"\s+(-\d+\.\d+)\s+.*?\s+(NM\??|M\??)\s*$")
    simon = []
    for line in pdf_text.splitlines():
        match = table_row.fullmatch(line)
        if match:
            name, ra, dec, flag = match.groups()
            simon.append((name.replace("−", "-"), float(ra), float(dec), flag))
    if len(simon) != 59:
        raise ValueError("Unexpected number of Simon table rows")

    all_rows = Counter()
    measured = Counter()
    ji_counts = Counter()
    matched_rv = Counter()
    matches = []
    for id_, ra, dec, plx, pmra, pmdec, g, color, rv, error, source, stream_id in ibata:
        stream = int(stream_id)
        has_velocity = int(source) > 0 and float(error) < 300
        all_rows[stream] += 1
        measured[stream] += has_velocity
        if id_ not in by_id:
            continue
        line = by_id[id_]
        ji_counts[stream] += 1
        matched_rv[stream] += has_velocity
        nearby = [(separation(float(ra), float(dec), r, d), name, flag)
                  for name, r, d, flag in simon]
        nearby = sorted(item for item in nearby if item[0] <= 1)
        if len(nearby) > 1:
            raise ValueError(f"Ambiguous Simon match for {id_}")
        if stream != 14:
            continue
        sep, name, flag = nearby[0] if nearby else (float("nan"), "", "")
        matches.append(dict(gaia_source_id=id_,ra_deg=ra,dec_deg=dec,
                            ibata_v_kms=rv,ibata_err_kms=error,ibata_source=source,
                            ji_id=int(line[:3]),ji_member=int(line[71]),ji_candidate=int(line[73]),
                            ji_binary=int(line[75]),ji_v_kms=line[45:51].strip(),
                            ji_err_kms=line[52:57].strip(),ji_elliptical_rh=line[40:44].strip(),
                            simon_name=name,simon_member=flag,simon_sep_arcsec=f"{sep:.3f}"))
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output/"matches.csv", matches)
    summary = [dict(stream_id=s, rows=all_rows[s], published_rvs=measured[s],
                    ji_id_matches=ji_counts[s], rvs_among_matches=matched_rv[s],
                    remaining_rvs=measured[s]-matched_rv[s]) for s in sorted(all_rows)]
    write_csv(output/"coverage.csv", summary)
    print(f"{len(ibata)} atlas rows, {len(all_rows)} streams, {len(ji)} Ji targets, "
          f"{len(simon)} Simon targets. {len(matches)} matches in stream 14.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.cache, args.output)
