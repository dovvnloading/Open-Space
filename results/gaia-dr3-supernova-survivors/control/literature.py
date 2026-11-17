"""Extract compact control data from the hashed primary source using Poppler."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import tempfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
URL = "https://arxiv.org/pdf/2606.11293"
DIGEST = "245ca54d647c0844b2f37efe63af532a457a3f0585b7eb3adee08c3828259920"


def write_table(name, fields, rows):
    with (HERE / name).open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def extract(pdf):
    data = pdf.read_bytes()
    if hashlib.sha256(data).hexdigest() != DIGEST:
        raise ValueError("The El-Badry et al. PDF differs from the analysed version")
    raw = subprocess.check_output(["pdftotext", "-layout", str(pdf), "-"], text=True)
    raw = raw.replace("−", "-")
    controls = list(csv.DictReader((HERE / "controls.csv").open()))
    control_ids = {r["source_id"] for r in controls}
    # Table 2: reported central value +/- RV error and median +hi -lo distance.
    # In -layout extraction the next line starts with the distance's lower error.
    kin_pattern = (r"^\s*(?:D6-\d|J\d{4}|LP 40-365|US 708)\s+(\d{17,19})\s+"
                   r"[\d.]+\s+-?[\d.]+\s+(-?\d+)\s*±\s*(\d+)\s+"
                   r"([\d.]+)\+([\d.]+)\s*\n *-([\d.]+)")
    kin = []
    for ident, rv, err, med, plus, minus in re.findall(kin_pattern, raw, re.M):
        kin.append(dict(source_id=ident, rv_kms=rv, rv_err_kms=err, distance_med_kpc=med,
                        distance_lo_kpc=f"{float(med)-float(minus):.2f}",
                        distance_hi_kpc=f"{float(med)+float(plus):.2f}"))
    assert len(kin) == 20 and {r["source_id"] for r in kin} == control_ids
    write_table("published_kinematics.csv", kin[0].keys(), kin)

    # Table 4: the source ID and spectroscopic/imaging verdict are unambiguous.
    appendix = raw.split("TABLE 4: All Gaia sources", 1)[1]
    known_verdicts = r"cool D6|hot D6|LP 40-365|blended|DA WD|DB WD|DC WD|DA|MS|sdO/B|WB"
    targets = []
    for line in appendix.splitlines():
        ident = re.match(r"^\s*(\d{17,19})\s+", line)
        if not ident:
            continue
        matched = re.search(r"\s(" + known_verdicts + r")(?:\s+|$)", line)
        if not matched:
            raise ValueError(f"Could not read verdict for {ident.group(1)}")
        targets.append(dict(source_id=ident.group(1), verdict=matched.group(1)))
    assert len({r["source_id"] for r in targets}) == len(targets) == 92
    assert len({r["source_id"] for r in targets} & control_ids) == 13
    write_table("published_targets.csv", targets[0].keys(), targets)
    print(len(kin), "controls and", len(targets), "literature targets")
    receipt = {"url": URL, "title": "A systematic survey for hypervelocity runaways from thermonuclear supernovae",
               "doi": "10.33232/001c.164326", "pdf_version": "June 30, 2026 (title page)",
               "bytes": len(data), "sha256": DIGEST, "tables": [2,4],
               "extraction": "Poppler pdftotext -layout, exact Gaia DR3 identifier join"}
    (HERE / "literature_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdf", type=Path)
    args = parser.parse_args()
    if args.pdf:
        extract(args.pdf)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            target = Path(tmp) / "article.pdf"
            urllib.request.urlretrieve(URL, target)
            extract(target)
