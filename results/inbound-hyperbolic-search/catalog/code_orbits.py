"""Cross-match preferred published CODE solutions as a separate check of the JPL screen."""
import argparse
import csv
import gzip
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import time
from urllib.parse import urljoin, urlsplit, parse_qs, urlencode
from urllib.request import urlopen, Request

BASE = "https://apollo.astro.amu.edu.pl/CODE/"
UA = "Open-Space-orbit-validation/1.0 (+https://github.com/dovvnloading/Open-Space)"


class Tables(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.rows, self.links, self.row, self.cell, self.link = [], [], [], None, None

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == "tr":
            self.row = []
        if tag in ("th", "td"):
            self.cell = ""
        if tag == "a" and "href" in attrs:
            self.link = [attrs["href"], ""]

    def handle_data(self, data):
        if self.cell is not None:
            self.cell += data
        if self.link is not None:
            self.link[1] += data

    def handle_endtag(self, tag):
        if tag == "a" and self.link:
            self.links.append(self.link)
            self.link = None
        elif tag in ("th", "td") and self.cell is not None:
            self.row.append(" ".join(self.cell.split()))
            self.cell = None
        elif tag == "tr" and self.row:
            self.rows.append(self.row)
            self.row = []


def download(url):
    req = Request(url, headers={"User-Agent": UA})
    with urlopen(req, timeout=30) as response:
        raw = response.read()
        stamp = response.headers.get("Date")
    doc = Tables()
    decoded = raw.decode("utf-8")
    doc.feed(decoded)
    doc.barycentric = "barycentric ecliptic" in decoded.lower()
    return doc, {"url":url,"bytes":len(raw),"sha256":hashlib.sha256(raw).hexdigest(),"http_date":stamp}


def des(name):
    match = re.search(r"\b([CPAXD]/\d{4} [A-Z]{1,2}\d*(?:-[A-Z]+)?)\b", " ".join(name.split()))
    return match.group(1) if match else ""


def main(out):
    jpl = json.loads(gzip.decompress((out / "sbdb.json.gz").read_bytes()))
    src = {}
    for values in jpl["data"]:
        r = dict(zip(jpl["fields"],values))
        ident = des(r["full_name"])
        if ident:
            src[ident] = r["spkid"]
    index, index_receipt = download(BASE)
    models = {}
    for link, name in index.links:
        params = parse_qs(urlsplit(link).query)
        ident = des(name)
        if ident in src and params.get("orb")==["osculating"] and "int" in params:
            models[ident] = params["int"][0]
    if not models:
        raise ValueError(f"No preferred models matched: sample {index.links[:15]}")
    print("CODE matched",len(models),"preferred orbit IDs",flush=True)
    outputs = []
    receipts = [index_receipt]
    for i, (ident, model) in enumerate(sorted(models.items())):
        url = urljoin(BASE,"orbit.php") + "?" + urlencode({"int":model,"orb":"original"})
        try:
            doc, rec = download(url)
            receipts.append(rec)
            flat = [x for row in doc.rows for x in row]
            if not doc.barycentric and not any("barycentric" in x.lower() for x in flat):
                raise ValueError(f"Original barycentric header missing: {doc.rows[-10:]}")
            energy = next(row for row in doc.rows if row and "reciprocal semi-major axis" in row[0].lower())
            def value(s):
                return float(s.replace(",","").replace("±","").replace("−","-").strip())
            if len(energy) < 3:
                raise ValueError(f"Unexpected original energy columns: {energy}")
            inva, sig = value(energy[1]), value(energy[2])
            if sig <= 0:
                raise ValueError(f"Invalid published sigma: {sig}")
            outputs.append({"designation":ident,"spkid":src[ident],"model":model,
                            "published_original_1_over_a_micro_au_inv":inva,"sigma":sig,
                            "url":url})
        except Exception as exc:
            print("CODE issue",ident,str(exc)[:450],flush=True)
            if i==0:
                print("CODE fields",doc.rows[:20] if 'doc' in locals() else "",flush=True)
            receipts.append({"url":url,"designation":ident,"error":str(exc)})
            if i>=4 and not outputs:
                break
        time.sleep(0.15)
    if outputs:
        with (out/"code-original.csv").open("w",newline="") as f:
            w = csv.DictWriter(f,fieldnames=list(outputs[0]),lineterminator="\n")
            w.writeheader()
            w.writerows(outputs)
    (out/"code-receipts.json").write_text(json.dumps(receipts,indent=2)+"\n")
    summary={"matched":len(models),"fetched_with_uncertainty":len(outputs),
             "negative_1_over_a":sum(r["published_original_1_over_a_micro_au_inv"]<0 for r in outputs),
             "negative_beyond_five_sigma":[r for r in outputs if r["published_original_1_over_a_micro_au_inv"]+5*r["sigma"]<0]}
    if models and not outputs:
        raise ValueError("No usable original orbit rows")
    (out/"code-summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("CODE summary",json.dumps(summary),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("out",type=Path)
    main(p.parse_args().out)
