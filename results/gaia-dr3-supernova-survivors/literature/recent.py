"""Date-bounded arXiv abstract inventory and exact identifier checks in primary studies."""
import argparse
import csv
import hashlib
import json
import re
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen

HERE=Path(__file__).resolve().parent
QUERY=('submittedDate:[202601010000 TO 202611172359] AND (abs:hypervelocity OR abs:"LP 40-365" '
       'OR abs:"partially burned" OR abs:"US 708" OR (abs:"supernova" AND abs:"runaway" AND abs:"white dwarf"))')
BASE="https://export.arxiv.org/api/query?"
PAPERS=["1804.11163","1902.05061","2306.03914","2311.13388","2408.08397","2506.08081",
        "2602.23900","2606.11293","2209.09915","1910.13053"]


def main(cache):
    cache.mkdir(exist_ok=True,parents=True)
    def fetch(name,url):
        f=cache/name
        if not f.exists():
            f.write_bytes(urlopen(url,timeout=100).read())
        data=f.read_bytes()
        return f,dict(url=url,sha256=hashlib.sha256(data).hexdigest(),bytes=len(data))
    url=BASE+urlencode(dict(search_query=QUERY,sortBy="submittedDate",sortOrder="descending",start=0,max_results=180))
    atom,receipt=fetch("arxiv_recent.atom",url)
    ns={"a":"http://www.w3.org/2005/Atom","op":"http://a9.com/-/spec/opensearch/1.1/"}
    root=ET.parse(atom).getroot()
    all_entries=root.findall("a:entry",ns)
    if len(all_entries)!=int(root.findtext("op:totalResults",namespaces=ns)):
        raise ValueError("abstract results require pagination")
    with (HERE/"recent_arxiv.csv").open("w",newline="") as out:
        w=csv.writer(out);w.writerow(["arxiv","title"])
        for e in all_entries:
            w.writerow([e.findtext("a:id",namespaces=ns).split("/")[-1],
                        " ".join(e.findtext("a:title",namespaces=ns).split())])
    with (HERE/"matches.csv").open() as f:
        ids=set(r["source_id"] for r in csv.DictReader(f))
    for name in ["metals/targets.csv","hot/tail.csv","sdss/tail.csv","spectra/lamost_high_velocity.csv"]:
        with (HERE.parent/name).open() as stream:
            for row in csv.DictReader(stream):
                sid = row.get("source_id",row.get("GaiaEDR3",""))
                if sid.isdigit():
                    ids.add(sid)
    hot = json.loads((HERE.parent/"hot/desi_summary.json").read_text())
    ids.update(hot["raw_ge_700_ids"])
    papers=[]
    for name in PAPERS:
        pdf,info=fetch("paper_"+name+".pdf","https://arxiv.org/pdf/"+name)
        text=subprocess.check_output(["pdftotext","-layout",str(pdf),"-"],text=True)
        in_text=set(re.findall(r"(?<!\d)\d{17,19}(?!\d)",text))
        papers.append(dict(**info,arxiv=name,matched_ids=sorted(ids&in_text)))
    (HERE/"recent_receipt.json").write_text(json.dumps(dict(abstract_search=receipt,
        query=QUERY,returned=len(all_entries),checked_ids=sorted(ids),primary_pdfs=papers),indent=2)+"\n")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("cache",type=Path)
    main(p.parse_args().cache)
