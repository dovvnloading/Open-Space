"""Compare independent published CODE original orbits against frozen JPL osculating solutions."""
import base64
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re
import time
from urllib.parse import parse_qs, urlencode, urljoin, urlparse
from urllib.request import Request, urlopen

from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parent
OUT = ROOT/"evidence"
BASE = "https://apollo.astro.amu.edu.pl/CODE/"
OUT.mkdir(exist_ok=True)


def fetch(url):
    with urlopen(Request(url,headers={"User-Agent":"Open-Space-orbit-study/1"}),timeout=18) as resp:
        raw = resp.read()
        return BeautifulSoup(raw,"html.parser"), {
            "url":resp.url, "date":resp.headers.get("Date"), "sha256":hashlib.sha256(raw).hexdigest()}


def original(pair):
    des, internal = pair
    url = urljoin(BASE,"orbit.php?"+urlencode({"int":internal,"orb":"original"}))
    info = {"pdes":des,"solution_id":internal,"url":url}
    try:
        soup, receipt = fetch(url)
        fields = {}
        for tr in soup.find_all("tr"):
            cells = [" ".join(c.stripped_strings) for c in tr.find_all("td",recursive=False)]
            if len(cells)>=2:
                fields[cells[0].lower()] = cells[1:]
        key = next(k for k in fields if "reciprocal semi-major axis" in k)
        nominal, sigma = fields[key][:2]
        def number(s):
            return float(s.replace("±","").replace("−","-").replace(" ",""))
        model = fields["type of model of motion"][0].split(" - ")[0]
        val, err = number(nominal), number(sigma)
        return {**info,**receipt,"inverse_a_micro_au":val,"sigma_inverse_a_micro_au":err,
                "z_unbound":-val/err if err else None,"motion_model":model,
                "published_class":fields["orbit quality class"][0]}
    except Exception as exc:
        return {**info,"error":str(exc)}
    finally:
        time.sleep(0.15)


def main():
    raw = gzip.decompress((ROOT.parent/"catalog/sbdb.json.gz").read_bytes())
    cat = json.loads(raw)
    candidates = {r[cat["fields"].index("pdes")]:r[0] for r in cat["data"]}
    soup, index = fetch(BASE)
    links = {}
    for a in soup.find_all("a",href=True):
        matches = re.search(r"([0-9]{4})\s+([A-Z]{1,2}[0-9]+(?:-[A-Z])?)"," ".join(a.stripped_strings))
        sol = parse_qs(urlparse(a["href"]).query).get("int",[])
        if matches and sol:
            des = " ".join(matches.groups())
            if des in candidates:
                if des in links and links[des]!=sol[0]:
                    raise RuntimeError("ambiguous preferred solution for "+des)
                links[des]=sol[0]
    if not links:
        raise RuntimeError("no preferred solution links matched the catalogue")
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(original, sorted(links.items())))
    for r in rows:
        r["spkid"]=candidates[r["pdes"]]
    dest = OUT/"code_original.csv"
    with dest.open("w", newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for r in rows for k in r)),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    valid=[r for r in rows if "inverse_a_micro_au" in r]
    result={"index":index,"jpl_sha256":hashlib.sha256(raw).hexdigest(),"matched":len(rows),
            "parsed":len(valid),"nominal_bound":sum(r["inverse_a_micro_au"]>0 for r in valid),
            "formal_unbound_z_ge_3":[r["pdes"] for r in valid if (r["z_unbound"] or 0)>=3],
            "formal_unbound_z_ge_5":[r["pdes"] for r in valid if (r["z_unbound"] or 0)>=5],
            "models":dict(collections.Counter(r["motion_model"] for r in valid))}
    summary=OUT/"code_summary.json"
    summary.write_text(json.dumps(result,indent=2)+"\n")
    print("CODE",result,flush=True)
    for f in (dest,summary):
        print("RESULT",f.name,base64.b64encode(f.read_bytes()).decode())
    if not valid:
        raise RuntimeError("no usable original orbit rows; inspect the CSV errors")

if __name__ == "__main__":
    main()
