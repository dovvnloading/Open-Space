"""Retrieve LSXPS first-detection epochs and observation IDs. Run in the receipt cache."""
import csv,gzip,concurrent.futures,urllib.request,re,json,hashlib
from lxml import html
rows=[x for x in csv.DictReader(gzip.open("swift11.gz","rt")) if x["Classification"] in ("1","3") and "2024-01"<=x["DiscoveryDate"]<"2026-11-17T12:09:02"]
def get(r):
    url="https://www.swift.ac.uk/LSXPS/transients/"+r["TransientID"]
    try:
        with urllib.request.urlopen(url,timeout=30) as resp:
            raw=resp.read()
        page=html.fromstring(raw)
        dd={d.text_content().strip():d.getnext() for d in page.xpath("//dt") if d.getnext() is not None}
        first=next(v for k,v in dd.items() if k.startswith("First detected"))
        obs=next((v for k,v in dd.items() if 'Discovery obs' in k or 'Detection dataset' in k),first)
        obsid=re.findall(r"/(?:Fields|datasets?)/(\d+)",html.tostring(obs,encoding="unicode"),re.I)
        if not obsid:obsid=re.findall(r"/LSXPS/Fields/(\d+)",raw.decode())
        epoch=re.search(r"\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}",first.text_content())[0]

        return [r["TransientID"],epoch,obsid[0] if obsid else "",hashlib.sha256(raw).hexdigest(),url,first.text_content()[:90]]
    except Exception as e:return [r["TransientID"],"","",str(e),url,""]
with concurrent.futures.ThreadPoolExecutor(5) as p:
    out=list(p.map(get,rows))
for i,item in enumerate(out):
    if not item[1]:out[i]=get(rows[i])
if not all(i[1] for i in out):raise RuntimeError("Missing public first-detection page")
json.dump(out,open("swift_details.json","w"))
print('detail pages with date',sum(bool(x[1]) for x in out),len(out),flush=True)
