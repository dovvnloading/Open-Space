"""Freeze observational inputs for independent orbit fits.

The MPC comet bulk file is filtered by packed permanent or provisional designation.
The optional MPEC source is a point-in-time fallback, flagged per object.
"""
from datetime import datetime, timezone
from pathlib import Path
import gzip
import hashlib
import html
import json
import re
import urllib.error
import urllib.parse
import urllib.request

DATA = Path(__file__).parent / "data"
CACHE = Path(__file__).parent / "cache"
DATA.mkdir(exist_ok=True)
CACHE.mkdir(exist_ok=True)
RECEIPTS = DATA / "acquisition.json"
NAMES = {
    "1I": (["0001I", "K17U010", "K17UI1I"], "2017-U181"),
    "2I": (["0002I", "K19Q040"], "2019-R106"),
    "3I": (["0003I", "K25N010"], "2025-N12"),
    "2024-L5": (["CK24L050"], "2024-O15"),
    "1980-E1": (["CJ80E010"], None),
    "1997-P2": (["CJ97P020"], None),
    "2008-J4": (["CK08J040"], None),
}
def digest(b):
    return hashlib.sha256(b).hexdigest()

def retrieve(url, name, receipts, body=None):
    p = CACHE / name
    if p.exists():
        return p.read_bytes()
    headers = {"User-Agent": "Open-Space orbit research/1.0"}
    if body:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=body, method="GET", headers=headers)
    rec = {"requested_url": url, "utc": datetime.now(timezone.utc).isoformat()}
    if body:
        rec["get_json_body"] = json.loads(body)
    receipts.append(rec)
    try:
        with urllib.request.urlopen(req, timeout=75) as r:
            b = r.read()
            rec.update(status=r.status, final_url=r.url, bytes=len(b), sha256=digest(b),
                       last_modified=r.headers.get("Last-Modified"))
        p.write_bytes(b)
        return b
    except (OSError, ValueError) as exc:
        rec["error"] = str(exc)
        raise

def select(source, keys):
    return [s[:80] for s in source if len(s)>=80 and any(k in s[:12] for k in keys)
            and re.match(r"\d{4} \d\d \d\d",s[15:25])]

def acquire():
    receipts = []
    selected = {}
    urls = [
        "https://www.minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz",
        "https://minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz",
    ]
    catalog = []
    for k,url in enumerate(urls):
        try:
            catalog = gzip.decompress(retrieve(url,f"cmt-{k}.gz",receipts)).decode("ascii","replace").splitlines()
            if not catalog:
                raise ValueError("empty MPC bulk file")
            break
        except (OSError, ValueError, EOFError):
            continue
    for obj,(keys,circ) in NAMES.items():
        lines = select(catalog,keys)
        src = "MPC bulk CmtObs"
        if not lines:
            desig = obj if obj.endswith("I") else "C/"+obj.replace("-"," ")
            body = json.dumps({"desigs":[desig],"output_format":["OBS80"]}).encode()
            try:
                response = json.loads(retrieve("https://data.minorplanetcenter.net/api/get-obs",
                                      obj+".json",receipts,body))
                field = response[0]["OBS80"]
                lines = select(field.splitlines() if isinstance(field,str) else field,keys)
                src = "MPC Observations API"
            except (ValueError, OSError, KeyError, TypeError, IndexError):
                pass
        if not lines and circ:
            year,letter = circ.split("-")
            seq = int(letter[1:]); lead = str(seq//10) if seq<100 else chr(87+seq//10).upper()
            tail = lead+str(seq%10)
            url = f"https://minorplanetcenter.net/mpec/K{year[2:]}/K{year[2:]}{letter[0]}{tail}.html"
            try:
                web = retrieve(url,obj+".html",receipts).decode("utf-8","replace")
                lines = select(html.unescape(re.sub(r"<[^>]*>","",web)).splitlines(),keys)
                src = "MPC discovery circular "+circ
            except OSError:
                pass
        if lines:
            p = DATA/(obj+".obs")
            p.write_text("\n".join(lines)+"\n")
            selected[obj] = {"source":src,"records":len(lines),"sha256":digest(p.read_bytes()),
                             "first_date":min(s[15:32].strip() for s in lines),
                             "last_date":max(s[15:32].strip() for s in lines)}
        else:
            selected[obj] = {"records":0,"status":"no authenticated astrometric rows retrieved"}
        print(obj, selected[obj], flush=True)
    for file,urls in {
        "stations.txt":["https://minorplanetcenter.net/iau/lists/ObsCodes.html",
                        "https://www.projectpluto.com/obscodes.htm"],
        "de440s.bsp":["https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de440s.bsp",
                      "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp"],
    }.items():
        for k,url in enumerate(urls):
            try:
                b=retrieve(url,f"{k}-{file}",receipts)
                if file.endswith("bsp") and digest(b)!="c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2":
                    raise ValueError("DE440s checksum differs from prior repository verification")
                if file=="stations.txt" and b"Greenwich" not in b:
                    raise ValueError("missing observatory table")
                (CACHE/file).write_bytes(b)
                break
            except (ValueError,OSError) as exc:
                print(file, str(exc),flush=True)
    output = {"inputs":receipts,"objects":selected}
    RECEIPTS.write_text(json.dumps(output,indent=2)+"\n")
    print("ACQUIRE_RECEIPT",json.dumps(output),flush=True)
    if not any(s.get("records") for s in selected.values()):
        raise RuntimeError("No observation records retrieved")

if __name__=="__main__":
    acquire()
