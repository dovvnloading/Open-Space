"""Refit MPC comet observations by designation using a separately compiled Find_Orb."""
import argparse
import base64
import collections
import concurrent.futures
import csv
import gzip
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import re
import subprocess
import tempfile
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "evidence"
OUT.mkdir(exist_ok=True)
SOURCE = "https://www.minorplanetcenter.net/iau/ECS/MPCAT-OBS/CmtObs.txt.gz"
FITTER = Path.home() / "bin" / "fo"
ALPHABET = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"
DAY_PER_YEAR = 365.25


def cat():
    spec = importlib.util.spec_from_file_location("audit", ROOT.parent / "catalog/audit.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    frozen = ROOT.parent / "catalog/sbdb.json.gz"
    if frozen.is_file():
        raw = gzip.decompress(frozen.read_bytes())
    else:
        with urlopen(module.URL, timeout=90) as response:
            raw = response.read()
    data = json.loads(raw)
    assert len(data["data"]) == int(data["count"])
    return [dict(zip(data["fields"], r)) for r in data["data"]], module.URL, hashlib.sha256(raw).hexdigest()


def aliases(row):
    """Match legacy MPC designation bytes, preserving component identifiers."""
    raw = row["pdes"].strip()
    out = set()
    for s in (raw, row["full_name"].strip().split(" (")[0]):
        if "/" in s:
            a, b = s.split("/", 1)
            if len(a) > 1 and a[:-1].isdigit() and a[-1] in "IPD":
                out.add(f"{int(a[:-1]):04}{a[-1]}")
            s = b
        m = re.fullmatch(r"(\d{4}) ([A-Z])([A-Z]?)(\d*)(?:-([A-Z]))?", s)
        if m:
            yr, half, second, serial, fragment = m.groups()
            n = int(serial or 0)
            if n // 10 < len(ALPHABET):
                out.add(ALPHABET[int(yr)//100] + yr[2:] + half + ALPHABET[n//10]
                        + str(n%10) + (second or (fragment or "0").lower()))
    for iso, alias in {"1I":"K17U010", "2I":"K19Q040", "3I":"K25N010"}.items():
        if raw == iso or row["full_name"].strip().startswith(iso+"/"):
            out.update([alias, f"{int(iso[:-1]):04}I"])
    return out


def distance_and_direction(el):
    # Solve for radial distance at the reported epoch from its osculating conic.
    e, q, dt = el["e"], el["q"], el["epoch"] - el["Tp"]
    mu = 0.0002959122082855911 * 1.001342
    if e == 1:
        B = dt * math.sqrt(mu / (2*q**3))
        D = 2*math.sinh(math.asinh(1.5*B)/3)
        return q*(1+D*D), D < 0
    a = q / abs(1-e)
    M = dt * math.sqrt(mu/a**3)
    if e < 1:
        M = math.remainder(M, 2*math.pi)
    m = abs(M)
    lo, hi = 0., (math.pi if e < 1 else math.asinh(m/e)+3.)
    for _ in range(65):
        u = (lo+hi)/2
        value = u-e*math.sin(u) if e < 1 else e*math.sinh(u)-u
        if value < m:
            lo = u
        else:
            hi = u
    u = (lo+hi)/2
    return a*((1-e*math.cos(u)) if e<1 else (e*math.cosh(u)-1)), M < 0


def fit(part):
    row, lines = part
    raw = b"".join(lines)
    meta = {"spkid": row["spkid"], "name": row["full_name"].strip(),
            "input_lines": len(lines), "astrometry_sha256": hashlib.sha256(raw).hexdigest() if lines else ""}
    if not lines:
        return {**meta, "status": "no unique matching records in bulk file"}
    if len(lines) < 3:
        return {**meta, "status": "under 3 lines"}
    tp = float(row["tp"])
    # The kernel covers 1550–2650. Many ancient catalogue solutions predate it.
    if not 2287924 <= tp <= 2688977:  # range starts 1551, ends before 2650
        return {**meta, "status": "perihelion outside kernel range"}
    target = max(2287924., tp-300*DAY_PER_YEAR)
    if target > tp-3650:
        return {**meta, "status":"insufficient pre-perihelion ephemeris"}
    with tempfile.TemporaryDirectory() as d:
        stage = Path(d)
        inp = stage/"observations.mpc"
        unpenalized=os.environ.get("FO_UNPENALIZED","0")=="1"
        staged=(b"#interstellar\n" if unpenalized else b"")+raw
        inp.write_bytes(staged)
        meta.update(unpenalized_initial_hyperbola=unpenalized,staged_sha256=hashlib.sha256(staged).hexdigest())
        # Unique output directory avoids reuse of another object's results.
        cmd = [str(FITTER), str(inp), "-i", "-h-1", "-tEjd"+str(target),
               "-O"+d, "-r18,20", "PERTURBERS=7fe", "JSON_SHORT_ELEMENTS="+d+"/solution.json",
               "JSON_ELEMENTS_NAME="+d+"/long.json", "JSON_COMBINED_NAME="+d+"/combined.json"]
        try:
            run = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=27)
            dest = stage/"solution.json"
            if run.returncode != 0 or not dest.is_file():
                return {**meta, "status":"fit timeout" if run.returncode in (-9, -15, -24) else "fitter did not return elements",
                        "diagnostic":run.stdout[-350:].decode(errors="replace")}
            answer = json.loads(dest.read_text())
            obj = next(iter(answer["objects"].values()))
            el, obs = obj["elements"], obj["observations"]
            r, inbound = distance_and_direction(el)
            center = el["central body"]
            return {**meta, "status":"fit returned", "epoch":el["epoch"], "requested_epoch":target,
                    "center":center, "r_au":r, "inbound":inbound, "e":el["e"], "sigma_e":el.get("e sigma"),
                    "q_au":el["q"], "inv_a_au":(1-el["e"])/el["q"], "rms_arcsec":el["rms_residual"],
                    "used":obs["used"],"total":obs["count"],"first_used_jd":obs["earliest_used"],"last_used_jd":obs["latest_used"]}
        except subprocess.TimeoutExpired:
            return {**meta,"status":"fit timeout"}
        except (KeyError, OSError, ValueError, ArithmeticError) as exc:
            return {**meta,"status":"unusable fit output","diagnostic":str(exc)}


def download_source(directory):
    directory.mkdir(exist_ok=True,parents=True)
    request=Request(SOURCE, headers={"User-Agent":"Open-Space-public-astrometry-research/1"})
    path=directory/"CmtObs.txt.gz"
    with urlopen(request,timeout=150) as response,path.open("wb") as out:
        digest=hashlib.sha256()
        for chunk in iter(lambda:response.read(1<<20),b""):
            out.write(chunk)
            digest.update(chunk)
        receipt={"mpc_url":response.url,"http_date":response.headers.get("Date"),
                 "last_modified":response.headers.get("Last-Modified"),"bytes":out.tell(),
                 "sha256_compressed":digest.hexdigest()}
    (directory/"source.json").write_text(json.dumps(receipt,indent=2)+"\n")


def main(source_dir,shard,parts):
    rows, caturl, cathash = cat()
    keys = collections.defaultdict(set)
    for row in rows:
        for k in aliases(row):
            keys[k].add(row["spkid"])
    grouped = collections.defaultdict(list)
    source=source_dir/"CmtObs.txt.gz"
    receipt=json.loads((source_dir/"source.json").read_text())
    if hashlib.sha256(source.read_bytes()).hexdigest()!=receipt["sha256_compressed"]:
        raise ValueError("MPC snapshot checksum mismatch")
    n=ambiguous=0
    with gzip.open(source,"rb") as records:
        for line in records:
            n+=1
            num,prov=line[:5].decode("ascii","replace").strip(),line[5:12].decode("ascii","replace").strip()
            ids=keys.get(prov,set()) if prov else keys.get(num,set())
            if len(ids)==1:
                grouped[next(iter(ids))].append(line)
            elif len(ids)>1:
                ambiguous+=1
    receipt.update(total_input_lines=n,ambiguous_input_lines=ambiguous,matched_ids=len(grouped),
                   catalogue_url=caturl,catalogue_sha256=cathash,
                   kernel_sha256=hashlib.sha256((Path.home()/".find_orb/linux_p1550p2650.430t").read_bytes()).hexdigest())
    print("INPUT", json.dumps(receipt), flush=True)
    (OUT/"inputs.json").write_text(json.dumps(receipt,indent=2)+"\n")
    # Examine strongly hyperbolic and usable recent arcs first; every row receives a disposition.
    rows.sort(key=lambda r: -((float(r["e"])-1)/float(r["q"])))
    ids = [(r,grouped[r["spkid"]]) for r in rows[shard::parts]]
    output = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        for index, obj in enumerate(pool.map(fit,ids),1):
            output.append(obj)
            if index <= 12 or index%100 == 0:
                print("PROGRESS", index, json.dumps(obj),flush=True)
    fields = list(dict.fromkeys(k for r in output for k in r))
    with (OUT/"refits.csv").open("w",newline="") as fp:
        w = csv.DictWriter(fp,fieldnames=fields,lineterminator="\n")
        w.writeheader()
        w.writerows(output)
    states = collections.Counter(r["status"] for r in output)
    print("STATUS",dict(states),flush=True)
    for p in (OUT/"inputs.json",OUT/"refits.csv"):
        # The compressed tables can be deposited from a read-only runner.
        data = gzip.compress(p.read_bytes(),mtime=0)
        encoded=base64.b64encode(data).decode()
        for start in range(0,len(encoded),4000):
            print("ARCHPART",p.name+".gz",start,encoded[start:start+4000])


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-source",type=Path)
    parser.add_argument("--source",type=Path)
    parser.add_argument("--shard",type=int,default=0)
    parser.add_argument("--parts",type=int,default=1)
    args=parser.parse_args()
    if args.download_source:
        download_source(args.download_source)
    elif not 0<=args.shard<args.parts:
        parser.error("invalid shard number")
    elif args.source:
        main(args.source,args.shard,args.parts)
    else:
        with tempfile.TemporaryDirectory() as tmp:
            download_source(Path(tmp))
            main(Path(tmp),args.shard,args.parts)
