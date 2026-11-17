"""Refit every frozen SBDB target above the explicit 1 km/s osculating threshold."""
import argparse
import base64
import bisect
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import date
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
from urllib.request import Request, urlopen

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from refit import aliases, cat, distance_and_direction, FITTER

MU_KM=132712440041.93938
AU_KM=149597870.7
KERNEL=Path.home()/".find_orb/linux_p1550p2650.430t"
KNOWN={"50788063":"1I","1003639":"2I","1004083":"3I"}
ENDPOINT="https://data.minorplanetcenter.net/api/get-obs"


def hashbytes(data):
    return hashlib.sha256(data).hexdigest()


def vinf(e,q):
    return math.sqrt(max(0.,float(e)-1)*MU_KM/(float(q)*AU_KM))


def matches(idstr, allowed):
    clean=idstr.replace(" ","")
    return clean in allowed or (len(clean)==8 and clean[0] in "ACDPXI" and clean[1:] in allowed) or (
        len(clean)==12 and (clean[:5] in allowed or clean[5:] in allowed))


def jd(line):
    try:
        year,month,day=line[15:32].split()
        d=float(day)
        return date(int(year),int(month),int(d)).toordinal()+1721424.5+d%1
    except (ValueError,OverflowError):
        return None


def api_if_missing(row, allowed):
    ident=row["pdes"]
    prefix=row["full_name"].strip().split(" (")[0]
    names=list(dict.fromkeys([KNOWN.get(str(row["spkid"])),prefix,ident]+sorted(allowed)))
    attempts=[]
    for name in filter(bool,names):
        request={"desigs":[name],"output_format":["OBS80"]}
        rec={"request":request}
        attempts.append(rec)
        try:
            with urlopen(Request(ENDPOINT,data=json.dumps(request).encode(),method="GET",
                                headers={"Content-Type":"application/json","User-Agent":"Open-Space orbit research/1"}),timeout=25) as resp:
                body=resp.read()
                rec.update(response_sha256=hashbytes(body),date=resp.headers.get("Date"),status=resp.status)
            lines=[l[:80].encode("ascii")+b"\n" for l in json.loads(body)[0]["OBS80"].splitlines()
                   if len(l)>=80 and matches(l[:12],allowed)]
            if len(lines)>=3:
                return lines,attempts
        except Exception as exc:
            rec["error"]=str(exc)
    return [],attempts


def run_one(row, raw):
    allowed=aliases(row)
    sk=str(row["spkid"])
    out={"spkid":sk,"pdes":row["pdes"],"name":row["full_name"].strip(),
         "catalogue_vinf_kms":vinf(row["e"],row["q"]),"catalogue_e":float(row["e"]),
         "catalogue_first":row["first_obs"],"catalogue_last":row["last_obs"],"catalogue_n":row["n_obs_used"],
         "catalogue_non_grav":{k:row[k] for k in ("A1","A2","A3","DT","S0") if row[k] is not None},
         "aliases":sorted(allowed),"bulk_input_lines":len(raw),"source":"frozen CmtObs"}
    if len(raw)<3:
        frozen=HERE.parent.parent/"orbits/data"/(KNOWN.get(sk,"unavailable")+".obs.gz")
        if frozen.is_file():
            data=gzip.decompress(frozen.read_bytes())
            raw=[line for line in data.splitlines(keepends=True) if matches(line[:12].decode("ascii"),allowed)]
            out.update(source="archived individual astrometry",individual_source=frozen.name,individual_sha256=hashbytes(data))
        if len(raw)<3:
            raw,attempts=api_if_missing(row,allowed)
            out.update(source="MPC get-obs fallback",api=ENDPOINT,attempts=attempts)
    content=b"".join(raw)
    times=defaultdict(list)
    for line in raw:
        t=jd(line)
        if t is not None:
            times[line[77:80].decode("ascii").strip()].append(t)
    alltimes=sorted(d for values in times.values() for d in values)
    for values in times.values():
        values.sort()
    out.update(input_lines=len(raw),raw_sha256=hashbytes(content),input_stations=len(times),
               input_first_jd=min(alltimes,default=None),input_last_jd=max(alltimes,default=None))
    if len(alltimes)<3:
        return {**out,"status":"insufficient associated observations"},content
    per=float(row["tp"])
    target=max(2287924.,per-300*365.25)
    cpu=120
    with tempfile.TemporaryDirectory() as tmp:
        stage=Path(tmp)
        inp=stage/"observations.mpc"
        staged=b"#interstellar\n"+content
        inp.write_bytes(staged)
        out.update(cpu_seconds=cpu,wall_limit_s=300,unpenalized=True,
                   staged_sha256=hashbytes(staged),requested_epoch=target)
        args=[str(FITTER),str(inp),"-i","-h-1","-tEjd"+str(target),"-O"+tmp,"-r120,122",
              "PERTURBERS=7fe","LINUX_JPL_FILENAME="+str(KERNEL),
              "JSON_SHORT_ELEMENTS="+tmp+"/short-%p.json","JSON_ELEMENTS_NAME="+tmp+"/long-%p.json",
              "JSON_COMBINED_NAME="+tmp+"/combined.json"]
        try:
            p=subprocess.run(args,cwd=tmp,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=300)
            out["returncode"]=p.returncode
            files=sorted(stage.glob("long-*.json"))
            if p.returncode or len(files)!=1:
                return {**out,"status":"cpu timeout" if p.returncode in (-9,-24,-15) else "no unique completed solution",
                        "element_files":[path.name for path in files],"diagnostic":p.stdout[-250:].decode(errors="replace")},content
            data=json.loads(files[0].read_bytes())
            if len(data["objects"])!=1:
                return {**out,"status":"multiple fitter identities"},content
            solved=next(iter(data["objects"].values()))
            el,obs=solved["elements"],solved["observations"]
            resid=obs["residuals"]
            missing=[]
            badcodes=set()
            for record in resid:
                if not matches(record["packed"],allowed):
                    badcodes.add(record["packed"])
                epochs=times.get(record["obscode"],[])
                k=bisect.bisect_left(epochs,record["JD"])
                if not any(abs(t-record["JD"])<=.00012 for t in epochs[max(0,k-1):k+1]):
                    missing.append((record["obscode"],record["JD"]))
            identified=(matches(solved["packed"],allowed) and not badcodes and not missing and len(resid)==obs["count"])
            r,inbound=distance_and_direction(el)
            de=el["e"]-1
            sig=el.get("e sigma")
            actual_span=obs["latest_used"]-obs["earliest_used"]
            raw_span=alltimes[-1]-alltimes[0]
            out.update(status="identity and observation epochs verified" if identified else "output identity or astrometry mismatch",
                       returned_name=solved["object"],returned_packed=solved["packed"],packed_residual_mismatches=sorted(badcodes),
                       station_time_mismatches=len(missing),mismatch_examples=missing[:4],
                       returned_epoch=el["epoch"],center=el["central body"],r_au=r,inbound=inbound,
                       e=el["e"],q_au=el["q"],sigma_e=sig,inv_a=-de/el["q"],
                       vinf_bary_kms=math.sqrt(1.001341838645975)*vinf(el["e"],el["q"]) if de>0 else None,
                       formal_e_excess_sigma=de/sig if sig and sig>0 else None,
                       rms_per_coordinate_arcsec=el["rms_residual"],weighted_rms=el["weighted_rms_residual"],
                       used=obs["used"],total=obs["count"],used_fraction=obs["used"]/obs["count"],
                       first_used_jd=obs["earliest_used"],last_used_jd=obs["latest_used"],
                       used_days=actual_span,input_days=raw_span,used_span_fraction=actual_span/raw_span if raw_span else None,
                       solution_sha256=hashbytes(files[0].read_bytes()))
            out["far_inbound_100au"]=bool(identified and abs(target-el["epoch"])<=.0001 and inbound and r>=100
                                         and (el["central body"]=="Object_-1" or "bary" in el["central body"].lower()))
            out["far_inbound_250au"]=bool(out["far_inbound_100au"] and r>=250)
        except subprocess.TimeoutExpired:
            out["status"]="wall timeout"
        except (OSError,ValueError,KeyError,ArithmeticError) as exc:
            out.update(status="unusable fitter output",diagnostic=str(exc))
    return out,content


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source",type=Path)
    parser.add_argument("shard",type=int,choices=range(6))
    parser.add_argument("out",type=Path)
    args=parser.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    rows,url,sha=cat()
    keys=defaultdict(set)
    for row in rows:
        for alias in aliases(row):
            keys[alias].add(str(row["spkid"]))
    selected=sorted((r for r in rows if vinf(r["e"],r["q"])>=1),key=lambda r:(-vinf(r["e"],r["q"]),str(r["spkid"])))
    assert len(selected)==67
    subset=selected[args.shard::6]
    ids={str(r["spkid"]) for r in subset}
    collected=defaultdict(list)
    source=args.source/"CmtObs.txt.gz"
    receipt=json.loads((args.source/"source.json").read_text())
    with source.open("rb") as fp:
        assert hashlib.file_digest(fp,"sha256").hexdigest()==receipt["sha256_compressed"]
    with gzip.open(source,"rb") as fp:
        for line in fp:
            compact=line[:12].decode("ascii","replace")
            num,prov=compact[:5].strip(),compact[5:].strip()
            # A fragment-bearing provisional ID takes priority over its parent number.
            found=keys.get(prov,set()) if prov else keys.get(num,set())
            if len(found)==1:
                sk=next(iter(found))
                if sk in ids:
                    collected[sk].append(line)
    receipt.update(snapshot_run_id=44029645452,snapshot_artifact_id=14304767729,
                   catalog_url=url,catalog_sha256=sha,selection_count=len(selected),shard=args.shard,parts=6,
                   fitter_sha256=hashbytes(FITTER.read_bytes()),kernel_sha256=hashbytes(KERNEL.read_bytes()),
                   flags=["-i","-h-1","-tEjd<target>","-r120,122","PERTURBERS=7fe","LINUX_JPL_FILENAME=<DE430t>"],
                   input_header="#interstellar",uses_catalogue_seed=False)
    records=[]
    observation_map={}
    with ThreadPoolExecutor(max_workers=2) as pool:
        fs=[pool.submit(run_one,r,collected[str(r["spkid"])]) for r in subset]
        for f in fs:
            result,raw=f.result()
            records.append(result)
            observation_map[result["spkid"]]=raw.decode("ascii")
            print("ROW",result["pdes"],result["status"],result.get("vinf_bary_kms"),result.get("formal_e_excess_sigma"),flush=True)
    assert len(records)==len(subset)
    report={"receipt":receipt,"status_counts":dict(Counter(r["status"] for r in records)),"records":records}
    payloads={f"shard-{args.shard}.json":json.dumps(report,indent=2)+"\n",
              f"observations-{args.shard}.json":json.dumps(observation_map,separators=(",",":"))+"\n"}
    for filename,body in payloads.items():
        raw=body.encode()
        if filename.startswith("observations-"):
            (args.out/(filename+".gz")).write_bytes(gzip.compress(raw,mtime=0))
        else:
            (args.out/filename).write_bytes(raw)
        code=base64.b64encode(gzip.compress(raw,mtime=0)).decode()
        for i in range(0,len(code),2700):
            print(f"PUBLIC_DATA {filename} {i//2700} {code[i:i+2700]}")


if __name__=="__main__":
    main()
