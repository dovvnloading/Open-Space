"""Deposit exact EVT2 transport bytes only for manually retained, published candidates."""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from urllib.parse import unquote, urlsplit

import archive
import campaign
import publication

PART = 24 * 1024**2
BASE = campaign.HERE / "evidence/candidate_inputs"
FORMAT = "chandra-candidate-inputs-v1"
DIGEST = re.compile(r"[a-f0-9]{64}")
FILENAME = re.compile(r"acisf(\d+)(?:_[a-zA-Z0-9]+)*?N(\d+)[a-zA-Z0-9_]*_evt2\.fits(?:\.gz)?",re.I)


def need(test, detail):
    if not test:
        raise ValueError(detail)


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def encoded(o):
    return (json.dumps(o, sort_keys=True, indent=2, allow_nan=False)+"\n").encode()


def receipt(src):
    name = src["filename"]
    match = FILENAME.fullmatch(name)
    need(match and int(match[1]) == int(src["obsid"]) and src["bytes"] > 0 and "revision" in src and
         DIGEST.fullmatch(src["sha256"]) and DIGEST.fullmatch(src["fits_sha256"]), "Incomplete native EVT2 receipt")
    archive.check_archive_url(src["url"])
    need(unquote(urlsplit(src["url"]).path.split("/")[-1]) == name,"Receipt URL does not identify the original filename")
    return int(match[2])


def verify(path, src):
    receipt(src)
    need(path.name == src["filename"] and path.stat().st_size == src["bytes"] and sha(path) == src["sha256"],
         "EVT2 archive name, length or SHA-256 changed")
    h = archive.inspect_product(path,int(src["obsid"]))
    for key in ("fits_sha256","revision","rows","readmode","grating"):
        need(key in src and h[key] == src[key],f"Reconstructed EVT2 {key} differs from mission receipt")


def split(path, src, output):
    target = output / "objects" / src["sha256"]
    target.mkdir(parents=True)
    parts=[]
    with path.open("rb") as f:
        while raw := f.read(PART):
            part_hash=hashlib.sha256(raw).hexdigest()
            filename=f"{len(parts):08d}-{part_hash}.evt2part"
            (target/filename).write_bytes(raw)
            parts.append({"file":str((target/filename).relative_to(output)),"sha256":part_hash,"bytes":len(raw)})
    need(parts and sum(p["bytes"] for p in parts)==src["bytes"],"Input length changed during deposit")
    return parts


def restore(entry, store, directory):
    src, parts = entry["archive_receipt"], entry["parts"]
    receipt(src)
    directory.mkdir(parents=True,exist_ok=True)
    target=directory/src["filename"]
    if target.exists():
        verify(target,src)
        return target
    with tempfile.TemporaryDirectory(dir=directory) as tmp:
        path=Path(tmp)/target.name
        with path.open("xb") as out:
            for i,p in enumerate(parts):
                h=p["sha256"]
                need(DIGEST.fullmatch(h) and p["file"]==f"objects/{src['sha256']}/{i:08d}-{h}.evt2part" and
                     0 < p["bytes"] <= PART and (i == len(parts)-1 or p["bytes"] == PART),"Invalid or unordered part manifest")
                blob=store/p["file"]
                need(blob.resolve().is_relative_to(store.resolve()) and blob.stat().st_size==p["bytes"],"Part location or length changed")
                raw=blob.read_bytes()
                need(hashlib.sha256(raw).hexdigest()==h,"Candidate part lost its original bytes")
                out.write(raw)
        verify(path,src)
        try:
            os.link(path,target)
        except FileExistsError:
            verify(target,src)
    return target


def manifest(store):
    index=store/"manifest.json"
    data=campaign.read(index)
    need(data["format"]==FORMAT and data["part_bytes_max"]==PART,"Unsupported candidate archive format")
    addressed=store/"manifests"/(sha(index)+".json")
    need(addressed.is_file() and sha(addressed)==sha(index),"Content-addressed manifest is missing or changed")
    return data


def capture(src, retained, previous, work):
    notes=[]
    h, name = src["sha256"], src["filename"]
    if previous and h in previous["products"]:
        try:
            path=restore(previous["products"][h],BASE,work)
            verify(path,src)
            return path,{"kind":"previous_verified_deposit","prior_manifest_sha256":sha(BASE/"manifest.json")},notes
        except (OSError,ValueError,KeyError,EOFError) as exc:
            notes.append({"source":"previous_deposit","error":str(exc)})
    for folder in retained:
        folder=folder.resolve()
        if not folder.is_dir():
            notes.append({"source":folder.name,"error":"retained artifact directory unavailable"})
            continue
        for path in sorted(folder.rglob(name)):
            if path.is_symlink() or not path.is_file() or not path.resolve().is_relative_to(folder):
                continue
            part_id={"artifact_root":folder.name,"member":str(path.relative_to(folder))}
            try:
                verify(path,src)
                return path,{"kind":"retained_actions_artifact",**part_id},notes
            except (OSError,ValueError,KeyError,EOFError) as exc:
                notes.append({**part_id,"error":str(exc)})
    pinned=work/"pinned-receipt.json"
    oid=int(src["obsid"])
    urls=[src["url"],*(f"{base}/{oid%10}/{oid}/primary/{name}" for base in archive.ROOTS)]
    for url in dict.fromkeys(urls):
        try:
            archive.dump(pinned,[{**src,"url":url}])
            archive.replay(pinned,work)
            path=work/name
            verify(path,src)
            return path,{"kind":"pinned_archive_replay","url":url},notes
        except (OSError,ValueError,KeyError,EOFError,archive.http.client.HTTPException) as exc:
            notes.append({"source":url,"error":str(exc)})
    raise RuntimeError(f"No exact retained or official-archive copy of {name}: {notes}")


def reviewed(decisions, supplement):
    pub=campaign.read(campaign.EVIDENCE/"publication.json")
    document=campaign.read(decisions)
    need(pub["decisions_sha256"]==sha(decisions) and document["reviewed_inputs_sha256"]==pub["reviewed_inputs_sha256"],
         "The manual decisions do not match the finalized publication receipt")
    all_ids=[publication.key(d) for d in document["decisions"]]
    chosen=[publication.key(d) for d in document["decisions"] if d["decision"]=="candidate"]
    counts=pub["dispositions"]
    need(len(set(all_ids))==len(all_ids) and len(chosen)==sum(g.get("candidate",0) for g in counts.values()) and
         len(all_ids)==sum(sum(g.values()) for g in counts.values()),"Publication and manual counts disagree")
    if not chosen:
        return {},pub
    _,leads,sig,_=publication.load(supplement)
    confirmed,_=publication.adjudicate(document,leads,sig)
    products={}
    for k in chosen:
        need(confirmed[k][0]=="candidate","Unreviewed input requested")
        r=leads[k]
        ledger=r["group"]["root"]/"observations"/f"{k[0]}.json.gz"
        obj=campaign.read_obs(ledger)
        src,=[p for p in obj["acquired"] if p["sha256"]==k[1]]
        actual,=[p["input"] for p in obj["products"] if p["input"]["sha256"]==k[1]]
        need(src==actual and src["filename"]==r["lead"]["input_name"] and src["url"]==r["receipt"]["url"]
             and src["fits_sha256"]==r["receipt"]["fits_sha256"] and int(src["bytes"])==int(r["receipt"]["bytes"]),
             "The published product differs from the acquired EVT2 receipt")
        entry=products.setdefault(k[1],{"archive_receipt":src,"filename_revision":receipt(src),"candidates":[],
            "ledger":str(ledger.relative_to(campaign.HERE)),"ledger_sha256":sha(ledger),"cohort":r["group"]["name"]})
        need(entry["archive_receipt"]==src and entry["ledger_sha256"]==sha(ledger),"Conflicting receipts for one checksum")
        entry["candidates"].append({"obsid":k[0],"input_sha256":k[1],"trigger_key":k[2]})
    return products,pub


def preserve(decisions,supplement,retained):
    products,pub=reviewed(decisions,supplement)
    previous=manifest(BASE) if BASE.exists() else None
    publication_hash=sha(campaign.EVIDENCE/"publication.json")
    if (previous and previous["publication_sha256"]==publication_hash and
            previous["decisions_sha256"]==sha(decisions) and set(previous["products"])==set(products)):
        try:
            for h,entry in products.items():
                need(all(entry[k]==previous["products"][h][k] for k in entry),"Published candidate provenance changed")
                with tempfile.TemporaryDirectory() as tmp:
                    restore(previous["products"][h],BASE,Path(tmp))
            print(f"Verified existing candidate store, {len(products)} unique products.")
            return previous
        except (OSError, ValueError, KeyError, EOFError):
            # Rebuild a damaged copy without replacing the current store until
            # every product has passed a complete reassembly and FITS check.
            pass
    repo_prefix="results/chandra_acis_fxt"
    BASE.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=BASE.parent) as tmp:
        temp=Path(tmp)
        stage=temp/"new"
        stage.mkdir()
        for h,entry in sorted(products.items()):
            src=entry["archive_receipt"]
            with tempfile.TemporaryDirectory(dir=temp) as work:
                raw,method,attempts=capture(src,retained,previous,Path(work))
                entry.update(parts=split(raw,src,stage),acquisition=method,unavailable_sources=attempts,
                             transport_encoding="gzip" if src["filename"].endswith(".gz") else "identity")
                restore(entry,stage,Path(work)/"verified")
            filename=src["filename"]
            entry["concatenate_from_repository_root"]=(f"LC_ALL=C cat -- {repo_prefix}/evidence/candidate_inputs/objects/{h}/"
                                                       f"[0-9]*.evt2part > \"$RUNNER_TEMP/{filename}\"")
            entry["restore_and_verify_command"]=f'python {repo_prefix}/candidate_inputs.py restore {h} "$RUNNER_TEMP/restored"'
            entry["exact_product_search_replay_command"]=f'python {repo_prefix}/candidate_inputs.py replay {h} "$RUNNER_TEMP/replay"'
        obj={"format":FORMAT,"part_bytes_max":PART,"products":products,"candidate_count":sum(len(e["candidates"]) for e in products.values()),
             "unique_product_count":len(products),"total_transport_bytes":sum(e["archive_receipt"]["bytes"] for e in products.values()),
             "reviewed_inputs_sha256":pub["reviewed_inputs_sha256"],"decisions_sha256":sha(decisions),
             "publication_sha256":publication_hash,"preservation_run_id":os.environ.get("GITHUB_RUN_ID"),
             "preservation_run_attempt":os.environ.get("GITHUB_RUN_ATTEMPT"),
             "package_implementation_sha256":sha(Path(__file__))}
        data=encoded(obj)
        (stage/"manifests").mkdir()
        (stage/"manifests"/(hashlib.sha256(data).hexdigest()+".json")).write_bytes(data)
        (stage/"manifest.json").write_bytes(data)
        manifest(stage)
        if BASE.exists():
            BASE.rename(temp/"old")
        try:
            stage.rename(BASE)
        except BaseException:
            if (temp/"old").exists():
                (temp/"old").rename(BASE)
            raise
    print(f"{len(products)} candidate products, {obj['total_transport_bytes']:,} archive bytes, manifest {sha(BASE/'manifest.json')}")
    return obj


def replay(entry,directory):
    from search import search_file
    src=entry["archive_receipt"]
    path=restore(entry,BASE,directory)
    ledger=publication.local(entry["ledger"])
    need(sha(ledger)==entry["ledger_sha256"],"Source ledger changed")
    p,=[p for p in campaign.read_obs(ledger)["products"] if p["input"]["sha256"]==src["sha256"]]
    old=copy.deepcopy(p["search"])
    seed=int(hashlib.sha256(("chandra-acis-fxt-v1:"+src["sha256"]).encode()).hexdigest()[:16],16)>>1
    new=search_file(path,seed=seed,n_controls=campaign.N_CONTROL,n_injections=campaign.N_INJECT)
    for result in (old,new):
        for scan in [result,*result["controls"]]:
            for event in scan["triggers"]:
                for annotation in ("obsid","catalog_cache","vet","flux","diagnostics","trigger_key"):
                    event.pop(annotation,None)
    need(new==old,"Complete raw triggers, controls or injection results differ on exact stored EVT2")
    print(f"Reproduced raw search, all {campaign.N_CONTROL} scrambles and {campaign.N_INJECT} injection trials for {src['filename']}.")


def check_index():
    m=manifest(BASE)
    expected={BASE/"manifest.json", BASE/"manifests"/(sha(BASE/"manifest.json")+".json")}
    for entry in m["products"].values():
        for piece in entry["parts"]:
            p=BASE/piece["file"]
            need(p.stat().st_size==piece["bytes"]<=PART and sha(p)==piece["sha256"],"Candidate part verification failed")
            expected.add(p)
    repo=campaign.HERE.parents[1]
    def git(*args):
        return subprocess.check_output(["git","-C",str(repo),*args])
    rows=git("ls-files","--stage","-z","--",str(BASE.relative_to(repo))).split(b"\0")
    staged={}
    for row in filter(None,rows):
        meta,name=row.split(b"\t",1)
        mode,oid,index=meta.decode().split()
        need(mode=="100644" and index=="0","Unexpected file mode or unresolved index stage in byte store")
        staged[repo/os.fsdecode(name)]=oid
    need(set(staged)==expected,"Some candidate bytes are ignored, unstaged, stale or absent")
    for path,oid in staged.items():
        need(git("hash-object","--no-filters",str(path)).decode().strip()==oid,
             "Staged Git blob differs from original bytes, possibly a filter or LFS pointer")
    print(f"Git contains {len(staged)} verified blobs in the candidate store; total native event bytes {m['total_transport_bytes']:,}. "
          "The repository verifier materializes and copies Git blobs from the entire committed tree. "
          "Allow memory for several times the tracked tree size. A small individual part does not limit the entire push size; "
          "batch object directories before the final manifest when needed. Run both repository gates on the final commit.")


if __name__=="__main__":
    need(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")=="dovvnloading/Open-Space" and
         os.getenv("GITHUB_REF")=="refs/heads/research/chandra-acis-fast-transients","Use the Chandra research workflow")
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest="mode",required=True)
    c=commands.add_parser("preserve")
    c.add_argument("--decisions",required=True)
    c.add_argument("--supplement",required=True)
    c.add_argument("--retained",type=Path,action="append",default=[])
    for mode in ("restore","replay"):
        c=commands.add_parser(mode)
        c.add_argument("sha256")
        c.add_argument("directory",type=Path)
    commands.add_parser("check-index")
    args=parser.parse_args()
    if args.mode=="preserve":
        preserve(publication.local(args.decisions),(campaign.HERE/args.supplement).resolve(),args.retained)
    elif args.mode=="check-index":
        check_index()
    else:
        entry=manifest(BASE)["products"][args.sha256]
        if args.mode=="restore":
            print(restore(entry,BASE,args.directory))
        else:
            replay(entry,args.directory)
