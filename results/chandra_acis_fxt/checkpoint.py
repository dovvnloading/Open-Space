"""Deposit available primary evidence in bounded Git pushes and publish measured progress."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

import campaign
import ingest_primary as ingest
import progress
import report


def command(*args):
    subprocess.run([str(a) for a in args],check=True)


def push(paths, message):
    command("git","add","--",*paths)
    result=subprocess.run(["git","diff","--cached","--quiet"])
    if result.returncode==1:
        command("git","commit","-qm",message)
    elif result.returncode!=0:
        raise RuntimeError("Cannot inspect staged evidence")
    command("git","push","origin","HEAD:refs/heads/"+ingest.BRANCH)


def collect():
    root=campaign.EVIDENCE
    path=campaign.HERE/"ingest_primary.py"
    command(sys.executable,path,"index")
    index=campaign.read(root/ingest.HISTORY/"source.json")
    current=index["artifacts"]
    groups=[str(n) for n in range(192) if f"chandra-evidence-{n}" in current
            and not ingest.receipt_file(root,str(n)).exists()]
    if ingest.SUP in current and not ingest.receipt_file(root,"supplement").exists():
        groups.insert(0,"supplement")
    # Keep related replies together while limiting each Git pack sent to the remote.
    batch=[]
    compressed=0
    for group in groups:
        obj=ingest.SUP if group=="supplement" else "chandra-evidence-"+group
        size=current[obj]["size_in_bytes"]
        if batch and (compressed+size>145*1024**2 or len(batch)>=10):
            deposit_group(path,root,batch)
            batch,compressed=[],0
        batch.append(group)
        compressed+=size
    if batch:
        deposit_group(path,root,batch)
    command(sys.executable,path,"audit","--allow-partial")
    s=progress.run(root)
    hist=ingest.HISTORY.as_posix()
    headline=s["catalog_headline"]
    early=f"evidence/provisional/moderate-count-v1/{ingest.RUN}-{ingest.ATTEMPT}"
    page=campaign.HERE/"README.md"
    lines=["# Chandra ACIS fast X-ray transient search","",f"**{headline}**","",
        f"The deposited scans produced {s['real_trigger_rows']:,} real trigger rows. "
        f"{s['detector_pass_rows']:,} passed the instrumental screen before counterpart review "
        f"and cross-product deduplication. Required catalog checks and photometry finished "
        f"in {s['complete_followup_obsids']:,} observations. Worker completion depends on field complexity, "
        "catalog response times and calibration. A global discovery count awaits the full search and review.","",
        f"[Current exposure by duration, coverage and input checksums](evidence/{hist}/PROGRESS.md) · "
        "[Measured priority flares](PRIORITY.md) · [Selection and reproduction commands](METHODS.md)","",
        "The prespecified 128-observation cohort contributed 1.5912 Ms of eligible 60-second start time. "
        "Two of 18 real trigger rows passed the instrument cuts. Both had one equal-or-larger instrument-stage "
        "maximum in 24 time scrambles (permutation estimate 0.08). Of 1,297 attempted synthetic flares, "
        f"939 were recovered. [Frozen cohort results]({early}/summary.json) · [Counterpart review](early_notes.md)","",
        "The separate sample of 14 observations with unknown release dates contained 175 real trigger rows, "
        "zero retained events, and 369.4 ks of eligible 60-second start time. One control event remained "
        "across 24 scrambles, and 139 of 165 attempted injections were recovered. "
        f"[Undated exposure and evidence](evidence/supplemental_search/{ingest.RUN}-{ingest.ATTEMPT}/README.md)","",
        "[Prospectively registered held-out test](HELD_OUT.md) · [Catalog reply extraction](CATALOG_STORAGE.md)"]
    page.write_text("\n".join(lines)+"\n")
    report.update_catalog(headline,"Measured")
    push([root,page,"README.md"],"Record measured Chandra coverage and reconcile catalog totals")


def deposit_group(program, root, batch):
    label=",".join(batch)
    print("Importing and verifying primary evidence groups "+label,flush=True)
    command(sys.executable,"-u",program,"batch","--groups",label,"--pack-catalogs")
    push([root],"Deposit checked Chandra primary evidence groups "+label)


def csc():
    root=campaign.EVIDENCE
    out=root/"csc_variability"
    pending=[]
    for path in sorted((root/"observations").glob("*.json.gz")):
        oid=int(path.name.split(".")[0])
        completed=out/"observations"/(str(oid)+".json.gz")
        row=campaign.read_obs(path)
        if row["status"] not in {"complete","followup_incomplete"}:
            continue
        if completed.is_file():
            old=campaign.read_obs(completed)
            if old["ledger_sha256"]==ingest.sha(path) and all(v.get("status")=="ok" for v in old["query_receipts"].values()):
                continue
        pending.append(oid)
    if not pending:
        return
    # Static catalogs and per-observation results live under a single cache.
    # Keep each independent batch summary, including failed or empty responses.
    for i in range(0,len(pending),200):
        ids=pending[i:i+200]
        command(sys.executable,"-u",campaign.HERE/"csc_variability.py","--evidence",root,
                "--out",out,"--obsids",",".join(map(str,ids)))
        result=out/"csc_audit.json"
        if result.exists():
            batches=out/"batches"
            batches.mkdir(exist_ok=True)
            (batches/f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}-{i//200}.json").write_bytes(result.read_bytes())
    push([out],"Deposit same-observation CSC variability checks in real and scrambled triggers")


def optical():
    root = campaign.EVIDENCE
    out = root/"lead_optical"
    command(sys.executable, "-u", campaign.HERE/"lead_optical.py", "--from-ledgers",
            "--evidence", root, "--out", out)
    push([out], "Deposit optical and infrared counterpart measurements for retained real triggers")


if __name__=="__main__":
    if (os.environ.get("GITHUB_ACTIONS")!="true" or os.environ.get("GITHUB_REPOSITORY")!=ingest.REPO
            or os.environ.get("GITHUB_REF")!="refs/heads/"+ingest.BRANCH):
        raise RuntimeError("Use the Chandra workflow")
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("operation",choices=("collect","csc","optical"))
    operations = {"collect": collect, "csc": csc, "optical": optical}
    operations[parser.parse_args().operation]()
