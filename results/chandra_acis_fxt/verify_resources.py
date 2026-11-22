"""Time the original committed-tree verifier or its tests, preserving logs and reproducing their exit code."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys

HERE=Path(__file__).resolve().parent
REPO=HERE.parents[1]
COMMANDS = {"verify":[sys.executable,"-m","open_space.verify"],
            "tests":[sys.executable,"-m","unittest","discover","-s","tests"],
            "survey-tests":[sys.executable,"-m","unittest","discover","-s",str(HERE),"-p","test_*.py"]}


def tree_size():
    by_path = {"all":{"files":0,"bytes":0}, "survey":{"files":0,"bytes":0}, "survey_json_gz":{"files":0,"bytes":0}}
    distinct={}
    proc = subprocess.Popen(["git","ls-tree","-r","-l","-z","HEAD"],stdout=subprocess.PIPE,cwd=REPO)
    tail=b""
    for block in iter(lambda:proc.stdout.read(1<<16),b""):
        parts=(tail+block).split(b"\0")
        tail=parts.pop()
        for record in parts:
            meta,name=record.split(b"\t",1)
            _,typ,oid,size=meta.split()
            if typ != b"blob":
                continue
            size=int(size)
            distinct[oid]=size
            survey=name.startswith(b"results/chandra_acis_fxt/")
            for group in ["all"]+(["survey"] if survey else [])+(["survey_json_gz"] if survey and name.endswith(b".json.gz") else []):
                by_path[group]["files"]+=1
                by_path[group]["bytes"]+=size
    if proc.wait() or tail:
        raise ValueError("Incomplete committed tree inventory")
    return {"by_tracked_path":by_path,"distinct_git_blobs":len(distinct),"distinct_blob_bytes":sum(distinct.values())}


def machine():
    limits={}
    for n in ["memory.current","memory.max","memory.peak","memory.swap.max","memory.events","memory.events.local"]:
        p=Path("/sys/fs/cgroup")/n
        if p.is_file():
            limits[n]=p.read_text().strip()
    volume=shutil.disk_usage(REPO)
    return {"cgroup":limits,"disk_free_bytes":volume.free,"disk_total_bytes":volume.total,
            "utc":datetime.now(timezone.utc).isoformat(),"nproc_available":len(os.sched_getaffinity(0))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("operation",choices=COMMANDS)
    p.add_argument("output",type=Path)
    a=p.parse_args()
    if os.getenv("GITHUB_ACTIONS")!="true" or os.getenv("GITHUB_REPOSITORY")!="dovvnloading/Open-Space":
        raise ValueError("Use the Chandra Actions workflow")
    output=a.output.resolve()
    if output.is_relative_to(REPO):
        raise ValueError("Keep resource receipts in RUNNER_TEMP and upload an Actions artifact")
    output.mkdir(parents=True,exist_ok=True)
    before=machine()
    receipt={"operation":a.operation,"command":COMMANDS[a.operation],"run_id":os.getenv("GITHUB_RUN_ID"),
             "attempt":os.getenv("GITHUB_RUN_ATTEMPT"),"checked_head":subprocess.check_output(
             ["git","rev-parse","HEAD"],cwd=REPO,text=True).strip(),"event_head":os.getenv("GITHUB_SHA"),
             "platform":platform.platform(),"tree":tree_size(),"before":before}
    report=output/(a.operation+".json")
    report.write_text(json.dumps(receipt,indent=2)+"\n")
    timing=output/(a.operation+".time.txt")
    proc = subprocess.run(["/usr/bin/time","-v","-o",str(timing),*COMMANDS[a.operation]],cwd=REPO)
    receipt.update(after=machine(),returncode=proc.returncode)
    report.write_text(json.dumps(receipt,indent=2)+"\n")
    raise SystemExit(proc.returncode)


if __name__ == "__main__":
    main()
