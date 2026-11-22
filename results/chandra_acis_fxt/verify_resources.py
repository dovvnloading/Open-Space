"""Time the original committed-tree verifier or its tests, preserving logs and reproducing their exit code."""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile

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
    legacy=Path("/sys/fs/cgroup/memory")
    for name in ("memory.limit_in_bytes","memory.max_usage_in_bytes","memory.usage_in_bytes",
                 "memory.memsw.limit_in_bytes","memory.failcnt"):
        p=legacy/name
        if p.is_file():
            limits["v1/"+name]=p.read_text().strip()
    meminfo={}
    for line in Path("/proc/meminfo").read_text().splitlines():
        key,value=line.split(":",1)
        if key in {"MemTotal","MemAvailable","SwapTotal","SwapFree"}:
            meminfo[key]=value.strip()
    volume=shutil.disk_usage(REPO)
    return {"cgroup":limits,"proc_meminfo":meminfo,"disk_free_bytes":volume.free,"disk_total_bytes":volume.total,
            "utc":datetime.now(timezone.utc).isoformat(),"nproc_available":len(os.sched_getaffinity(0))}


def reserve_for_tree(tree, before):
    """Provision temporary swap only when the measured verifier footprint approaches runner capacity."""
    gib = 1024**3
    size = tree["by_tracked_path"]["all"]["bytes"]
    # Run 45008126305 used 9,643,280 KiB for a 2,432,852,215-byte Git tree.
    # Git's batch output and the per-path Python bytes coexist in the unmodified verifier.
    expected = (size*9+1)//2 + gib
    info = before["proc_meminfo"]
    available = sum(int(info[key].split()[0])*1024 for key in ("MemAvailable","SwapFree"))
    record = {"estimated_bytes_including_headroom":expected,"available_ram_and_swap_bytes":available,
              "calibration_run":45008126305}
    if available >= expected:
        return {**record,"status":"sufficient"}
    if (os.getenv("GITHUB_ACTIONS") != "true" or
        os.getenv("GITHUB_WORKFLOW_REF","").split("@",1)[0] != "dovvnloading/Open-Space/.github/workflows/chandra-acis-fxt.yml" or
        os.getenv("GITHUB_REF") != "refs/heads/research/chandra-acis-fast-transients" or not os.getenv("RUNNER_TEMP")):
        return {**record,"status":"outside_research_runner"}
    extra = min(24, math.ceil((expected-available)/gib)+2)*gib
    if before["disk_free_bytes"] < extra+12*gib:
        return {**record,"status":"insufficient_free_disk"}
    fd, name = tempfile.mkstemp(prefix="chandra-verify-swap-",dir=os.environ["RUNNER_TEMP"])
    os.close(fd)
    try:
        for cmd in (["fallocate","-l",str(extra),name],["sudo","-n","mkswap",name],
                    ["sudo","-n","swapon",name]):
            r = subprocess.run(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=False)
            if r.returncode:
                Path(name).unlink(missing_ok=True)
                return {**record,"status":"unavailable","operation":cmd[-2] if cmd[0]=="sudo" else cmd[0],
                        "returncode":r.returncode,"detail":r.stdout[-800:]}
    except OSError as error:
        Path(name).unlink(missing_ok=True)
        return {**record,"status":"unavailable","detail":str(error)}
    # The hosted VM is destroyed after the job. Keeping this file active also covers its test run.
    return {**record,"status":"active_on_ephemeral_runner","additional_bytes":extra,"after":machine()}


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
    tree=tree_size()
    reserve=reserve_for_tree(tree,before) if a.operation in {"verify","tests"} else {"status":"unused"}
    receipt={"operation":a.operation,"command":COMMANDS[a.operation],"run_id":os.getenv("GITHUB_RUN_ID"),
             "attempt":os.getenv("GITHUB_RUN_ATTEMPT"),"checked_head":subprocess.check_output(
             ["git","rev-parse","HEAD"],cwd=REPO,text=True).strip(),"event_head":os.getenv("GITHUB_SHA"),
             "platform":platform.platform(),"tree":tree,"before":before,"memory_reserve":reserve}
    report=output/(a.operation+".json")
    report.write_text(json.dumps(receipt,indent=2)+"\n")
    timing=output/(a.operation+".time.txt")
    proc = subprocess.run(["/usr/bin/time","-v","-o",str(timing),*COMMANDS[a.operation]],cwd=REPO)
    receipt.update(after=machine(),returncode=proc.returncode)
    report.write_text(json.dumps(receipt,indent=2)+"\n")
    raise SystemExit(proc.returncode)


if __name__ == "__main__":
    main()
