"""Save the input identities and log for a replay of a deposited observation."""
import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile

import campaign
import ingest_primary as ingest


def run(obsid):
    if (os.getenv("GITHUB_ACTIONS")!="true" or os.getenv("GITHUB_REPOSITORY")!=ingest.REPO or
            os.getenv("GITHUB_REF")!="refs/heads/"+ingest.BRANCH):
        raise RuntimeError("Use the Chandra research workflow")
    source=campaign.EVIDENCE/"observations"/f"{obsid}.json.gz"
    before=ingest.sha(source)
    ledger=campaign.read_obs(source)
    if ledger["obsid"]!=obsid or not ledger["products"]:
        raise ValueError("Missing or mismatched reference observation")
    root=campaign.EVIDENCE/"replays"/str(obsid)/(os.environ["GITHUB_RUN_ID"]+"-"+os.environ["GITHUB_RUN_ATTEMPT"])
    root.mkdir(parents=True,exist_ok=True)
    cmd=[sys.executable,"-u",str(campaign.HERE/"catalog_pack.py"),"replay",str(obsid)]
    started=datetime.now(timezone.utc).isoformat()
    with tempfile.TemporaryFile() as log:
        proc=subprocess.run(cmd,stdout=log,stderr=subprocess.STDOUT)
        log.seek(0)
        data=log.read()
    if ingest.sha(source)!=before:
        raise RuntimeError("The reference observation changed on disk")
    path=root/"replay.log.gz"
    path.write_bytes(gzip.compress(data,mtime=0))
    cov=[{"filename":p["input"]["filename"],"url":p["input"]["url"],"bytes":p["input"]["bytes"],
          "sha256":p["input"]["sha256"],"fits_sha256":p["input"]["fits_sha256"],
          "real_rows":len(p["search"]["triggers"]),"scrambles":len(p["search"]["controls"]),
          "injections":len(p["search"]["injections"])} for p in ledger["products"]]
    campaign.dump(root/"replay.json",{"obsid":obsid,"started_utc":started,
        "checked_commit":subprocess.check_output(["git","rev-parse","HEAD"],text=True).strip(),
        "observation_ledger_sha256":before,"inputs":cov,"command":"python results/chandra_acis_fxt/catalog_pack.py replay "+str(obsid),
        "implementation":campaign.identity(),"capture_program_sha256":ingest.sha(Path(__file__)),
        "replay_program_sha256":ingest.sha(campaign.HERE/"reproduce.py"),
        "catalog_restore_program_sha256":ingest.sha(campaign.HERE/"catalog_pack.py"),
        "uncompressed_log_sha256":hashlib.sha256(data).hexdigest(),"log_gzip_sha256":ingest.sha(path),
        "status":"reproduced" if proc.returncode==0 else "failed","exit_code":proc.returncode})
    print(data.decode(errors="replace")[-3400:],flush=True)
    return proc.returncode


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("obsid",type=int)
    raise SystemExit(run(p.parse_args().obsid))
