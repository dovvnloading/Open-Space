"""Fixed-position follow-up of an exploratory flare in the GW170817 field.

Inside the Chandra workflow run field_followup.py after checkpoint collection,
CSC matching and optical matching. Selection does not change blind trial counts.
"""
from argparse import ArgumentParser, Namespace
import os
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory

import alerts
import campaign
import field_replay
import ingest_primary as ingest
import priority
import target_history
from archive import dump, replay
from vet import CAMPAIGN_REF


OBSID = 27089
PRODUCT = "dfc736671520bc9f5d88db7df50d02bdb21434803ec8643925057c0315fa5f2f"
ARTIFACT = 14747454573
ZIP_SHA = "sha256:ab9b7930c7b5aebb21f716fd590a78022974090dc97f373997b6ea222f8d6cb3"
TRIGGER = "real/0"


def select(root):
    receipt = campaign.read(ingest.receipt_file(root, "80"))
    members = root / ingest.HISTORY / "members/80.json.gz"
    if (receipt["group"], receipt["artifact_id"], receipt["zip_digest"], receipt["code_commit"],
            receipt["member_manifest_sha256"]) != ("80", ARTIFACT, ZIP_SHA, ingest.CODE, ingest.sha(members)):
        raise ValueError("Group 80 producer receipt does not match the reviewed evidence")
    ledger = root / "observations" / f"{OBSID}.json.gz"
    relative = ledger.relative_to(root).as_posix()
    original, = [m for m in campaign.read_obs(members) if m["source"] == relative]
    if (original["destination"], original["sha256"]) != (relative, ingest.sha(ledger)):
        raise ValueError("Original observation ledger differs from the reviewed file")
    observation = campaign.read_obs(ledger)
    product, = [p for p in observation["products"] if p["input"]["sha256"] == PRODUCT and p["complete"]]
    event, = [t for t in product["search"]["triggers"] if t["trigger_key"] == TRIGGER and t["instrument_pass"]]
    if event["obsid"] != observation["obsid"] or event["obsid"] != OBSID:
        raise ValueError("Source observation identity differs")
    if event["diagnostics"]["files"]["diagnostics.png"] != "032c0c6469ca3780f603158cb047b451988489df2d54fa1f097537c6ae16bc17":
        raise ValueError("Source diagnostics differ from the reviewed nine-panel image")
    fields = ("obsid", "trigger_key", "ra_deg", "dec_deg", "event_mjd", "position_error_arcsec",
              "trial_duration_s", "duration_s")
    lead = {key:event[key] for key in fields}
    lead.update(input_sha256=PRODUCT, input_name=product["input"]["filename"],
                input_url=product["input"]["url"])
    return product["input"], event, lead


def run():
    if (os.getenv("GITHUB_ACTIONS") != "true" or os.getenv("GITHUB_REF") != CAMPAIGN_REF
            or os.getenv("GITHUB_REPOSITORY") != ingest.REPO):
        raise RuntimeError("Use the Chandra research workflow")
    root = campaign.EVIDENCE
    out = root / "field_followup" / str(OBSID) / f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    out.mkdir(parents=True, exist_ok=True)
    source, event, lead = select(root)
    dump(out/"leads.json", [lead])
    dump(out/"event_input.json", [source])
    optical = root/"lead_optical"/"observations"/str(OBSID)/PRODUCT/(TRIGGER+".json.gz")
    csc = root/"csc_variability"/"observations"/f"{OBSID}.json.gz"
    record = {"primary_run":ingest.RUN, "source_artifact_id":ARTIFACT, "source_zip_digest":ZIP_SHA,
              "source_workflow_commit":ingest.CODE, "followup_commit":os.environ["GITHUB_SHA"],
              "selection":"Exploratory off-target ten-photon concentration in a field with repeated visits.",
              "blind_statistics_changed":False, "stages":{}, "inputs":{}}
    for path in (optical, csc):
        record["inputs"][path.relative_to(root).as_posix()] = (
            {"status":"available", "sha256":ingest.sha(path)} if path.is_file() else {"status":"unavailable"})
    dump(out/"execution.json", record)

    def photons():
        with TemporaryDirectory(dir=os.environ["RUNNER_TEMP"]) as tmp:
            replay(out/"event_input.json", Path(tmp))
            priority.comparison(out, event, Path(tmp)/source["filename"])

    def circulars():
        # Use the archived GCN snapshot already deposited by the two-flare follow-up.
        saved = root/"priority"/str(ingest.RUN)/"alerts/source"
        with TemporaryDirectory(dir=os.environ["RUNNER_TEMP"]) as work:
            folder = Path(work)
            shutil.copyfile(out/"leads.json", folder/"leads.json")
            (folder/"alerts").mkdir()
            (folder/"alerts/source").symlink_to(saved, target_is_directory=True)
            try:
                alerts.main(Namespace(evidence=folder, replay=True, self_test=False))
            finally:
                result = folder/"alerts/review.json"
                if result.is_file():
                    shutil.copyfile(result, out/"gcn_review.json")
                    field_replay.copy_members(folder/"alerts", out, campaign.read(result))
                record["gcn_source_directory"] = saved.relative_to(root).as_posix()

    def visits():
        target_history.run([lead], out/"history", "2026-11-21")
        history = campaign.read(out/"history/history.json")
        target, = history["targets"].values()
        record["visit_coverage"] = {"pointing_query_status":target["pointing_query_status"],
                                    "nearby_obsids":len(target["nearby_obsids"]),
                                    "supported_products":len(target["supported"]),
                                    "coverage_gaps":len(target["gaps"])}
        if target["pointing_query_status"] != "ok":
            raise RuntimeError("Independent pointing query failed; its raw response is retained")

    for name, operation in (("photons", photons), ("gcn", circulars), ("all_years_visits", visits)):
        try:
            operation()
            record["stages"][name] = {"status":"finished"}
        except Exception as error:
            record["stages"][name] = {"status":"failed", "error":f"{type(error).__name__}: {error}"}
        dump(out/"execution.json", record)
    if any(row["status"] == "failed" for row in record["stages"].values()):
        raise RuntimeError("One or more fixed-position checks failed; see execution.json")


if __name__ == "__main__":
    cli=ArgumentParser(description=__doc__)
    cli.add_argument("--backfill",metavar="ORIGINAL_RUN",help="Restore notices and summarize the pinned field products offline")
    args=cli.parse_args()
    if args.backfill:
        if (os.getenv("GITHUB_ACTIONS")!="true" or os.getenv("GITHUB_REF")!=CAMPAIGN_REF or
                os.getenv("GITHUB_REPOSITORY")!=ingest.REPO or
                os.getenv("GITHUB_WORKFLOW_REF","").split("@")[0]!=ingest.REPO+"/"+ingest.FLOW):
            raise RuntimeError("Use the Chandra research workflow")
        field_replay.backfill(campaign.EVIDENCE,args.backfill)
    else:
        run()
