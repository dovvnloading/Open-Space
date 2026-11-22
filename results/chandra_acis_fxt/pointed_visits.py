"""All-mission forced photometry for two exploratory high-latitude primary leads."""
from pathlib import Path
import os
from tempfile import TemporaryDirectory

import campaign
import followup
import ingest_primary as ingest
import priority
import target_history
from archive import dump


# Both positions come from the original blind search and were selected for follow-up
# after viewing the discovery images. Neither field selection nor follow-up is a new trial.
POINTS = [
    (101, 14749756165, 22929, "9a02d9ed96de9146bd86ca7a9bbfcb2e420d487e1a49f39b184daefd84fd6fbd",
     "8e2479c63c3771d7470b5a31d92a5b56ecb3325d2ac2ad70248741949983f0aa",
     "Eight-photon window in the Spiderweb field with sparse flanks and weak galaxy associations."),
    (99, 14749108711, 23779, "49987e33e7bc66ca583baa1b7aab2ea59c067b58889807ec78dd851fd752585b",
     "7ba0b13c534f625ba11d2c1cee1fa11eb12d9b452b503ad3cd9f3a6f0f889d88",
     "Nine-photon window in the NSA 103757 field beside a persistently bright sky region."),
]


def original(root, group, artifact, oid, h, image):
    receipt = campaign.read(ingest.receipt_file(root, str(group)))
    members = root/ingest.HISTORY/"members"/f"{group}.json.gz"
    ingest.require(receipt["artifact_id"] == artifact and receipt["code_commit"] == ingest.CODE and
                   receipt["member_manifest_sha256"] == ingest.sha(members), "Unverified primary source")
    path = root/"observations"/f"{oid}.json.gz"
    relative = path.relative_to(root).as_posix()
    member, = [m for m in campaign.read_obs(members) if m["source"] == relative]
    ingest.require(member["destination"] == relative and member["sha256"] == ingest.sha(path), "Changed primary ledger")
    obs = campaign.read_obs(path)
    product, = [p for p in obs["products"] if p.get("complete") and p["input"]["sha256"] == h]
    event, = [t for t in product["search"]["triggers"] if t["trigger_key"] == "real/0"]
    ingest.require(obs["obsid"] == event["obsid"] == oid and event["instrument_pass"] and
                   event["diagnostics"]["files"]["diagnostics.png"] == image, "Changed measured trigger")
    return product["input"], event, receipt, member


def run():
    ingest.require(os.getenv("GITHUB_ACTIONS") == "true" and os.getenv("GITHUB_REPOSITORY") == ingest.REPO and
                   os.getenv("GITHUB_REF") == "refs/heads/"+ingest.BRANCH and
                   os.getenv("GITHUB_WORKFLOW_REF","").split("@",1)[0] == ingest.REPO+"/"+ingest.FLOW,
                   "Use the Chandra research workflow")
    root = campaign.EVIDENCE
    out = root/"pointed_visits"/f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    out.mkdir(parents=True, exist_ok=True)
    leads = []
    attribution = []
    errors = []
    for group, artifact, oid, h, png, why in POINTS:
        p, ev, receipt, member = original(root, group, artifact, oid, h, png)
        lead = {k:ev[k] for k in ("obsid","trigger_key","ra_deg","dec_deg","position_error_arcsec","event_mjd",
                                   "trial_duration_s","duration_s")}
        lead.update(input_sha256=h, input_name=p["filename"], input_url=p["url"])
        leads.append(lead)
        dest = out/str(oid)
        dest.mkdir(exist_ok=True)
        dump(dest/"trigger.json", ev)
        dump(dest/"input.json", [p])
        attribution.append({"obsid":oid, "selection":why, "group":group, "source_artifact":artifact,
                            "source_zip_digest":receipt["zip_digest"], "ledger":member, "image_sha256":png})
        attempts = []
        try:
            with TemporaryDirectory(dir=os.environ["RUNNER_TEMP"]) as temp:
                inp = followup.event_bytes(p, oid, Path(temp), [], attempts)
                priority.comparison(dest,ev,Path(inp["path"]))
        except Exception as exc:
            errors.append({"stage":"photons","obsid":oid,"error":f"{type(exc).__name__}: {exc}"})
        dump(dest/"acquisition.json", attempts)
    dump(out/"leads.json", leads)
    dump(out/"sources.json", {"source_run":ingest.RUN,"source_commit":ingest.CODE,
                             "execution_commit":os.environ["GITHUB_SHA"],"targets":attribution,
                             "selection_role":"Exploratory characterization after primary triggering."})
    try:
        target_history.run(leads, out/"history", "2026-11-21")
        history = campaign.read(out/"history/history.json")
        for k, target in history["targets"].items():
            if target["pointing_query_status"] != "ok":
                errors.append({"stage":"observations","identity":k,"status":target["pointing_query_status"]})
    except Exception as exc:
        errors.append({"stage":"history","error":f"{type(exc).__name__}: {exc}"})
    dump(out/"execution.json", {"source_run":ingest.RUN,"execution_run":os.environ["GITHUB_RUN_ID"],
         "implementation_sha256":ingest.sha(Path(__file__)), "errors":errors,
         "held_out_selection_changed":False})
    if errors:
        raise RuntimeError("One or more field-specific checks failed; see execution.json")


if __name__ == "__main__":
    run()
