"""All-mission forced photometry of selected primary leads; optional time-resolved energy follow-up."""
import argparse
from pathlib import Path
import os
import re
from tempfile import TemporaryDirectory

import campaign
import followup
import ingest_primary as ingest
import priority
import pointed_bands
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

SECOND="offtarget-22558"
COHORTS={"high-latitude":[(*p[:5],"real/0",p[5]) for p in POINTS],SECOND:[
    (113,14750167042,22558,"0dadf392e5987516e910002b7e5caf124893a2d65b3474625ed5a1b9d1828ed6",
     "ce06ccbdfd22f88e216f05c0863cd1f6918399a46cd028d062271f277ed3531c","real/0",
     "Fifteen-photon off-target flare selected for comparison of the trial and both quiet flanks.")
]}


def archived_member(root, members, name):
    [member]=[m for m in members if m["source"]==name]
    dest=followup.safe(root,member["destination"])
    candidates=[dest]
    if member["destination"]==name and (not dest.is_file() or ingest.sha(dest)!=member["sha256"]):
        candidates.extend((root/followup.SCOPED).glob("*/before/"+name))
    for p in candidates:
        if p.is_file() and p.stat().st_size==member["bytes"] and ingest.sha(p)==member["sha256"]:
            return p,member
    raise ValueError("Original member no longer matches: "+name)


def original(root, group, artifact, oid, h, image, tag="real/0"):
    ingest.require(re.fullmatch(r"real/(0|[1-9]\d*)",tag),"Expected an exact real trigger key")
    receipt = campaign.read(ingest.receipt_file(root, str(group)))
    members = root/ingest.HISTORY/"members"/f"{group}.json.gz"
    ingest.require(receipt["artifact_id"] == artifact and receipt["code_commit"] == ingest.CODE and
                   receipt["artifact_name"] == "chandra-evidence-"+str(group) and receipt["plan_commit"]==ingest.PLAN and
                   receipt["member_manifest_sha256"] == ingest.sha(members), "Unverified primary source")
    path,member=archived_member(root,campaign.read_obs(members),f"observations/{oid}.json.gz")
    obs = campaign.read_obs(path)
    product, = [p for p in obs["products"] if p.get("complete") and p["input"]["sha256"] == h]
    event, = [t for t in product["search"]["triggers"] if t["trigger_key"] == tag]
    ingest.require(obs["obsid"] == event["obsid"] == oid and event["instrument_pass"] and
                   event["diagnostics"]["files"]["diagnostics.png"] == image, "Changed measured trigger")
    return product["input"], event, receipt, member


def optical_context(root,dest,lead):
    import lead_optical
    rel=Path("observations")/str(lead["obsid"])/lead["input_sha256"]/(lead["trigger_key"]+".json.gz")
    existing=root/"lead_optical"
    selection=dest/"optical_selection.json"
    dump(selection,[lead])
    if (existing/rel).is_file():
        area=existing
        origin="deposited optical record"
    else:
        area=dest/"optical"
        origin="field-specific archive queries"
        lead_optical.run(root,area,selection,False)
    # The optical replay validates source-response hashes and reproduces the saved rows without HTTP.
    lead_optical.run(root,area,selection,True)
    p=area/rel
    obj=campaign.read_obs(p)
    ingest.require(all(obj["lead"][k]==lead[k] for k in lead_optical.FIELDS),"Optical target identity differs")
    saved={"origin":origin,"path":p.relative_to(root).as_posix(),"sha256":ingest.sha(p),"complete":obj["complete"],
           "source_directory":(area/"source").relative_to(root).as_posix(),
           "catalogs":{k:{"status":v["status"],"rows":len(v["rows"])} for k,v in obj["catalogs"].items()}}
    dump(dest/"optical_context.json",saved)
    ingest.require(obj["complete"],"Photometric catalog evidence remains partial")


def run(cohort="high-latitude"):
    ingest.require(os.getenv("GITHUB_ACTIONS") == "true" and os.getenv("GITHUB_REPOSITORY") == ingest.REPO and
                   os.getenv("GITHUB_REF") == "refs/heads/"+ingest.BRANCH and
                   os.getenv("GITHUB_WORKFLOW_REF","").split("@",1)[0] == ingest.REPO+"/"+ingest.FLOW,
                   "Use the Chandra research workflow")
    root = campaign.EVIDENCE
    out = root/"pointed_visits"/f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    if cohort!= "high-latitude":
        out=out/cohort
    out.mkdir(parents=True, exist_ok=True)
    leads = []
    attribution = []
    errors = []
    trace={"cohort":cohort,"source_run":ingest.RUN,"execution_run":os.environ["GITHUB_RUN_ID"],
           "implementation_sha256":ingest.sha(Path(__file__)),"energy_code_sha256":ingest.sha(Path(pointed_bands.__file__)),
           "targets":{},"errors":errors,"held_out_selection_changed":False,"primary_control_replicas":24}
    dump(out/"execution.json",trace)
    for group, artifact, oid, h, png, tag, why in COHORTS[cohort]:
        tracekey=f"{oid}/{h}/{tag}"
        stages=trace["targets"][tracekey]={"source":"checking original artifact members"}
        dump(out/"execution.json",trace)
        try:
            p, ev, receipt, member = original(root, group, artifact, oid, h, png, tag)
        except Exception as exc:
            stages["source"]="failed"
            errors.append({"stage":"source","identity":tracekey,"error":f"{type(exc).__name__}: {exc}"})
            dump(out/"execution.json",trace)
            continue
        stages["source"]="verified"
        lead = {k:ev[k] for k in ("obsid","trigger_key","ra_deg","dec_deg","position_error_arcsec","event_mjd",
                                   "trial_duration_s","duration_s")}
        lead.update(input_sha256=h, input_name=p["filename"], input_url=p["url"])
        leads.append(lead)
        dest = out/str(oid)
        if cohort!= "high-latitude":
            dest=dest/h/tag
        dest.mkdir(parents=True,exist_ok=True)
        dump(dest/"trigger.json", ev)
        dump(dest/"input.json", [p])
        attribution.append({"obsid":oid, "selection":why, "group":group, "source_artifact":artifact,"trigger_key":tag,
                            "source_zip_digest":receipt["zip_digest"], "ledger":member, "image_sha256":png})
        attempts = []
        if cohort==SECOND:
            try:
                entries=campaign.read_obs(root/ingest.HISTORY/"members"/f"{group}.json.gz")
                prefix=f"photometry/{oid}/{p['filename']}/{tag}/"
                _,image_member=archived_member(root,entries,prefix+"diagnostics.png")
                rowfile,row_member=archived_member(root,entries,prefix+"event_rows.csv.gz")
                ingest.require(image_member["sha256"]==png and ev["diagnostics"]["input_sha256"]==h and
                               row_member["sha256"]==ev["diagnostics"]["files"]["event_rows.csv.gz"],
                               "Primary event or image receipt mismatch")
                if oid==22558 and tag=="real/0":
                    ingest.require(ev["on_counts"]==15,"Fifteen-photon discovery-window identity differs")
                ingest.install(rowfile,dest/"event_rows.csv.gz")
                energy=pointed_bands.profile(dest/"event_rows.csv.gz",ev)
                energy.update(obsid=oid,trigger_key=tag,source_group=group,artifact_id=artifact,source_zip_digest=receipt["zip_digest"],
                              ledger_member=member,primary_rows_member=row_member,primary_image_member=image_member,
                              primary_controls=24,original_heldout_selection_changed=False)
                dump(dest/"energy_bins.json",energy)
                (dest/"ENERGY.md").write_text(pointed_bands.markdown(energy))
                stages["energy_bins"]="complete"
            except Exception as exc:
                stages["energy_bins"]="failed"
                errors.append({"stage":"energy_bins","identity":tracekey,"error":f"{type(exc).__name__}: {exc}"})
            dump(out/"execution.json",trace)
        stages["photon_replay"]="acquiring pinned EVT2"
        dump(out/"execution.json",trace)
        try:
            with TemporaryDirectory(dir=os.environ["RUNNER_TEMP"]) as temp:
                inp = followup.event_bytes(p, oid, Path(temp), [], attempts)
                priority.comparison(dest,ev,Path(inp["path"]))
                stages["photons"]="complete"
                if cohort==SECOND and stages.get("energy_bins")=="complete":
                    checked=pointed_bands.regenerate(Path(inp["path"]),ev,energy["event_rows_sha256"],temp)
                    compare=campaign.read(dest/"photon_check.json")
                    for stage,phase in (("before","before"),("flare","trigger"),("after","after")):
                        e=energy["counts_by_phase"][phase]
                        ingest.require(compare["source_annulus_counts"][stage]==e["annulus"]["total"] and
                                       compare["source_counts"][stage]==e["aperture"]["total"],"Native independent count differs")
                    stages["regenerated_event_rows_sha256"]=checked
                stages["photon_replay"]="complete"
        except Exception as exc:
            stages["photon_replay"]="failed"
            errors.append({"stage":"photons","identity":tracekey,"error":f"{type(exc).__name__}: {exc}"})
        dump(dest/"acquisition.json", attempts)
        dump(out/"execution.json",trace)
        if cohort==SECOND:
            try:
                optical_context(root,dest,lead)
                stages["optical"]="complete"
            except Exception as exc:
                stages["optical"]="partial"
                errors.append({"stage":"optical","identity":tracekey,"error":f"{type(exc).__name__}: {exc}"})
            dump(out/"execution.json",trace)
    dump(out/"leads.json", leads)
    dump(out/"sources.json", {"source_run":ingest.RUN,"source_commit":ingest.CODE,
                             "execution_commit":os.environ["GITHUB_SHA"],"targets":attribution,
                             "selection_role":"Exploratory characterization after primary triggering."})
    trace["history"]="querying the all-mission pointing and CSC archives"
    dump(out/"execution.json",trace)
    try:
        target_history.run(leads, out/"history", "2026-11-21")
        history = campaign.read(out/"history/history.json")
        for k, target in history["targets"].items():
            folder=out/"history"/target["folder"]
            trace["targets"][k]["all_mission_visits"]={"status":target["pointing_query_status"],
                "nearby_obsids":len(target["nearby_obsids"]),"supported_products":len(target["supported"]),"gaps":len(target["gaps"])}
            csc=folder/"csc.json"
            if csc.is_file():
                obj=campaign.read(csc)
                trace["targets"][k]["csc"]={"path":csc.relative_to(root).as_posix(),"sha256":ingest.sha(csc),
                    "statuses":{name:obj.get(name,{}).get("status",obj.get("status","unavailable")) for name in ("master","detections")}}
                if cohort==SECOND and any(s!="ok" for s in trace["targets"][k]["csc"]["statuses"].values()):
                    errors.append({"stage":"csc","identity":k,"statuses":trace["targets"][k]["csc"]["statuses"]})
            if target["pointing_query_status"] != "ok":
                errors.append({"stage":"observations","identity":k,"status":target["pointing_query_status"]})
        trace["history"]="finished; individual coverage and service statuses are in targets"
    except Exception as exc:
        trace["history"]="partial"
        errors.append({"stage":"history","error":f"{type(exc).__name__}: {exc}"})
    dump(out/"execution.json",trace)
    if cohort==SECOND:
        for l in leads:
            key=f"{l['obsid']}/{l['input_sha256']}/{l['trigger_key']}"
            d=out/key
            if not (d/"energy_bins.json").is_file():
                continue
            links=[("Source acquisition and status",out/"execution.json"),
                   ("Aperture and nearby comparison",d/"photon_check.json"),("All-mission independent visits",out/"history/history.json")]
            opt=d/"optical_context.json"
            if opt.is_file():
                info=campaign.read(opt)
                links.extend([("Gaia, Pan-STARRS1 and AllWISE receipt",opt),("Per-catalog rows and positions",root/info["path"])])
            if "csc" in trace["targets"][key]:
                links.append(("CSC master and observation responses",root/trace["targets"][key]["csc"]["path"]))
            body=pointed_bands.markdown(campaign.read(d/"energy_bins.json"))
            (d/"ENERGY.md").write_text(body+"\n"+" · ".join(f"[{label}]({os.path.relpath(path,d)})"
                       for label,path in links if path.is_file())+"\n")
    if errors:
        raise RuntimeError("One or more field-specific checks failed; see execution.json")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort",choices=tuple(COHORTS),default="high-latitude")
    run(parser.parse_args().cohort)
