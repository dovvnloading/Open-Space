"""Compare ACIS-S4 frame occupancy in identical real and scrambled trigger ledgers."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import os
from pathlib import Path

import campaign
import ingest_primary as ingest
import progress
import report


SELECTION = "ACIS CCD 8, instrument pass, at least three selected trial photons with identical EVENTS.TIME"
STREAMS = ["real", *[f"control-{n:02}" for n in range(24)]]


def row_pattern(path, event):
    """Compare CHIPY rows of selected aperture photons sharing one recorded TIME."""
    checksum = (event.get("diagnostics") or {}).get("files",{}).get("event_rows.csv.gz")
    if not checksum or not path.is_file():
        return {"row_audit":"unavailable"}
    ingest.require(ingest.sha(path)==checksum, "Original real-event table hash differs")
    times = defaultdict(list)
    n_peak = 0
    with gzip.open(path,"rt",encoding="ascii",newline="") as f:
        for p in csv.DictReader(f):
            if p["phase"]=="trigger" and p["region"]=="aperture":
                x,y=int(p["chipx"]),int(p["chipy"])
                ingest.require(1<=x<=1024 and 1<=y<=1024,"Invalid ACIS detector coordinate")
                times[float(p["time_met_s"])].append((x,y))
                n_peak+=int(p["in_peak"])
    counts=[len(events) for events in times.values()]
    ingest.require((sum(counts),len(times),max(counts),n_peak)==
                   (event["on_counts"],event["distinct_frames"],event["max_counts_per_frame"],event["peak_counts"]),
                   "Real source-row aggregation disagrees with the frozen trigger")
    maximum = 0
    repeated = 0
    for group in times.values():
        node_row = Counter(((x-1)//256,y) for x,y in group)
        best = max(node_row.values())
        maximum = max(maximum,best)
        repeated += best>=2
    return {"row_audit":"measured","event_rows_sha256":checksum,"max_same_node_row_and_time":maximum,
            "time_tags_with_repeated_node_row":repeated}


def collect_product(root, oid, product, summaries, ledger_sha):
    search = product["search"]
    controls = sorted(search["controls"],key=lambda v:v["replica"])
    ingest.require([v["replica"] for v in controls]==list(range(24)), "Missing or repeated control replica")
    domains = [search["triggers"],*[c["triggers"] for c in controls]]
    rows = []
    for stream,events in zip(STREAMS,domains,strict=True):
        totals=summaries[stream]
        for e in events:
            if e["ccd_id"]!=8:
                continue
            totals["raw_s4"]+=1
            if not e["instrument_pass"]:
                continue
            totals["instrument_pass_s4"]+=1
            if e["max_counts_per_frame"]<3:
                continue
            totals["selected_s4"]+=1
            status=(e.get("vet") or {}).get("status","catalog_incomplete")
            retained=status not in report.FINAL_VETO
            totals["selected_after_catalog"]+=retained
            key=e["trigger_key"]
            inp=product["input"]
            line={"obsid":oid,"stream":stream,"trigger_key":key,"input_sha256":inp["sha256"],
                  "input_name":inp["filename"],"ledger_sha256":ledger_sha,"ra_deg":e["ra_deg"],
                  "dec_deg":e["dec_deg"],"on_counts":e["on_counts"],"peak_counts":e["peak_counts"],
                  "distinct_time_tags":e["distinct_frames"],"max_counts_per_time_tag":e["max_counts_per_frame"],
                  "trial_duration_s":e["trial_duration_s"],"seed_score":e["seed_score"],
                  "catalog_status":status,"survives_catalog":retained}
            if stream=="real":
                directory=root/"photometry"/str(oid)/inp["filename"]/key
                line.update(row_pattern(directory/"event_rows.csv.gz",e))
            else:
                line["row_audit"]="no archived control photon table"
            rows.append(line)
    return rows


def run(root, name):
    out=root/"frame_forensics"/"population"/name
    out.mkdir(parents=True,exist_ok=True)
    inventory=root/"inventory/inventory.json"
    ids=[r["obsid"] for r in campaign.read(inventory)]
    ingest.require(len(ids)==len(set(ids)), "Duplicate dated inventory observation")
    summaries={s:Counter(raw_s4=0,instrument_pass_s4=0,selected_s4=0,selected_after_catalog=0) for s in STREAMS}
    selected=[]
    ledgers=[]
    for oid in sorted(ids):
        path=root/"observations"/f"{oid}.json.gz"
        if not path.is_file():
            continue
        obs=campaign.read_obs(path)
        if not progress.describe(obs)["scanned_all_products"]:
            continue
        sha=ingest.sha(path)
        ingest.require(obs["obsid"]==oid,"Observation and dated inventory differ")
        ledgers.append({"obsid":oid,"sha256":sha})
        for prod in obs["products"]:
            if prod.get("complete"):
                selected.extend(collect_product(root,oid,prod,summaries,sha))
    fields="obsid stream trigger_key input_sha256 input_name ledger_sha256 ra_deg dec_deg on_counts peak_counts distinct_time_tags max_counts_per_time_tag trial_duration_s seed_score catalog_status survives_catalog row_audit event_rows_sha256 max_same_node_row_and_time time_tags_with_repeated_node_row".split()
    report.rows_to_csv(out/"flagged_windows.csv", fields, selected)
    progress.csv_gz(out/"ledgers.csv.gz",["obsid","sha256"],ledgers)
    report.rows_to_csv(out/"streams.csv",["stream","raw_s4","instrument_pass_s4","selected_s4","selected_after_catalog"],
                       ({"stream":s,**summaries[s]} for s in STREAMS))
    cells=[summaries[s]["selected_s4"] for s in STREAMS[1:]]
    missing=sum(r["stream"]=="real" and r["row_audit"]=="unavailable" for r in selected)
    summary={"selection":SELECTION,"time_key":"Identical TIME in the source CCD and original uninterrupted GTI",
             "dated_inventory":len(ids),"complete_scans":len(ledgers),"inventory_sha256":ingest.sha(inventory),
             "stream_counts":summaries,"real_missing_photon_tables":missing,"implementation_sha256":ingest.sha(Path(__file__)),
             "selected_windows_sha256":ingest.sha(out/"flagged_windows.csv"),"ledgers_sha256":ingest.sha(out/"ledgers.csv.gz"),
             "streams_sha256":ingest.sha(out/"streams.csv"),"run_id":os.environ["GITHUB_RUN_ID"],
             "attempt":os.environ["GITHUB_RUN_ATTEMPT"],"checkout":os.environ["GITHUB_SHA"]}
    ingest.save(out/"summary.json",summary)
    (out/"README.md").write_text("\n".join(["# ACIS-S4 frame occupancy", "",
        f"**{summaries['real']['selected_s4']} real instrument-passing S4 windows had at least three photons at one time tag.** "
        f"The 24 scrambled streams had {min(cells)}–{max(cells)} such windows each, totaling {sum(cells)}.",
        "",f"The same count criterion was applied to all 25 streams across {len(ledgers):,} fully scanned dated observations. "
        "Counts here are per event-file window, before cross-product deduplication. Multiple control realizations share "
        "the same images and event populations. Shuffling also disrupts real instrumental time correlations.", "",
        f"For the real windows, the photon tables additionally summarize how many photons shared both an identical time "
        f"and a CHIPY row within a 256-column CHIPX node. {missing} selected real windows lacked a deposited photon table. "
        "That row measurement is available only for the real events. An extended same-row pattern in the same exposure "
        "number can support an S4 readout explanation. A shared timestamp alone has several physical explanations. "
        "The comparisons were chosen after inspection of two events and do not change primary or held-out trial statistics.",
        "", "[Selected windows and their exact inputs](flagged_windows.csv) · [Counts for all 25 streams](streams.csv) · "
        "[Contributing ledgers](ledgers.csv.gz) · [Measurement receipt](summary.json) · "
        "[CXC S4 readout description](https://cxc.cfa.harvard.edu/ciao/why/destreak.html)", ""]))
    return summary


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--final",action="store_true")
    args=p.parse_args()
    ingest.require(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")==ingest.REPO and
                   os.getenv("GITHUB_REF")=="refs/heads/"+ingest.BRANCH, "Use the Chandra research workflow")
    name="final" if args.final else f"{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    run(campaign.EVIDENCE,name)
