"""Independent SPHEREx DR1 referee audit, using frozen files (and optional archive refetch).

Run from the result directory: python referee_review.py offline
                               python referee_review.py fresh-cutouts dr1|qr3
The offline mode checks independently calculated counts/claims, hashes and
individual receipt consistency, then regenerates analysis products and compares
them to the deposited versions. The fresh-cutouts modes re-fetch original bytes
and verify all recorded SHA-256 values, independently of visit photometry code.
"""
import csv
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from urllib.request import Request, urlopen
from concurrent.futures import ThreadPoolExecutor, as_completed

ROOT = Path(__file__).resolve().parent
HEX64 = re.compile(r"^[0-9a-f]{64}$")
RECEIPTS = ("nep", "wide", "wide_remaining")


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def digest(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def load(path):
    return json.loads(path.read_text())


def require(flag, explanation):
    if not flag:
        raise AssertionError(explanation)


def approx(value, target, tolerance, name):
    require(abs(float(value)-target) <= tolerance,
            f"{name}: found {value}, expected {target} within {tolerance}")


def record_inputs():
    manifest = rows(ROOT / "input_hashes.csv")
    require(len(manifest) == 12, f"Input manifest has {len(manifest)} lines")
    checked, cube_receipts = 0, 0
    for entry in manifest:
        field = ROOT if entry["field"] == "nep" else ROOT / "wide"
        if entry["name"] == "cube":
            receipt = load(field / "receipts/cube.json")
            require(receipt["sha256"] == entry["sha256"], "Cube SHA receipt mismatch")
            require(int(receipt["size"]) == int(entry["bytes"]), "Cube length receipt mismatch")
            require(len(receipt["ranges"]) > 30, "Insufficient byte-range receipts")
            offset = 0
            for segment in receipt["ranges"]:
                require(segment["first"] == offset, "Cube ranges have a gap")
                require(segment["last"] >= segment["first"], "Empty range")
                require(HEX64.fullmatch(segment["sha256"]), "Missing range hash")
                offset = segment["last"] + 1
            require(offset == int(entry["bytes"]), "Cube ranges do not cover FITS")
            cube_receipts += 1
        else:
            item = field / entry["file"]
            require(item.stat().st_size == int(entry["bytes"]), f"Size: {item}")
            require(digest(item) == entry["sha256"], f"SHA-256: {item}")
            checked += 1
    require((checked, cube_receipts) == (10, 2), "Input verification incomplete")
    return dict(verified_catalog_bytes=checked, cube_metadata_range_receipts=cube_receipts,
                caveat="Cube FITS bytes need the separate fresh cube job")


def check_spectra():
    stats = {}
    total_injections = 0
    for fld, relative, run, expected in (
        ("wide", "wide", "full", (59, 59, 59, 52, 8944)),
        ("nep", ".", "nep_two_band", (148, 54, 48, 5, 2208)),
    ):
        folder = ROOT / relative
        analysis = ROOT / ("wide/full" if fld == "wide" else run)
        src = rows(folder / "sources.csv")
        spectrum = rows(folder / "spectra.csv")
        ranks = rows(analysis / "ranks.csv")
        trials = rows(analysis / "injections.csv")
        selected = {x["source_id"] for x in src if x["disposition"] == "extracted"}
        ranked = {x["source_id"] for x in ranks}
        hits = [x for x in ranks if float(x["score"]) >= 5]
        calc = (len(src), len(selected), len(ranked), len(hits), len(trials))
        require(calc == expected, f"{fld} counts {calc}, expected {expected}")
        require(len(spectrum) == len(selected)*102, f"{fld} lacks source-channel rows")
        require({x["source_id"] for x in spectrum} == selected, "Unselected spectrum")
        require(len({(x["source_id"], x["channel"]) for x in spectrum}) == len(spectrum),
                "Duplicate source/channel")
        require(len({x["source_id"] for x in ranks}) == len(ranks), "Duplicate ranks")
        am = load(analysis / "analysis.json")
        require(am["n_injections"] == len(trials), "Injection file length / analysis mismatch")
        model_hashes = am["input_sha256"]
        for key, src_path in (
            ("spectra.csv", folder / "spectra.csv"),
            ("relative_fnu.csv", ROOT / "templates/relative_fnu.csv"),
            ("analyze.py", ROOT / "analyze.py"),
        ):
            require(digest(src_path) == model_hashes[key], f"Wrong model input {key}")
        if fld == "wide":
            require(digest(folder / "catalogs/model_redshifts.csv") ==
                    model_hashes["model_redshifts.csv"], "Model redshift hash mismatch")
        grouping = {}
        for trial in trials:
            frac = trial["amplitude_frac"]
            n, found = grouping.get(frac, (0, 0))
            grouping[frac] = (n+1, found+int(trial["recovered"]))
        expected_grid = (
            {"0.1": (2236,190), "0.25": (2236,961), "0.5": (2236,1665), "1": (2236,1913)}
            if fld == "wide" else
            {"0.1": (552,82), "0.25": (552,437), "0.5": (552,517), "1": (552,525)}
        )
        require(grouping == expected_grid, f"{fld} recovery mismatch: {grouping}")
        stats[fld] = {"catalog_positions": len(src), "extracted": len(selected),
                      "ranked": len(ranked), "raw_score_ge5": len(hits),
                      "recovery": grouping}
        total_injections += len(trials)
    require(total_injections == 11152, "Incorrect total trial denominator")
    return {"fields": stats, "positions_extracted": 113, "ranked": 107,
            "raw_triggers": 57, "injection_trials": total_injections}


def collected_receipts():
    receipts = []
    pairs = set()
    for directory in RECEIPTS:
        folder = ROOT / "visits" / directory
        for target in rows(folder / "targets.csv"):
            source = target["source_id"]
            channel = int(target["channel"])
            pair = (source, channel)
            require(pair not in pairs, f"Duplicate follow-up pair {pair}")
            pairs.add(pair)
            stem = source + ("_" + target["label"] if target["label"] else "")
            recs = rows(folder / (stem + ".csv"))
            summ = load(folder / (stem + ".json"))
            require(summ["source_id"] == source, "Wrong source in visit summary")
            require(len(recs) == summ["attempted"], f"Incomplete ledger {stem}")
            for rec in recs:
                require(rec["release"] in ("dr1","qr3"), "Unknown release")
                require(HEX64.fullmatch(rec["sha256"]), "Empty or malformed cutout hash")
                require(rec["url"].startswith("https://irsa.ipac.caltech.edu/ibe/data/spherex/"),
                        "Unapproved cutout host")
                require("/" + rec["release"] + "/" in rec["url"], "Wrong release URL")
                require(rec["status"] in ("ok", "few_blanks", "aperture_flags", "off_bin", "edge")
                        or rec["status"].startswith(("read_error", "archive_bytes_changed")),
                        f"Unknown visit outcome {rec['status']}")
                rec["_source"] = source
                rec["_channel"] = channel
            receipts.extend(recs)
    return pairs, receipts


def check_receipts():
    pairs, receipts = collected_receipts()
    require(len(pairs) == 63, "Expected 57 original plus six extra feature locations")
    require(len(receipts) == 923, f"Only {len(receipts)} native cutout receipts")
    require(sum(x["status"] == "ok" for x in receipts) == 542, "Wrong clean cutout count")
    require(not [x for x in receipts if x["status"].startswith(("read_error", "archive_bytes_changed"))],
            "Unresolved cutout read error")
    ids = {(x["release"], x["obs_id"]) for x in receipts}
    require(len({x["obs_id"] for x in receipts}) == 339,
            "Distinct observation-ID count differs from report")
    dr1 = {x["obs_id"] for x in receipts if x["release"] == "dr1"}
    qr3 = {x["obs_id"] for x in receipts if x["release"] == "qr3"}
    require(not dr1.intersection(qr3), "Observation IDs shared across releases")
    f = rows(ROOT / "followup_by_feature.csv")
    require(len(f) == 63, "Follow-up index has missing feature rows")
    original = [x for x in f if x["original_trigger"] == "True"]
    require(len(original) == 57, "Lost an original trigger")
    require(len([x for x in original if x["qr3_full_bracket"] == "True"]) == 29,
            "Incorrect independent three-bin follow-up count")
    require(len([x for x in original if x["qr3_full_bracket"] != "True"]) == 28,
            "Incorrect pending count")
    require(not [x for x in original if x["qr3_same_sign_ge3"] == "True"],
            "A signed template residual repeated")
    require({(x["source_id"], int(x["channel"])) for x in f} == pairs,
            "Follow-up index and individual target pairs disagree")
    raw_wide = rows(ROOT / "visits/wide_remaining/all_raw_trigger_checks.csv")
    require(len(raw_wide) == 104, "Wide original trigger ledger missing release rows")
    require(sum(x["three_bins"] == "1" and x["release"] == "qr3" for x in raw_wide) == 25,
            "Incorrect wide QR3 clean comparison count")
    require(sum(x["three_bins"] == "1" and x["release"] == "dr1" for x in raw_wide) == 42,
            "Incorrect wide DR1 complete count")
    require(load(ROOT / "visits/wide_remaining/all_raw_trigger_summary.json")["both_complete"] == 21,
            "Wrong wide matched-release count")
    counts = load(ROOT / "search_counts.json")
    require(counts["wide"]["target_exposure_tests"] +
            counts["nep"]["target_exposure_tests"] == 923,
            "Search-count table differs from individual exposure sum")
    return {"feature_locations": len(pairs), "original_triggers": 57,
            "independently_measured": 29, "pending": 28,
            "native_receipts": len(receipts), "accepted": 542,
            "distinct_observation_ids": len(dr1|qr3),
            "qr3_same_sign_template_repeats": 0}


def check_two_sources():
    wide = ROOT / "wide"
    name = "J135745.42+295804.4"
    w = [x for x in rows(wide / "full/ranks.csv") if x["source_id"] == name]
    require(len(w)==1 and w[0]["template_rank"]=="1", "Leading spike rank changed")
    approx(w[0]["peak_wave_um"],4.6,.0002,"Spike wavelength")
    require(w[0]["peak_channel"] == "91", "Spike channel")
    spec = {int(x["channel"]): x for x in rows(wide/"spectra.csv") if x["source_id"]==name}
    approx(spec[91]["flux"],16.44,.04,"4.6 micron central coadd flux")
    for ch, expected in ((90,6.72),(92,6.41)):
        approx(spec[ch]["flux"], expected, .04, f"flank {ch}")
    rows_19 = [x for x in rows(wide/"coadd_checks.csv") if x["source_id"] == name
               and x["mask"] == "strict" and x["radius_arcsec"] == "19"
               and x["sky_inner"] == "36"]
    require(len(rows_19) == 1, "Missing fixed 19arcsec independent coadd check")
    approx(rows_19[0]["fraction"],1.61,.03,"161% coadd local excess")
    spatial = [x for x in rows(wide/"spatial_localization.csv")
               if x["source_id"] == name and x["mask"] == "strict"]
    require(len(spatial)==1, "Missing spatial localization")
    approx(spatial[0]["single_pixel_fraction"],.87,.005,"Single pixel fraction")
    approx(spatial[0]["excess_to_cont_centroid_arcsec"],9.56,.05,"Centroid offset")
    require(spatial[0]["single_pixel_nhits"]=="14", "Spike NHITS 14")
    visit=load(ROOT/"visits/wide"/(name+".json"))
    require(visit["attempted"]==13, "Leading spike visit count")
    approx(visit["groups"]["qr3_all"]["fractional_contrast"],-.032,.001,
           "Leading spike held-out QR3")
    galaxy=load(ROOT/"visits/profile_fit.json")
    require(galaxy["source_id"]=="J135607.93+300452.8", "Wrong PAH source")
    for rel, expected, n in (("dr1",3.288,11),("qr3",3.292,7)):
        info=galaxy[rel]
        approx(info["rest_center_um"],expected,.001,"UGC rest center "+rel)
        require(info["accepted_in_window"]==n, "UGC native profile coverage")
    transfer=load(ROOT/"visits/wide/pah_transfer.json")
    require((transfer["training_n"],transfer["transfer_n"])==(11,7),
            "UGC train/holdout exposure totals")
    approx(transfer["linear_plus_pah"]["qr3_no_refit_chi2"],21.29,.03,"PAH QR3 chi-square")
    approx(transfer["linear"]["qr3_no_refit_chi2"],241.05,.03,"Continuum QR3 chi-square")
    return {"spike_um":4.6,"coadd_fraction":float(rows_19[0]["fraction"]),
            "spike_qr3_fraction":visit["groups"]["qr3_all"]["fractional_contrast"],
            "spatial_fraction":float(spatial[0]["single_pixel_fraction"]),
            "pah_rest_um":[galaxy["dr1"]["rest_center_um"],galaxy["qr3"]["rest_center_um"]],
            "pah_qr3_chisq":transfer["linear_plus_pah"]["qr3_no_refit_chi2"]}


def replay():
    expected = {}
    dirs = ("nep_two_band","nep_full","wide/full")
    outputs=("analysis.json","ranks.csv","injections.csv","recovery_by_shape.csv",
             "eligibility.csv","filter_scores.csv","residuals.csv")
    for folder in dirs:
        for file in outputs:
            p = ROOT/folder/file
            if p.is_file(): expected[str(p.relative_to(ROOT))]=digest(p)
    extra=("wide/full/peer_channel_check.csv","nep_two_band/peer_channel_check.csv",
           "wide/full/redshift_check.csv","search_counts.json","followup_by_feature.csv",
           "visits/wide_remaining/all_raw_trigger_summary.json")
    for file in extra:
        p=ROOT/file
        if p.is_file():expected[file]=digest(p)
    env=dict(os.environ,OPENBLAS_NUM_THREADS="1",OMP_NUM_THREADS="1")
    cmds=[
        [sys.executable,"analyze.py","--out","nep_full"],
        [sys.executable,"analyze.py","--two-band","--out","nep_two_band"],
        [sys.executable,"analyze.py","--input","wide/spectra.csv","--redshifts",
         "wide/catalogs/model_redshifts.csv","--out","wide/full"],
        [sys.executable,"channel_check.py","nep_two_band"],
        [sys.executable,"channel_check.py","wide/full"],
        [sys.executable,"redshift_check.py","wide"],
        [sys.executable,"visits/summarize_wide_triggers.py"],
        [sys.executable,"evidence.py"],
        [sys.executable,"visits/profile_fit.py"],
        [sys.executable,"visits/pah_transfer.py","."],
    ]
    for cmd in cmds:
        print("REPLAY", " ".join(cmd), flush=True)
        subprocess.run(cmd,cwd=ROOT,env=env,check=True,timeout=2300)
    checked, drift = 0, []
    for path, old in expected.items():
        actual=digest(ROOT/path)
        checked += 1
        if actual != old: drift.append(path)
    require(not drift, f"Regenerated output bytes differ: {drift}")
    return {"replay_commands":len(cmds),"byte_identical_outputs":checked}


def fresh_cutouts(release):
    """Count only independently downloaded, SHA-identical cutout bytes as verified.

    IRSA's IBE occasionally responds 503. Reattempt *only* failures, with
    decreasing concurrency, and publish a receipt audit even on an outage.
    """
    require(release in ("dr1", "qr3"), "Choose dr1 or qr3")
    _, receipts = collected_receipts()
    subset=[x for x in receipts if x["release"] == release]
    def check_one(row):
        last=None
        for attempt in range(4):
            try:
                request=Request(row["url"],headers={"User-Agent":"Open-Space-independent-referee/1.0"})
                with urlopen(request,timeout=70) as reply:
                    require(reply.status==200, f"Unexpected HTTP {reply.status}")
                    h=hashlib.sha256()
                    n=0
                    while True:
                        chunk=reply.read(1024*1024)
                        if not chunk:break
                        h.update(chunk);n+=len(chunk)
                require(n>0, "Empty cutout")
                require(h.hexdigest()==row["sha256"],
                        f"SHA MISMATCH {row['_source']} {row['url']}")
                return n
            except Exception as ex:
                last=ex
                if attempt != 3: time.sleep(2.0*(2**attempt))
        raise RuntimeError(f"{row['_source']} {row['url']}: {last}")

    pending=list(subset)
    successes=0
    size=0
    errors=[]
    for iteration,workers in enumerate((4,2,1),start=1):
        if not pending: break
        if iteration>1:
            time.sleep(10*iteration)
        failed=[]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures={pool.submit(check_one,row):row for row in pending}
            for n,future in enumerate(as_completed(futures),start=1):
                row=futures[future]
                try:
                    size += future.result()
                    successes+=1
                except Exception as e:
                    failed.append((row,str(e)))
                if n%100==0 or n==len(pending):
                    print("verified",release,"pass",iteration,n,"of",len(pending),
                          "remaining",len(failed),flush=True)
        errors=[reason for _,reason in failed]
        pending=[row for row,_ in failed]
        if pending: print("RETRY",release,"unavailable",len(pending),"pass",iteration,flush=True)
    report={"release":release,"receipts":len(subset),"fresh_hashed_cutouts":successes,
            "unverified_cutouts":len(pending),"total_fetched_bytes":size,
            "failures":errors}
    (ROOT/("referee_cutouts_"+release+".json")).write_text(json.dumps(report,indent=2)+"\n")
    require(not errors, f"{len(errors)} cutouts could not be freshly verified: {errors[:6]}")
    return report


def main():
    mode = sys.argv[1] if len(sys.argv)>1 else "offline"
    if mode=="offline":
        report={"input_hashes":record_inputs(),"sample":check_spectra(),
                "native_receipts":check_receipts(),"objects":check_two_sources()}
        report["replay"]=replay()
        # Run independent checks again, against regenerated files.
        check_spectra(); check_receipts(); check_two_sources()
        (ROOT/"referee_report.json").write_text(json.dumps(report,indent=2)+"\n")
    elif mode=="fresh-cutouts":
        report=fresh_cutouts(sys.argv[2])
        (ROOT/("referee_cutouts_"+sys.argv[2]+".json")).write_text(json.dumps(report,indent=2)+"\n")
    else:
        raise SystemExit("Use offline or fresh-cutouts dr1|qr3")
    print("REFEREE RESULT:",json.dumps(report,sort_keys=True),flush=True)


if __name__=="__main__":
    main()
