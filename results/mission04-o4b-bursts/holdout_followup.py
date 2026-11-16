"""Audit holdout calibration and follow the final loudest / threshold-crossing H/L seconds.

From this directory, run python holdout_followup.py calibration CACHE, then
python holdout_followup.py follow CACHE after results/holdout/campaign.json exists.
"""
import argparse
from dataclasses import asdict
import hashlib
import json
import math
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import h5py
import numpy as np

from acquire import DURATION, FS, DQ_REQUIRED, INJ_REQUIRED, fetch, load_detector, require_flags
from core import confirm
from detector_context import acquire_context, summarize_calibration
from diagnose import auxiliary, strain_data, centers_near, local_power
from holdout_data import acquire_holdout
import vetting

ROOT = Path(__file__).parent
RESULTS = ROOT/"results"
MANIFEST_SHA = "22e00840a0c9fc81fa5f98fbbf966fa5459264472b7fbb8ede3853f251b4b544"


def calibration(cache, blocks):
    sources = acquire_context(cache)
    rows = summarize_calibration(cache,blocks)
    out = {}
    for ifo in ("H1","L1"):
        selected = [r for r in rows if r["detector"] == ifo]
        values = {}
        for band in selected[0]["calibration"]:
            fields = {}
            for key in selected[0]["calibration"][band]:
                pairs = [r["calibration"][band][key] for r in selected]
                fields[key] = [min(p[0] for p in pairs), max(p[1] for p in pairs)]
            values[band] = fields
        epochs = sorted(set(t for r in selected for t in r["epochs_gps"]))
        out[ifo] = {"nearest_hourly_products":len(epochs),"epochs_gps":epochs,
                    "epochs_count_by_block":{str(r["block_gps"]):len(r["epochs_gps"]) for r in selected},
                    "max_distance_to_nearest_s":max(r["largest_epoch_distance_s"] for r in selected),
                    "block_with_largest_epoch_distance_gps":max(selected,key=lambda r:r["largest_epoch_distance_s"])["block_gps"],
                    "blocks_having_150hz_envelope_outside_0p92_to_1p08":[r["block_gps"] for r in selected if
                        r["calibration"]["150-150Hz"]["magnitude_envelope_16_84"][0] < 0.92 or
                        r["calibration"]["150-150Hz"]["magnitude_envelope_16_84"][1] > 1.08],
                    "extrema_across_frequencies_and_all_selected_seconds":values}
    return {"selection_sha256":MANIFEST_SHA,"sources":sources,"n_blocks":len(blocks),"detectors":out,
            "columns":"Median R_true/R_model and pointwise 16th–84th percentile bounds. "
                      "Envelope extrema are descriptive, not a simultaneous probability band."}


def virgo_delays(cache, block, gps, duration, band):
    """Allow terrestrial delay uncertainty; maximize same delay grid off source."""
    if "V1" not in block["files"]:
        return {"available":False,"reason":"no published file"}
    start = block["gps"]
    offsets = np.linspace(-0.028,0.028,7)
    second = int(gps-start)
    with h5py.File(cache/block["files"]["V1"]) as f:
        good = require_flags(f["quality/simple"],"DQmask","DQShortnames",DQ_REQUIRED,start)
        good &= require_flags(f["quality/injections"],"Injmask","InjShortnames",INJ_REQUIRED,start)
        if not good[second]:
            return {"available":False,"reason":"Virgo masks exclude target"}
        left, right = max(0,second-32),min(DURATION,second+33)
        bad_left, bad_right = np.flatnonzero(~good[left:second]), np.flatnonzero(~good[second:right])
        if bad_left.size:
            left += int(bad_left[-1])+1
        if bad_right.size:
            right = second+int(bad_right[0])
        centers = centers_near(gps,start+left,start+right,duration+0.056)
        if centers is None or len(centers)<8:
            return {"available":False,"reason":"insufficient contiguous Virgo data"}
        data = f["strain/Strain"][left*FS:right*FS]
    measured = [local_power(data,FS,start+left,centers+lag,duration,band) for lag in offsets]
    if any(not m["available"] for m in measured):
        return {"available":False,"reason":"nonfinite or degenerate Virgo power"}
    power = np.vstack([m["all_norm_energies"] for m in measured])
    maxima = np.max(power,axis=0)
    return {"available":True,"trial_lags_s":offsets.tolist(),"strongest_lag_s":float(offsets[np.argmax(power[:,0])]),
            "max_normalized_band_energy":float(maxima[0]),"offsource_ge":int(np.sum(maxima[1:]>=maxima[0])),
            "offsource_n":len(maxima)-1}


def follow(cache, block, target):
    gps, band, duration = target["gps"], target["band"], target["duration"]
    # Overlap public strain transfers; each file is verified independently.
    with ThreadPoolExecutor(max_workers=2) as pool:
        hl = pool.submit(acquire_holdout,cache,block)
        if "virgo_source" in block:
            r = block["virgo_source"]
            virgo = pool.submit(fetch,r["url"],cache/r["name"],r["md5"])
            virgo.result()
        paths = hl.result()
    files = dict(block["files"])
    if "virgo_source" in block:
        r = block["virgo_source"]
        files["V1"] = r["name"]
    block = {**block, "files":files}
    info = {"target":target, "survey_block_gps":block["gps"], "catalog_matches_within_32s":[],
            "catalog_receipt_sha256":hashlib.sha256((ROOT/vetting.RECEIPT_FILENAME).read_bytes()).hexdigest()}
    all_events = vetting.events()
    for source in ("GraceDB","LVK","GWOSC"):
        e = min((e for e in all_events if e["source"] == source),key=lambda e:abs(e["gps"]-gps))
        info.setdefault("nearest_published_by_source",{})[source] = {**e,"delta_s":e["gps"]-gps}
    info["catalog_matches_within_32s"] = vetting.matches(gps-32,gps+32)
    info["released_masks_in_winning_window"] = {}
    low = math.floor(gps-duration/2)-block["gps"]
    high = math.ceil(gps+duration/2)-block["gps"]
    for ifo,name in files.items():
        with h5py.File(cache/name) as source:
            flags = {}
            for group,mask,labels,required in (
                ("quality/simple","DQmask","DQShortnames",DQ_REQUIRED),
                ("quality/injections","Injmask","InjShortnames",(*INJ_REQUIRED,"NO_CW_HW_INJ")),
            ):
                for flag in required:
                    v = require_flags(source[group],mask,labels,(flag,),block["gps"])
                    flags[flag] = bool(np.all(v[low:high]))
            info["released_masks_in_winning_window"][ifo] = flags
    orig, good = [], []
    for ifo in ("H1","L1"):
        start, strain, mask = load_detector(paths[ifo])
        orig.append(strain);good.append(vetting.mask_seconds(start,mask))
    try:
        independent = confirm(*orig,FS,gps-block["gps"],band=band,duration=duration,good_h=good[0],good_l=good[1])
        info["independent_time_domain"] = {**asdict(independent), "peak_gps":block["gps"]+independent.peak_time}
    except ValueError as e:
        info["independent_time_domain"] = {"unavailable":str(e)}
    del orig, strain
    info["strain"] = strain_data(cache, block, gps, duration, band)
    if "V1" not in files:
        info["strain"]["V1"] = {"available":False, "reason":"no GWOSC file listed for block"}
    info["virgo_travel_time_check"] = virgo_delays(cache,block,gps,duration,band)
    info["auxiliary"] = auxiliary(cache,gps,duration,band)
    return info


def target_list(campaign):
    peaks = [campaign["top_epochs"][0]]
    for row in campaign["threshold_crossings"]:
        if row["shift_s"] != 0:
            continue
        t = row["H1_peak_gps"]
        if not any(t == old["gps"] for old in peaks):
            peaks.append({"gps":t,"score":row["score"],"band":[row["low_Hz"],row["high_Hz"]],
                          "duration":row["duration_s"], "lag_seconds":row["residual_lag_s"],
                          "coherence":row["coherence"], "sign":row["sign"]})
    return peaks


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("task",choices=("calibration","follow"))
    p.add_argument("cache",type=Path)
    p.add_argument("--results",type=Path,default=RESULTS)
    args = p.parse_args()
    args.cache.mkdir(parents=True,exist_ok=True)
    raw = (args.results/"holdout_selection.json").read_bytes()
    if hashlib.sha256(raw).hexdigest() != MANIFEST_SHA:
        raise ValueError("Wrong holdout selection")
    blocks = json.loads(raw)["blocks"]
    out = args.results/"holdout_followup"
    out.mkdir(parents=True,exist_ok=True)
    if args.task == "calibration":
        objects = [(out/"calibration.json",calibration(args.cache, blocks))]
    else:
        acquire_context(args.cache)
        campaign = json.loads((args.results/"holdout"/"campaign.json").read_text())
        if campaign["selection_sha256"] != MANIFEST_SHA or campaign["n_frames"] != len(blocks):
            raise ValueError("Incomplete or inconsistent campaign summary")
        objects = []
        for t in target_list(campaign):
            block = next(b for b in blocks if b["gps"]<=t["gps"]<b["gps"]+DURATION)
            result = follow(args.cache,block,t)
            result["selection_sha256"] = MANIFEST_SHA
            objects.append((out/f"foreground_{math.floor(t['gps'])}.json",result))
    for path, data in objects:
        path.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
        print(path,flush=True)


if __name__ == "__main__":
    main()
