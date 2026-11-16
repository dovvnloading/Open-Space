"""Run the frozen disjoint O4b campaign, verify per-block checkpoints, and package ledgers.

OPENBLAS_NUM_THREADS=1 python holdout.py RAW_CACHE CHECKPOINTS --workers 5
Monthly ZIPs contain one folder per GPS epoch. Foreground CSV includes every
live second, including score zero. Sparse rows contain every shifted or
unshifted one-second maximum >=3. Masks are half-open runs in block seconds.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
from dataclasses import asdict
import hashlib
import io
import json
import subprocess
from pathlib import Path
import time
from zipfile import ZipFile, ZIP_DEFLATED, ZipInfo

import numpy as np

import acquire
import core
import holdout_data
import holdout_sensitivity
import survey
import vetting

HERE = Path(__file__).resolve().parent
CAPTURE = 3.0
SEED_FILES = ("core.py", "holdout_data.py", "holdout_sensitivity.py", "sensitivity.py",
              "acquire.py", "survey.py", "vetting.py", "public_events.csv")


def csv_bytes(header, rows):
    buf = io.StringIO(newline="")
    writer = csv.writer(buf, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buf.getvalue().encode()


def write_member(archive, name, data):
    entry = ZipInfo(name)
    entry.compress_type = ZIP_DEFLATED
    entry.external_attr = 0o100644 << 16
    archive.writestr(entry, data)


def good_runs(good):
    points = np.flatnonzero(np.diff(np.r_[0, good.astype("i1"), 0]))
    return points.reshape(-1, 2).tolist()


def trigger_row(gps, shift, tr):
    return [shift, gps+tr.time, gps+tr.time+shift, tr.score, *tr.band, tr.duration,
            tr.lag_seconds, tr.coherence, tr.sign]


TRIGGER_COLS = ["shift_s", "H1_peak_gps", "L1_tile_gps", "score", "low_Hz", "high_Hz",
                "duration_s", "residual_lag_s", "coherence", "sign"]


def calculate(cache, block, draws):
    gps = block["gps"]
    inputs, prepared = [], []
    masks = {}
    for detector in ("H1", "L1"):
        raw, good, counts = survey.load(cache, block, detector, excise=True)
        p = core.prepare(raw, acquire.FS, good)
        inputs.append((raw, good))
        prepared.append(p)
        masks[detector] = dict(cumulative_seconds=counts, post_veto_runs=good_runs(good),
                               searchable_runs=good_runs(p.live))
    live = prepared[0].live & prepared[1].live
    report = {"gps": gps, "month": block["month"], "live_s": int(live.sum()), "masks": masks,
              "excisions": vetting.matches(gps-vetting.MARGIN, gps+acquire.DURATION+vetting.MARGIN),
              "slide_columns": survey.EPOCH_FIELDS, "slides": []}
    baseline = np.zeros(len(live))
    tails, fg_rows = [], []
    for shift in survey.SHIFTS:
        slide = core.scan(*prepared, shifts=(shift,))[shift]
        top = max(slide.triggers, key=lambda tr: tr.score, default=None)
        ge20 = 0
        for tr in slide.triggers:
            if tr.score >= CAPTURE:
                tails.append(trigger_row(gps, shift, tr))
            if tr.score >= survey.THRESHOLD:
                ge20 += 1
        report["slides"].append([shift, slide.live_seconds, slide.nominal_tiles, slide.gated_tiles,
                                 len(slide.triggers), ge20, top.score if top else 0.])
        if shift:
            continue
        if slide.live_seconds != live.sum():
            raise ValueError("Foreground exposure mismatch")
        report["top_foreground"] = None if top is None else {**asdict(top), "gps": gps+top.time}
        by_sec = {int(tr.time): tr for tr in slide.triggers}
        for sec in np.flatnonzero(live):
            tr = by_sec.get(int(sec))
            if tr is not None:
                baseline[sec] = tr.score
            fg_rows.append([int(sec), float(baseline[sec])])
    injections = holdout_sensitivity.run(block, inputs, prepared, live, baseline, draws)
    return report, fg_rows, tails, injections



def verified_pair(cache, block):
    try:
        return holdout_data.acquire_holdout(cache, block)
    except (OSError, ValueError, subprocess.CalledProcessError):
        # OSDF serves the same immutable release; always use the pinned MD5.
        mirror = {**block, "sources": {}}
        for det, source in block["sources"].items():
            directory, name = source["url"].rsplit("/", 2)[1:]
            url = (f"https://origin.ligo.caltech.edu:8443/gwdata/O4b/{acquire.RELEASE}/"
                   f"STRAIN_HDF/{det}/{directory}/{name}")
            mirror["sources"][det] = {**source, "url": url}
        return holdout_data.acquire_holdout(cache, mirror)


def execute(cache, checkpoints, block, draws, fingerprint, retain):
    gps = block["gps"]
    path = checkpoints / f"{gps}.zip"
    if path.exists():
        with ZipFile(path) as z:
            if z.testzip() is not None or json.loads(z.read("epoch.json"))["configuration_sha256"] != fingerprint:
                raise ValueError(f"Invalid/incompatible checkpoint {path}")
        return gps, 0.
    start = time.monotonic()
    files = verified_pair(cache, block)
    report, fg_rows, tails, injections = calculate(cache, block, draws)
    report["configuration_sha256"] = fingerprint
    tmp = path.with_suffix(".part")
    with ZipFile(tmp, "w") as z:
        write_member(z, "epoch.json", json.dumps(report, allow_nan=False, separators=(",", ":")))
        write_member(z, "foreground.csv", csv_bytes(["offset_second", "score"], fg_rows))
        write_member(z, "tail.csv", csv_bytes(TRIGGER_COLS, tails))
        write_member(z, "injections.json", json.dumps(injections, separators=(",", ":"), allow_nan=False))
    with ZipFile(tmp) as z:
        if z.testzip() is not None:
            raise OSError(f"Bad checkpoint CRC {gps}")
    tmp.replace(path)
    if not retain:
        for f in files.values():
            f.unlink(missing_ok=True)
    return gps, time.monotonic()-start



def monthly_shard(month, starts, checkpoints, out):
    dest = out / f"ledgers-{month}.zip"
    temp = dest.with_suffix(".part")
    with ZipFile(temp, "w") as merged:
        for gps in starts:
            with ZipFile(checkpoints / f"{gps}.zip") as frame:
                for member in sorted(frame.namelist()):
                    write_member(merged, f"{gps}/{member}", frame.read(member))
    temp.replace(dest)
    return {"file": dest.name, "bytes": dest.stat().st_size,
            "sha256": hashlib.sha256(dest.read_bytes()).hexdigest()}

def campaign(blocks, checkpoints, out, fingerprint, selection_digest):
    month_lists = defaultdict(list)
    rows, top, trials, crossings = [], [], [], []
    sparse_scores = []
    monitors = []
    for b in blocks:
        gps = b["gps"]
        month_lists[b["month"]].append(gps)
        with ZipFile(checkpoints / f"{gps}.zip") as z:
            rec = json.loads(z.read("epoch.json"))
            if rec["configuration_sha256"] != fingerprint or rec["gps"] != gps or z.testzip() is not None:
                raise ValueError(f"Bad epoch checkpoint for {gps}")
            fg = list(csv.DictReader(io.StringIO(z.read("foreground.csv").decode())))
            if len(fg) != rec["live_s"] or len(set(t["offset_second"] for t in fg)) != len(fg):
                raise ValueError("Foreground ledger exposure mismatch")
            if max((float(t["score"]) for t in fg), default=0.) != rec["slides"][0][-1]:
                raise ValueError("Foreground maximum mismatch")
            rows.append(rec["slides"])
            if rec["top_foreground"] is not None:
                top.append(rec["top_foreground"])
            trials.extend(json.loads(z.read("injections.json")))
            monitors.append({"gps": gps, "live_s": rec["live_s"], "foreground": rec["top_foreground"],
                             "foreground_ge20": rec["slides"][0][-2]})
            for t in csv.DictReader(io.StringIO(z.read("tail.csv").decode())):
                s, score = int(t["shift_s"]), float(t["score"])
                if s:
                    sparse_scores.append(score)
                if score >= survey.THRESHOLD:
                    crossings.append({k: int(v) if k in ("shift_s", "sign") else float(v) for k,v in t.items()})
    ranks = sorted(top, key=lambda tr: -tr["score"])
    matrices = [np.array(r, dtype=np.float64) for r in rows]
    all_slides = []
    for idx, shift in enumerate(survey.SHIFTS):
        parts = [m[idx] for m in matrices]
        if any(p[0] != shift for p in parts):
            raise ValueError("Shift mismatch")
        all_slides.append([shift, *(int(sum(p[j] for p in parts)) for j in range(1, 6)),
                           float(max(p[-1] for p in parts))])
    fg, bg = all_slides[0], all_slides[1:]
    louder = sum(s[-1] >= fg[-1] for s in bg)
    nbg = sum(s[-2] for s in bg)
    if len(crossings) != nbg + fg[-2]:
        raise ValueError("Threshold counts disagree with sparse ledger")
    summary = {"configuration_sha256": fingerprint, "selection_sha256": selection_digest,
               "n_frames": len(blocks), "scheduled_s": acquire.DURATION*len(blocks),
               "HL_live_s": fg[1], "background_s": sum(s[1] for s in bg),
               "foreground_ge20": fg[-2], "background_ge20": nbg,
               "foreground_loudest": fg[-1], "background_campaigns_ge_loudest": louder,
               "campaign_rank_fraction": (louder+1)/(len(bg)+1),
               "background_seconds_ge_loudest": sum(s >= fg[-1] for s in sparse_scores) if fg[-1] >= CAPTURE else None,
               "top_epochs": ranks[:10], "threshold_crossings": crossings,
               "slide_columns": survey.EPOCH_FIELDS, "campaign_slides": all_slides, "epochs": monitors}
    recovered = sum(t["recovered"] for t in trials)
    if len(trials) != len(set(t["trial"] for t in trials)):
        raise ValueError("Duplicate injection trials")
    summary["rate"] = holdout_sensitivity.summarize(summary["scheduled_s"], fg[-2], recovered, len(trials))
    shard_receipts = []
    for month, starts in sorted(month_lists.items()):
        shard_receipts.append(monthly_shard(month, starts, checkpoints, out))
    summary["monthly_ledgers"] = shard_receipts
    survey.write_json(out / "campaign.json", summary)
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    p.add_argument("checkpoints", type=Path)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--selection", type=Path, default=HERE/"results"/"holdout_selection.json")
    p.add_argument("--out", type=Path, default=HERE/"results"/"holdout")
    p.add_argument("--retain", action="store_true", help="retain raw HDF5 pairs after checkpoint")
    p.add_argument("--aggregate-only", action="store_true")
    p.add_argument("--pack-ready-months", action="store_true", help="write only completed monthly ledger shards")
    args = p.parse_args()
    if args.workers < 1:
        p.error("workers must be positive")
    for folder in (args.cache, args.checkpoints, args.out):
        folder.mkdir(exist_ok=True, parents=True)
    raw = args.selection.read_bytes()
    blocks = json.loads(raw)["blocks"]
    if args.pack_ready_months:
        for month in sorted(set(b["month"] for b in blocks)):
            starts = [b["gps"] for b in blocks if b["month"] == month]
            if all((args.checkpoints / f"{gps}.zip").exists() for gps in starts):
                print(monthly_shard(month, starts, args.checkpoints, args.out), flush=True)
        return
    selection_digest = hashlib.sha256(raw).hexdigest()
    digest = hashlib.sha256(raw)
    for name in SEED_FILES:
        digest.update((HERE/name).read_bytes())
    digest.update(str((survey.SHIFTS, survey.THRESHOLD, CAPTURE)).encode())
    fingerprint = digest.hexdigest()
    draws = holdout_sensitivity.draw_plan(blocks)
    if not args.aggregate_only:
        jobs = [(args.cache,args.checkpoints,b,draws[b["gps"]],fingerprint,args.retain) for b in blocks]
        with ProcessPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(execute, *job):job[2]["gps"] for job in jobs}
            for i, future in enumerate(as_completed(futures), 1):
                gps, dt = future.result()
                print(f"{i}/{len(jobs)} GPS {gps}, {dt:.1f} s", flush=True)
    summary = campaign(blocks, args.checkpoints, args.out, fingerprint, selection_digest)
    print({k: summary[k] for k in ("n_frames","HL_live_s","background_s","foreground_loudest",
                                    "foreground_ge20","background_ge20","campaign_rank_fraction","rate")}, flush=True)


if __name__ == "__main__":
    main()
