"""Execute the fixed O4b ten-frame H/L burst survey and its held-out GW250114 control."""
import argparse
import csv
from dataclasses import asdict
import gzip
import io
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.ndimage import binary_erosion

import acquire
import core
import vetting

SHIFTS = (0,) + tuple(s for magnitude in range(4, 257, 4) for s in (magnitude, -magnitude))
THRESHOLD = 20.0
EPOCH_FIELDS = ["shift_s", "live_s", "nominal_tiles", "gated_tiles", "seconds_with_gate", "seconds_ge20", "max_score"]


def write_json(path, obj):
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def load(cache, block, detector, excise):
    path = cache / block["files"][detector]
    gps, strain, good = acquire.load_detector(path)
    with h5py.File(path, "r") as f:
        cumulative = np.ones(acquire.DURATION, dtype=bool)
        counts = {}
        for group, mask, names, flags in (
            ("quality/simple", "DQmask", "DQShortnames", acquire.DQ_REQUIRED),
            ("quality/injections", "Injmask", "InjShortnames", acquire.INJ_REQUIRED),
        ):
            for flag in flags:
                cumulative &= acquire.require_flags(f[group], mask, names, [flag], gps)
                counts[flag] = int(cumulative.sum())
    counts["finite"] = int(good.sum())
    if excise:
        good = vetting.mask_seconds(gps, good)
    counts["after_public_event_mask"] = int(good.sum())
    counts["filter_guards"] = int(binary_erosion(good, structure=np.ones(3)).sum())
    return strain, good, counts


def describe(tr, gps):
    return {**asdict(tr), "gps": gps + tr.time}


def ledger(path, gps, live, triggers):
    columns = ["gps_second", "score", "peak_gps", "low_hz", "high_hz", "duration_s",
               "lag_s", "coherence", "sign"]
    by_second = {int(t.time): t for t in triggers}
    with path.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed:
        with io.TextIOWrapper(compressed, newline="", encoding="ascii") as stream:
            writer = csv.writer(stream)
            writer.writerow(columns)
            for t in np.flatnonzero(live):
                tr = by_second.get(int(t))
                if tr is None:
                    writer.writerow((gps + t, 0.0, "", "", "", "", "", "", ""))
                else:
                    writer.writerow((gps+t, tr.score, gps+tr.time, *tr.band, tr.duration,
                                     tr.lag_seconds, tr.coherence, tr.sign))


def run_epoch(cache, output, block):
    gps = block["gps"]
    control = block["role"] == "control"
    inputs = [load(cache, block, d, excise=not control) for d in ("H1", "L1")]
    pairs = [core.prepare(strain, acquire.FS, good) for strain, good, counts in inputs]
    live = pairs[0].live & pairs[1].live
    masks = {d: i[2] for d, i in zip(("H1", "L1"), inputs)}
    for d, p in zip(("H1", "L1"), pairs):
        masks[d]["psd_and_guards"] = int(p.live.sum())
    vgood = None
    if "V1" in block["files"]:
        vs, vgood, vm = load(cache, block, "V1", excise=not control)
        del vs
        masks["V1"] = vm
    report = {"gps_start": gps, "month": block["month"], "role": block["role"],
              "masks_cumulative_seconds": masks, "HL_live_s": int(live.sum()),
              "V1_after_public_event_within_HL_live": int((vgood & live).sum()) if vgood is not None else 0,
              "excisions": vetting.matches(gps-vetting.MARGIN, gps+acquire.DURATION+vetting.MARGIN) if not control else [],
              "slide_columns": EPOCH_FIELDS, "slides": []}
    for s in SHIFTS:
        slide = core.scan(*pairs, shifts=(s,))[s]
        best = sorted(slide.triggers, key=lambda t: -t.score)
        count = sum(t.score >= THRESHOLD for t in best)
        report["slides"].append([s, slide.live_seconds, slide.nominal_tiles, slide.gated_tiles,
                                 len(best), count, best[0].score if best else 0.])
        if s:
            continue
        ledger(output / f"{gps}.csv.gz", gps, live, best)
        report["top"] = [describe(tr, gps) for tr in best[:5]]
        print("zero-lag", gps, report["HL_live_s"], report["top"][:1], "ge20", count, flush=True)
        follow = best[:1] if not control else [t for t in best if abs(t.time + gps - block["event_gps"]) <= 8][:1]
        if not control and count:
            follow = best[:count]
        report["follow_up"] = []
        for tr in follow:
            info = describe(tr, gps)
            try:
                independent = core.confirm(inputs[0][0], inputs[1][0], acquire.FS, tr.time,
                                          band=tr.band, duration=tr.duration,
                                          good_h=inputs[0][1], good_l=inputs[1][1])
                info["independent"] = asdict(independent)
            except ValueError as e:
                info["independent_unavailable"] = str(e)
            info["V1_public_good_second"] = bool(vgood[int(tr.time)]) if vgood is not None else False
            if control:
                info.update({"event": block["event"], "reference_gps": block["event_gps"],
                             "delta_s": info["gps"] - block["event_gps"],
                             "within_1s": abs(info["gps"] - block["event_gps"]) <= 1,
                             "rank_in_control_frame": best.index(tr)+1})
            report["follow_up"].append(info)
    write_json(output / f"{gps}.json", report)
    return report


def loudest_background(cache, reports, threshold):
    """Count actual off-time seconds above the observed campaign loudest score.

    Recreate gated ranks from the fixed primary search. Only epoch/slide pairs
    whose stored maximum can reach the threshold need a second evaluation.
    """
    counts = {s: 0 for s in SHIFTS[1:]}
    threshold_crossings = []
    for report in reports:
        if report["role"] != "search":
            continue
        rows = [r for r in report["slides"][1:] if r[-1] >= min(threshold, THRESHOLD)]
        if not rows:
            continue
        gps = report["gps_start"]
        pairs = []
        for d in ("H1", "L1"):
            _, strain, good = acquire.load_detector(cache / f"{d[0]}-{d}_GWOSC_{acquire.RELEASE}-{gps}-4096.hdf5")
            pairs.append(core.prepare(strain, acquire.FS, vetting.mask_seconds(gps, good)))
            del strain
        for row in rows:
            shift = row[0]
            slide = core.scan(*pairs, shifts=(shift,))[shift]
            if max(t.score for t in slide.triggers) != row[-1]:
                raise ValueError("background replay changed ranks")
            counts[shift] += sum(t.score >= threshold for t in slide.triggers)
            for tr in slide.triggers:
                if tr.score >= THRESHOLD:
                    threshold_crossings.append({**describe(tr, gps), "epoch_start": gps, "shift_s": shift,
                                                "L1_tile_gps": gps+tr.time+shift})
    return counts, threshold_crossings


def aggregate(reports):
    survey = [r for r in reports if r["role"] == "search"]
    if len(survey) != 10 or len({r["month"] for r in survey}) != 10:
        raise ValueError("require all ten distinct fixed survey months")
    combined = []
    for i, s in enumerate(SHIFTS):
        rows = [r["slides"][i] for r in survey]
        if any(row[0] != s for row in rows):
            raise ValueError("inconsistent slide order")
        combined.append([s, *(sum(row[j] for row in rows) for j in range(1, 6)), max(row[-1] for row in rows)])
    fg = combined[0]
    bg = combined[1:]
    top = sorted((tr for r in survey for tr in r["top"]), key=lambda t: -t["score"])[:5]
    exceed = sum(row[-1] >= fg[-1] for row in bg)
    return {"slide_columns": EPOCH_FIELDS, "campaign_slides": combined,
            "HL_foreground_s": fg[1], "background_s": sum(row[1] for row in bg),
            "foreground_seconds_ge20": fg[-2], "background_seconds_ge20": sum(row[-2] for row in bg),
            "foreground_loudest": fg[-1], "background_slides_equal_or_louder": exceed,
            "empirical_campaign_tail_fraction": (exceed + 1) / (len(bg) + 1),
            "foreground_nominal_tiles": fg[2], "background_nominal_tiles": sum(row[2] for row in bg),
            "foreground_gated_tiles": fg[3], "background_gated_tiles": sum(row[3] for row in bg),
            "foreground_lag_evaluations": fg[3]*len(core.LAGS),
            "background_lag_evaluations": sum(row[3] for row in bg)*len(core.LAGS),
            "background_slide_exposure_range_s": [min(row[1] for row in bg), max(row[1] for row in bg)],
            "V1_after_public_event_within_HL_live": sum(r["V1_after_public_event_within_HL_live"] for r in survey),
            "top_foreground": top, "epoch_companions": [r["gps_start"] for r in reports]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache", type=Path)
    parser.add_argument("--out", type=Path, default=Path(__file__).parent / "results")
    parser.add_argument("--selection", type=Path, help="frozen manifest, default: CACHE/search_selection.json")
    parser.add_argument("--only-gps", type=int)
    parser.add_argument("--aggregate-only", action="store_true")
    parser.add_argument("--tail-counts", action="store_true",
                        help="also replay qualifying background shifts for a stationary FAR diagnostic")
    args = parser.parse_args()
    output = args.out / "epochs"
    output.mkdir(parents=True, exist_ok=True)
    selection = json.loads((args.selection or args.cache / "search_selection.json").read_text())
    for block in selection["blocks"]:
        if not args.aggregate_only and (args.only_gps is None or block["gps"] == args.only_gps):
            run_epoch(args.cache, output, block)
    if args.only_gps is None:
        reports = [json.loads((output / f"{b['gps']}.json").read_text()) for b in selection["blocks"]]
        summary = aggregate(reports)
        if args.tail_counts:
            counts, crossings = loudest_background(args.cache, reports, summary["foreground_loudest"])
            summary["background_ge20_triggers"] = crossings
            summary["background_seconds_ge_loudest"] = sum(counts.values())
            summary["background_ge_loudest_by_shift"] = counts
        write_json(args.out / "survey.json", summary)
        print({k: v for k, v in summary.items() if k not in ("campaign_slides", "top_foreground")}, flush=True)


if __name__ == "__main__":
    main()
