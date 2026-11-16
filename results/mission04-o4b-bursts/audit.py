"""Independent reviewer replay for Mission 04. No primary ranking or rate code at ledger level.

python audit.py --results results --output results/audit_results.json
Add --pilot-cache CACHE to verify all pilot H/L source digests, masks, shifted
exposures and 40 injections with newly estimated spectra from signal-plus-noise.
"""
import argparse
from collections import defaultdict
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
from zipfile import ZipFile

import numpy as np
from scipy.optimize import brentq
from scipy.stats import binom, poisson

SHIFTS = [0] + [sign * lag for lag in range(4, 257, 4) for sign in (1, -1)]
DURATION = 4096
CUT, FLOOR = 20.0, 3.0


def verify(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path, algorithm="sha256"):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, algorithm).hexdigest()


def read_json(path):
    return json.loads(path.read_bytes())


def read_csv(data):
    return list(csv.DictReader(io.StringIO(data.decode("utf-8"))))


def close(a, b):
    return math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10)


def limits(seconds, counts, n, k, alpha):
    # Solve P(X>=k | efficiency) = alpha and P(Y<=counts | mean) = alpha.
    low = brentq(lambda e: binom.sf(k-1, n, e)-alpha, 0., 1.) if k else 0.
    high = brentq(lambda m: poisson.cdf(counts, m)-alpha, 0., 4*(counts+10))
    return {"efficiency_lower": low, "count_upper": high, "effective_seconds": seconds*low,
            "upper_per_day": high*86400/(seconds*low) if low else None}


def load_runs(intervals):
    mask = np.zeros(DURATION, dtype=bool)
    end = -1
    for a, b in intervals:
        verify(0 <= a < b <= DURATION and a > end, "Invalid or duplicate half-open run")
        mask[a:b] = True
        end = b
    return mask



def eligible_from_veto(good):
    # A live second needs one full adjacent second on either side and at least
    # eight interior training seconds in the other half of the 32 s epoch.
    guarded = np.zeros(DURATION, dtype=bool)
    guarded[1:-1] = good[:-2] & good[1:-1] & good[2:]
    halves = guarded.reshape(-1, 16)
    train = np.sum(halves[:, 1:-1], axis=1)
    opposite = np.arange(len(train)) ^ 1
    return (halves & (train[opposite] >= 8)[:, None]).ravel()


def shifted_exposures(h, l):
    # Independent non-circular implementation. For s>0 pair H[i] with L[i+s].
    return [int(np.count_nonzero(h[max(-s, 0):min(DURATION, DURATION-s)] &
                                  l[max(s, 0):min(DURATION, DURATION+s)])) for s in SHIFTS]


def one_epoch(archive, prefix, block):
    def member(name):
        return archive.read(prefix + name)

    gps = block["gps"]
    rec = json.loads(member("epoch.json"))
    verify((rec["gps"], rec["month"]) == (gps, block["month"]), "Epoch identity mismatch")
    detector_masks = []
    for ifo in ("H1", "L1"):
        m = rec["masks"][ifo]
        before, after = (load_runs(m[x]) for x in ("post_veto_runs", "searchable_runs"))
        verify(np.array_equal(after, eligible_from_veto(before)), "Searchable mask differs from veto/PSD selection")
        verify(before.sum() == m["cumulative_seconds"]["after_public_event_mask"], "Veto count mismatch")
        detector_masks.append(after)
    h, l = detector_masks
    live = h & l
    exposure = shifted_exposures(h, l)
    slides = rec["slides"]
    verify(len(slides) == len(SHIFTS), "Incomplete shifts")
    for s, secs, row in zip(SHIFTS, exposure, slides):
        verify(row[:3] == [s, secs, secs*88], f"Shift exposure or tile mismatch {gps}/{s}")

    fg = read_csv(member("foreground.csv"))
    seconds = [int(r["offset_second"]) for r in fg]
    verify(np.array_equal(seconds, np.flatnonzero(live)), "Foreground seconds differ from H/L mask")
    scores = {sec: float(r["score"]) for sec, r in zip(seconds, fg)}
    verify(all(np.isfinite(t) and t >= 0 for t in scores.values()), "Invalid foreground score")
    verify(exposure[0] == rec["live_s"] and max(scores.values(), default=0) == slides[0][-1],
           "Foreground maximum/exposure mismatch")

    by_shift = defaultdict(dict)
    for row in read_csv(member("tail.csv")):
        s = int(row["shift_s"])
        time = float(row["H1_peak_gps"])
        sec = math.floor(time) - gps
        val = float(row["score"])
        verify(s in SHIFTS and 0 <= sec < DURATION and 0 <= sec+s < DURATION and
               h[sec] and l[sec+s] and float(row["L1_tile_gps"]) == time+s, "Tail outside live pair")
        verify(val >= FLOOR and np.isfinite(val) and sec not in by_shift[s], "Invalid sparse tail row")
        by_shift[s][sec] = val
    verify(by_shift[0] == {k:v for k,v in scores.items() if v >= FLOOR}, "Foreground sparse mismatch")
    for s, row in zip(SHIFTS, slides):
        tail = by_shift[s].values()
        verify(sum(v >= CUT for v in tail) == row[-2], "Threshold count mismatch")
        verify(max(tail, default=0.) == (row[-1] if row[-1] >= FLOOR else 0.), "Tail maximum mismatch")
    injections = json.loads(member("injections.json"))
    for t in injections:
        sec = t["second"]
        verify(0 <= sec < DURATION and t["gps_start"] == gps and
               t["live"] == bool(live[sec]) and t["baseline"] == scores.get(sec), "Invalid injection epoch")
        if t["live"]:
            injected = t["injected_score"]
            verify(np.isfinite(injected) and injected >= 0 and
                   (scores[sec] < CUT <= injected) == t["recovered"], "Invalid injection decision")
        else:
            verify(t["injected_score"] is None and not t["recovered"], "Recovered vetoed draw")
    return exposure, by_shift, injections, rec["configuration_sha256"]


def audit_holdout(results, selection, pilot_gps):
    root = results / "holdout"
    camp = read_json(root/"campaign.json")
    sel_path = results/"holdout_selection.json"
    verify(sha(sel_path) == camp["selection_sha256"], "Holdout selection checksum mismatch")
    blocks = {b["gps"]: b for b in selection["blocks"]}
    verify(len(blocks) == len(selection["blocks"]) == camp["n_frames"], "Duplicate or missing blocks")
    starts = sorted(blocks)
    verify(all(b-a >= DURATION for a,b in zip(starts, starts[1:])) and
           all(abs(g-p) >= DURATION for g in starts for p in pilot_gps), "Overlapping epochs")
    exp = np.zeros(len(SHIFTS), dtype=int)
    tails = defaultdict(list)
    trials, receipts = [], {}
    seen = set()
    for receipt in camp["monthly_ledgers"]:
        path = root / receipt["file"]
        digest = sha(path)
        verify(path.stat().st_size == receipt["bytes"] and digest == receipt["sha256"], "Monthly archive hash mismatch")
        verify(path.name not in receipts, "Duplicate monthly receipt")
        receipts[path.name] = digest
        with ZipFile(path) as z:
            names = z.namelist()
            frame_set = {int(name.split("/")[0]) for name in names}
            expected = {f"{g}/{name}" for g in frame_set for name in
                        ("epoch.json", "foreground.csv", "tail.csv", "injections.json")}
            verify(len(names) == len(expected) and set(names) == expected and z.testzip() is None,
                   "Missing, duplicate or corrupt archive member")
            for gps in sorted(frame_set):
                verify(gps in blocks and gps not in seen and path.name == f"ledgers-{blocks[gps]['month']}.zip",
                       "Duplicate, missing, or misplaced epoch")
                seen.add(gps)
                exposure, by_shift, records, fingerprint = one_epoch(z, f"{gps}/", blocks[gps])
                verify(fingerprint == camp["configuration_sha256"], "Configuration identifier mismatch")
                exp += exposure
                for s, values in by_shift.items():
                    tails[s].extend(values.values())
                trials.extend(records)
    verify(seen == set(blocks), "Unverified holdout epochs")
    by_trial = sorted(trials, key=lambda t: t["trial"])
    verify([t["trial"] for t in by_trial] == list(range(8192)), "Missing or duplicate random trial")
    rng = np.random.default_rng(404150)
    for record in by_trial:
        b, sec = divmod(int(rng.integers(len(starts)*DURATION)), DURATION)
        offset, z = rng.random(), rng.uniform(-1, 1)
        az, psi, phase = rng.random(3) * [2*np.pi, np.pi, 2*np.pi]
        verify((record["gps_start"], record["second"], record["offset"], record["cos_theta"],
                record["azimuth"], record["psi"], record["phase"]) == (starts[b], sec, offset, z, az, psi, phase),
               "Sky/time draw differs from fixed-seed schedule")
    high = max(tails[0], default=0.)
    verify(high >= FLOOR, "Sparse tail does not reach campaign foreground maximum")
    maxes = [max(tails[s], default=0.) for s in SHIFTS]
    verify(len(SHIFTS) == len(camp["campaign_slides"]), "Incomplete published campaign shifts")
    for i, row in enumerate(camp["campaign_slides"]):
        verify(row[0] == SHIFTS[i] and row[1] == exp[i] and (row[-1] if row[-1] >= FLOOR else 0.) == maxes[i],
               "Campaign shift differs from ledgers")
    louder = sum(x >= high for x in maxes[1:])
    count = sum(s >= CUT for s in tails[0])
    background = [v for shift in SHIFTS[1:] for v in tails[shift]]
    rank = (louder+1)/(len(SHIFTS))
    verify(camp["HL_live_s"] == exp[0] and camp["background_s"] == sum(exp[1:]) and
           camp["foreground_loudest"] == high and camp["background_seconds_ge_loudest"] == sum(s>=high for s in background)
           and camp["foreground_ge20"] == count and camp["background_ge20"] == sum(s>=CUT for s in background)
           and camp["background_campaigns_ge_loudest"] == louder and camp["campaign_rank_fraction"] == rank,
           "Campaign significance/exposure mismatch")
    scheduled = len(blocks)*DURATION
    k = sum(t["recovered"] for t in trials)
    rate = camp["rate"]
    verify(rate["scheduled_s"] == scheduled == camp["scheduled_s"] and rate["draws"] == len(trials) and
           rate["recovered"] == k and rate["foreground_seconds_ge20"] == count, "Invalid population counts")
    out = {"epochs": len(blocks), "foreground_seconds": int(exp[0]), "background_seconds": int(sum(exp[1:])),
           "scheduled_seconds": scheduled, "shifts_per_epoch": len(SHIFTS)-1, "shifted_seconds_at_threshold":
           sum(v >= CUT for v in background), "foreground_maximum": high, "shifted_maxima_ge_foreground": louder,
           "campaign_rank_fraction": rank, "total_draws": len(trials), "searchable_draws": sum(t["live"] for t in trials),
           "recovered": k, "foreground_seconds_at_threshold": count, "shard_sha256": receipts}
    for key, alpha in (("pointwise90", .05), ("two_campaign_simultaneous90", .025)):
        lim = limits(scheduled, count, len(trials), k, alpha)
        for a,b in (("effective_time_s_lower","effective_seconds"), ("count_upper","count_upper"),
                    ("efficiency_lower","efficiency_lower"), ("rate_per_day_upper","upper_per_day")):
            verify(close(rate[key][a], lim[b]), f"Holdout limit differs: {key}/{a}")
        out[key] = lim
    return out


def audit_pilot(results, pilot_manifest, cache):
    summ, sensitivity = (read_json(results / f"{x}.json") for x in ("survey", "sensitivity"))
    blocks = [b for b in pilot_manifest["blocks"] if b["role"] == "search"]
    scores, live_by_block, shiftrows = {}, {}, []
    for b in blocks:
        g = b["gps"]
        rec = read_json(results/"epochs"/f"{g}.json")
        with gzip.open(results/"epochs"/f"{g}.csv.gz", "rt") as stream:
            fg = list(csv.DictReader(stream))
        seconds = [int(row["gps_second"]) for row in fg]
        verify(len(seconds) == len(set(seconds)) == rec["HL_live_s"] and
               all(g <= s < g+DURATION for s in seconds), "Invalid pilot exposure")
        scores.update({int(r["gps_second"]): float(r["score"]) for r in fg})
        live_by_block[g] = set(seconds)
        verify([r[0] for r in rec["slides"]] == SHIFTS and
               rec["slides"][0][-1] == max(float(r["score"]) for r in fg), "Invalid pilot shift order or maximum")
        shiftrows.append(rec["slides"])
    maxes = [max(r[i][-1] for r in shiftrows) for i in range(len(SHIFTS))]
    fgmax = max(scores.values())
    bg = sum(row[1] for block in shiftrows for row in block[1:])
    rank = (1+sum(v>=fgmax for v in maxes[1:]))/len(SHIFTS)
    verify(len(scores) == summ["HL_foreground_s"] == sensitivity["HL_live_s"] and summ["background_s"] == bg and
           fgmax == summ["foreground_loudest"] and rank == summ["empirical_campaign_tail_fraction"],
           "Pilot summary mismatch")
    count = sum(v >= CUT for v in scores.values())
    verify(count == summ["foreground_seconds_ge20"] == sensitivity["foreground_seconds_ge_threshold"],
           "Pilot foreground count mismatch")
    with gzip.open(results/"injections.csv.gz", "rt") as stream:
        trials = list(csv.DictReader(stream))
    out = {"epochs": len(blocks), "foreground_seconds": len(scores), "background_seconds": bg,
           "foreground_maximum": fgmax, "campaign_rank_fraction": rank, "foreground_seconds_at_threshold": count, "cases": {}}
    targets = {"primary":2048, "lower":512, "middle":512, "higher_frequency":512}
    verify(set(sensitivity["cases"]) == set(targets) and len(trials) == sum(targets.values()),
           "Missing or unexpected pilot population")
    for case, summary in sensitivity["cases"].items():
        rows = [t for t in trials if t["case"] == case]
        verify(len(rows) == targets[case] == summary["n"] and
               sorted(int(r["trial"]) for r in rows) == list(range(targets[case])), "Pilot trial IDs mismatch")
        k = 0
        for r in rows:
            sec = int(r["block"])+int(r["second"])
            actual = scores[sec] < CUT <= float(r["injected_score"])
            verify(scores[sec] == float(r["baseline"]) and int(actual) == int(r["recovered"]), "Pilot injection mismatch")
            k += actual
        lim = limits(len(scores), count, len(rows), k, .05)
        verify(k == summary["recovered"] and close(lim["efficiency_lower"], summary["lower_efficiency_95"]) and
               close(lim["upper_per_day"], summary["rate_upper_per_day_pointwise90"]), f"Pilot rate mismatch: {case}")
        out["cases"][case] = {"n": len(rows), "k": k, **lim}
    if cache is not None:
        out["raw_replay"] = raw_pilot(cache, pilot_manifest, blocks, trials, live_by_block, shiftrows)
    return out


def raw_pilot(cache, manifest, blocks, trials, live_by_block, shiftrows):
    import acquire
    import core
    import vetting

    errors, source_hashes = [], {}
    for block, slides in zip(blocks, shiftrows):
        gps = block["gps"]
        rows = [r for r in trials if int(r["block"]) == gps]
        picks = [min(rows, key=lambda r: abs(float(r["injected_score"])-CUT)),
                 min(rows, key=lambda r: min(float(r["offset"]), 1-float(r["offset"]))),
                 min(rows, key=lambda r: min((int(r["second"])+float(r["offset"]))%16,
                                              16-(int(r["second"])+float(r["offset"]))%16)),
                 max(rows, key=lambda r: abs(float(r["delay_s"])))]
        arrays, masks, live = [], [], []
        for ifo in ("H1", "L1"):
            name = block["files"][ifo]
            path = cache / name
            digest = sha(path, "md5")
            verify(digest == manifest["manifest"][name]["md5"], "Pilot raw checksum mismatch")
            source_hashes[name] = digest
            origin, strain, good = acquire.load_detector(path)
            verify(origin == gps, "Pilot raw origin mismatch")
            good = vetting.mask_seconds(gps, good)
            live.append(eligible_from_veto(good))
            arrays.append(strain)
            masks.append(good)
        verify(set(gps+np.flatnonzero(live[0]&live[1])) == live_by_block[gps], "Pilot masks differ from ledger")
        verify(shifted_exposures(*live) == [row[1] for row in slides], "Pilot shifted raw masks differ")
        for r in picks:
            sec = int(r["second"])
            first, stop = max(0, 32*(sec//32)-2), min(DURATION, 32*(sec//32+1)+2)
            locally_prepared = []
            for i, (raw, good) in enumerate(zip(arrays, masks)):
                delay = 0. if i == 0 else float(r["delay_s"])
                wave = core.sine_gaussian(4096, stop-first, sec+float(r["offset"])-first+delay,
                            float(r["frequency"]), 9., float(r["hrss"]), float(r["phase"]))
                injected = raw[first*4096:stop*4096] + float(r[("F_H","F_L")[i]]) * wave
                # Recompute the complete local PSD from signal-plus-noise, as distinct from
                # the production short-crop injection optimized with a stored opposite-half PSD.
                locally_prepared.append(core.prepare(injected, 4096, good[first:stop], first_second=first))
            response = core.scan(*locally_prepared)[0]
            val = max((t.score for t in response.triggers if int(t.time) == sec), default=0.)
            err = abs(float(r["injected_score"])-val)
            verify(err < 1e-4 and (val >= CUT) == bool(int(r["recovered"])), "Raw pilot replay mismatch")
            errors.append(err)
    return {"source_md5":source_hashes, "raw_injection_trials":len(errors), "maximum_score_abs_error":max(errors),
            "detector_mask_shift_intersections_checked":len(blocks)*len(SHIFTS)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", type=Path, default=Path(__file__).parent/"results")
    parser.add_argument("--pilot-selection", type=Path)
    parser.add_argument("--pilot-cache", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    pilot = read_json(args.pilot_selection or args.results/"search_selection.json")
    holdout = read_json(args.results/"holdout_selection.json")
    out = {"pilot":audit_pilot(args.results,pilot,args.pilot_cache),
           "holdout":audit_holdout(args.results,holdout,[b["gps"] for b in pilot["blocks"]])}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(out, indent=2, allow_nan=False)+"\n")


if __name__ == "__main__":
    main()
