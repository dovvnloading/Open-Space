#!/usr/bin/env python3
"""Independent numerical audit of the frozen Sector 105 research record.

Checks deposited CSV/JSON and per-trial injections without calling scan.py's rates().
Optional replay directories allow comparison to re-downloaded original FITS.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent


def table(p):
    with p.open(newline="") as f:
        return list(csv.DictReader(f))


def compressed_table(p):
    with gzip.open(p, "rt", newline="") as f:
        return list(csv.DictReader(f))


def near(x, y, tol=1e-8):
    if not math.isclose(float(x), float(y), rel_tol=0, abs_tol=tol):
        raise AssertionError(f"Mismatch: {x} versus {y}")


def audit_cohort(folder, count, first_rank, trigger_count, expected_trials):
    products = table(folder / "products.csv")
    audit = table(folder / "audit.csv")
    triggers = table(folder / "triggers.csv")
    trials = compressed_table(folder / "injections.csv.gz")
    report = json.loads((folder / "report.json").read_text())
    assert len(products) == len(audit) == count, (folder, len(products), len(audit))
    assert len(triggers) == trigger_count and len(trials) == expected_trials
    assert [int(p["rank"]) for p in products] == list(range(first_rank, first_rank+count))
    assert all(p["status"] == "ok" and len(p["sha256"]) == 64 and int(p["bytes"]) > 0 for p in products)
    assert len({p["tic"] for p in products}) == count
    p_by_id = {p["tic"]: p for p in products}
    a_by_id = {a["tic"]: a for a in audit}
    assert set(p_by_id) == set(a_by_id) == set(t["tic"] for t in trials)
    for t in triggers:
        assert t["tic"] in p_by_id and int(t["sector"]) == 105
    for a in audit:
        assert int(a["n_triggers"]) == sum(t["tic"] == a["tic"] for t in triggers)
        if folder.name == "confirmation":
            assert a["sha256"] == p_by_id[a["tic"]]["sha256"]
    assert len({t["tic"] for t in triggers}) == (39 if count == 256 else 14)
    assert all(int(t["sector"]) == 105 for t in trials)
    wall_days = sum(float(a["span_d"]) for a in audit)
    grouped = defaultdict(list)
    for t in trials:
        grouped[(t["tic"], round(float(t["depth"]), 4), int(t["duration_h"]))].append(int(t["recovered"]))
    assert len(grouped) == count * 9
    assert all(len(x) == 8 and set(x) <= {0,1} for x in grouped.values())
    results = {(round(float(x["depth"]), 4), int(x["duration_h"])): x for x in report["limits"]}
    assert len(results) == 9
    # Rebuild effective exposure, a weighted Hoeffding bound, and a Poisson
    # upper count bound without using any pipeline functions.
    upper_count = float(chi2.ppf(.975, 2*(trigger_count+1))/2)
    for (depth, duration), r in results.items():
        sub = [(tic, hits) for (tic, d, h), hits in grouped.items()
               if d == depth and h == duration]
        assert len(sub) == count
        recovered = sum(sum(hits) for _, hits in sub)
        exposure = sum(float(a_by_id[tic]["span_d"])*sum(hits)/len(hits)
                       for tic, hits in sub)
        margin = math.sqrt(0.5*math.log(1/.025)*sum(
            float(a_by_id[tic]["span_d"])**2/len(hits) for tic, hits in sub))
        lower = max(0., exposure-margin)
        assert recovered == int(r["recovered"])
        assert len(sub)*8 == int(r["trials"])
        near(exposure, r["effective_star_days"], 1e-7)
        near(lower, r["exposure_lower_97p5"], 1e-7)
        near(upper_count/lower, r["rate_upper_95_per_star_day"], 1e-10)
        assert int(r["all_trigger_count"]) == trigger_count
    if count == 256:
        assert int(report["registered"]) == 256 and int(report["n_analyzed"]) == 256
        assert int(report["failures"]) == 0 and int(report["stars_triggered"]) == 39
        near(wall_days, 7108.144747191158, 1e-5)
        two=results[(.02, 12)]
        one=results[(.01, 12)]
        assert (int(two["recovered"]),int(two["trials"])) == (1149,2048)
        near(two["effective_star_days"], 3987.919098888011, 1e-6)
        near(two["rate_upper_95_per_star_day"], .021948754507137385, 1e-10)
        assert (int(one["recovered"]),int(one["trials"])) == (466,2048)
        near(100*466/2048, 22.75390625, 1e-7)
    else:
        assert int(report["n_star_sectors"]) == 128 and int(report["n_distinct_tics"]) == 128
        near(wall_days, 3554.1, .2)
    return {"stars":count,"wall_star_days":round(wall_days,4),
            "events":len(triggers),"injection_trials":len(trials),
            "recovery_2pct_12h":f"{results[(.02,12)]['recovered']}/{results[(.02,12)]['trials']}",
            "upper_rate_2pct_12h":results[(.02,12)]["rate_upper_95_per_star_day"]}


def events_match(before, after):
    assert len(before) == len(after)
    def key(row):
        return (int(row["tic"]), float(row["peak_btjd"]), float(row["depth"]),
                int(row["duration_h"]), float(row["ingress_fraction"]))
    for first, second in zip(sorted(map(key,before)), sorted(map(key,after))):
        assert first[0] == second[0] and first[3] == second[3]
        near(first[1],second[1],1e-7)
        near(first[2],second[2],1e-7)
        near(first[4],second[4],1e-10)


def replay_check(path, cohort):
    before = HERE / cohort
    old = json.loads((before/"report.json").read_text())
    new = json.loads((path/"report.json").read_text())
    assert int(new["all_trigger_count"]) == int(old["all_trigger_count"])
    assert new["seed"] == old["seed"] and new["trials_per_cell"] == old["trials_per_cell"]
    if (path/"products.csv").exists():
        expected={r["tic"]:(r["sha256"],int(r["bytes"])) for r in table(before/"products.csv")}
        obtained={r["tic"]:(r["sha256"],int(r["bytes"])) for r in table(path/"products.csv")}
        assert obtained==expected, "Replayed FITS bytes differ from the original SHA-256 receipts"
    events_match(table(before/"triggers.csv"), table(path/"triggers.csv"))
    new_limits = {(float(x["depth"]),int(x["duration_h"])):x for x in new["limits"]}
    for original in old["limits"]:
        fresh=new_limits[(float(original["depth"]),int(original["duration_h"]))]
        assert (fresh["trials"],fresh["recovered"]) == (original["trials"],original["recovered"])
        near(fresh["effective_star_days"],original["effective_star_days"],1e-6)
        near(fresh["rate_upper_95_per_star_day"],original["rate_upper_95_per_star_day"],1e-8)
    print("Re-extracted FITS agrees with saved scan and injection recovery:",cohort)


def receipts():
    availability=json.loads((HERE/"availability.json").read_text())
    seen={e["url"].split("/")[-1]:e for e in availability}
    assert seen["s0105.csv"]["status"] == 200
    assert seen["s0106.csv"]["status"] == seen["s0107.csv"]["status"] == 404
    sample=json.loads((HERE/"sample.json").read_text())
    assert sample["sector"] == 105 and sample["raw_rows"] == sample["unique"] == 399217
    assert sample["sha256"] == seen["s0105.csv"]["sha256"]
    assert len(sample["targets"]) == 128
    pilot=table(HERE/"products.csv")
    confirmed=table(HERE/"confirmation"/"products.csv")
    allrows=pilot+confirmed
    assert len(allrows) == 384 and len({r["tic"] for r in allrows}) == 384
    assert [str(x["tic"]) for x in sample["targets"]] == [r["tic"] for r in pilot]
    heldout=table(HERE/"confirmation"/"audit.csv")
    mags=[float(p["tmag"]) for p in pilot] + [float(a["tessmag"]) for a in heldout]
    assert len(mags) == 384 and sum(m < 13 for m in mags) == 252
    near(min(mags),6.15,.02)
    near(max(mags),14.99,.02)
    earlier=json.loads(gzip.open(HERE/"confirmation"/"earlier-visit-receipts.json.gz","rt").read())
    confirmed_ids={r["tic"] for r in table(HERE/"confirmation"/"triggers.csv")}
    assert len(earlier) == 39 and {str(s["tic"]) for s in earlier} == confirmed_ids
    has=[s for s in earlier if any(v["status"]=="ok" for v in s["visits"])]
    good=[v for s in earlier for v in s["visits"] if v["status"] == "ok"]
    assert len(has) == 36 and len(good) == 75
    assert all(len(v["sha256"])==64 for v in good)
    # Check fetched catalog response bytes against recorded SHA-256.
    totals={}
    for path,expected in [(HERE/"pilot-catalog-replay.jsonl.gz",56),
                          (HERE/"confirmation"/"catalog-receipts.jsonl.gz",156)]:
        with gzip.open(path,"rt") as inp:
            replies=list(map(json.loads,inp))
        assert len(replies) == expected
        for reply in replies:
            if "response" in reply:
                raw=reply["response"].encode("utf-8")
                assert hashlib.sha256(raw).hexdigest() == reply["sha256"]
                assert len(raw) == int(reply["bytes"])
        totals[path.name]=len(replies)
    return {"target_roster_sha256":sample["sha256"],
            "sample_stars":len(allrows),"tmag_under_13":sum(m<13 for m in mags),
            "earlier_products":len(good),"earlier_stars":len(has),
            "s106_status":404,"s107_status":404,"catalog_receipt_counts":totals}


def source():
    candidates=table(HERE/"confirmation"/"triggers.csv")
    ev=next(t for t in candidates if t["tic"]=="316937197" and
            abs(float(t["peak_btjd"])-4217.63214)<.02)
    near(ev["depth"],.17154443348306608,1e-9)
    cutouts=json.loads((HERE/"confirmation"/"pixels.json").read_text())
    assert len(cutouts)==39
    center=next(s for s in cutouts if str(s["tic"])=="316937197")
    pixel=next(e for e in center["events"] if abs(e["peak_btjd"]-4217.63214)<.02)
    near(pixel["aperture_fractional_depth"],.007054948844162463,1e-10)
    assert pixel["centroid_offset_px"] < .3 and len(center["fits_sha256"])==64
    assert sum(e["aperture_fractional_depth"]>=.01 for s in cutouts
               for e in s["events"]) == 11
    history=json.loads((HERE/"confirmation"/"lead-history.json").read_text())
    timings=[t for t in history["timings"] if t["accepted"]]
    assert len(timings)==13
    x=np.array([int(t["cycle"]) for t in timings],dtype=float)
    y=np.array([float(t["btjd"]) for t in timings])
    period=np.dot(x-x.mean(),y-y.mean())/np.dot(x-x.mean(),x-x.mean())
    epoch=y.mean()-period*x.mean()
    residual=np.sqrt(np.mean((y-epoch-period*x)**2))
    near(period,18.767053552317485,1e-10)
    near(period,history["refined_clock"]["period_d"],1e-10)
    near(residual,history["refined_clock"]["rms_residual_d"],1e-9)
    near(y.max()-y.min(),history["refined_clock"]["baseline_d"],1e-6)
    cats=table(HERE/"confirmation"/"catalogs.csv")
    lead=next(t for t in cats if t["tic"]=="316937197")
    eb=next(t for t in cats if t["tic"]=="116658378")
    assert lead["vsx_status"]=="ok" and int(lead["vsx_sources_63arcsec"])==0
    assert lead["gaia_status"]=="ok" and lead["non_single_star"]=="0"
    assert lead["best_class_name"]=="" and lead["simbad_otype"]=="*"
    assert eb["best_class_name"]=="ECL" and eb["simbad_otype"]=="EB*"
    near(lead["brightest_neighbor_sep_arcsec"],1.564,.01)
    return {"TIC316937197_period_d":period,
            "TIC316937197_qlp_depth":float(ev["depth"]),
            "TIC316937197_pixel_depth":pixel["aperture_fractional_depth"],
            "TIC316937197_eclipsing_host":"unresolved; no VSX ECL or Gaia DR3 ECL/NSS match",
            "other_pixel_depths_at_least_1pct":11}


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pilot-replay",type=Path)
    ap.add_argument("--confirmation-replay",type=Path)
    ap.add_argument("--out",type=Path)
    args=ap.parse_args()
    results={"pilot":audit_cohort(HERE,128,1,34,9216),
             "confirmation":audit_cohort(HERE/"confirmation",256,129,65,18432),
             "receipts":receipts(),"lead":source()}
    if args.pilot_replay:
        replay_check(args.pilot_replay,".")
    if args.confirmation_replay:
        replay_check(args.confirmation_replay,"confirmation")
    if args.out:
        args.out.write_text(json.dumps(results,indent=2)+"\n")
    print(json.dumps(results,indent=2))
