"""Independent, offline consistency audit of the 17 November 2026 deposited tables.

Checks the numerical census, input receipts, source provenance, paired injections,
positive-control data bytes, decisive candidate comparisons, and catalog claims.
No network is needed; run from any directory with Python 3.10+.
"""
from collections import Counter
from pathlib import Path
import csv
import gzip
import hashlib
import json
import math
import re

ROOT = Path(__file__).resolve().parent
SHA = re.compile(r"[0-9a-f]{64}\\Z")


def rows(rel, delimiter=","):
    path = ROOT / rel
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def obj(rel):
    return json.loads((ROOT / rel).read_text())


def near(actual, expected, tol=0.01):
    assert abs(actual - expected) <= tol, (actual, expected, tol)


def run():
    coverage = obj("tables/coverage.json")
    fits = rows("tables/grid-fits.csv.gz")
    assert len(fits) == sum(r["inventoried"] for r in coverage) == 10976
    successful = [r for r in fits if r["fit_status"] == "ok"]
    assert len(successful) == sum(r["fitted"] for r in coverage) == 10971
    assert sum(int(r["n"]) for r in successful) == sum(r["usable_epochs"] for r in coverage) == 63064338
    assert sum(int(r["fit_status"] != "ok") for r in fits) == 5
    bysurvey = Counter(r["survey"] for r in fits)
    assert bysurvey == {"KMT": 6789, "OGLE": 2429, "MOA": 699, "PRIME": 970, "Gaia": 89}
    assert all(r["year"] in ("2024", "2025") for r in fits if r["survey"] == "KMT")
    idx = obj("kmt/index-receipts.json")["seasons"]
    assert [(r["year"], r["listed"]) for r in idx] == [(2024,3441),(2025,3348),(2026,2248)]
    assert "proprietary" in idx[2]["release_label"].lower()
    assert len(rows("tables/survey-id-groups.csv.gz")) == 9035
    assert len(rows("other/known_roster.csv")) == 34
    matched = rows("other/known_matches.csv")
    assert len({r["known_id"] for r in matched if r["status"] == "known by published ID"}) == 5

    # Compare every deposited screen input SHA against its independent survey
    # source receipt, including KMT's pySIS 404 -> DIAPL fallback cases.
    receipts = {}
    for year in (2024, 2025):
        for r in rows(f"kmt/pySIS-{year}.tsv", "\\t"):
            if r["sha256"]:
                receipts[r["event"]] = r["sha256"]
    for rec in obj("kmt/fallback-receipts.json"):
        if rec.get("curve", {}).get("sha256") and rec["event"] not in receipts:
            receipts[rec["event"]] = rec["curve"]["sha256"]
    for r in rows("ogle/photometry.csv"):
        receipts["OGLE-" + r["event"]] = r["sha256"]
    for r in rows("other/photometry.csv"):
        receipts[r["event"]] = r["input_sha256"]
    for r in rows("other/gaia.csv"):
        receipts[r["event"]] = r["sha256"]
    bad = [(r["event"], r["input_sha256"], receipts.get(r["event"])) for r in fits
           if not SHA.fullmatch(r["input_sha256"]) or r["input_sha256"] != receipts.get(r["event"])]
    assert not bad, bad[:8]

    ogle = rows("ogle/photometry.csv")
    assert len(ogle) == 2429 and sum(int(r["n"]) for r in ogle) == 4339671
    assert sum(int(r["bytes"]) for r in ogle) == 165027725
    assert all(r["status"] == "ok" for r in ogle)

    injections = rows("tables/injections.csv.gz")
    assert len(injections) == 8640
    assert len({(r["event"],r["t0"],r["te_d"],r["rho"]) for r in injections}) == 8640
    summary = rows("tables/injection-summary.csv")
    assert len(summary) == 30 and sum(int(r["trials"]) for r in summary) == 8640
    assert sum(int(r["recovered"]) for r in summary) == sum(int(r["paired_pass"]) for r in injections)
    for te, rho, n in [(0.03,0,8),(.03,2,23),(.10,0,28),(.10,2,69)]:
        row = next(r for r in summary if r["survey"]=="KMT" and float(r["te_d"])==te and float(r["rho"])==rho)
        assert (int(row["recovered"]), int(row["trials"]), int(row["alerts"])) == (n,288,24)

    control = obj("controls/earth-recovery.json")
    raw = gzip.decompress((ROOT/"controls/earth-figure1.tsv.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == control["response_sha256"]
    assert control["input_rows"] == 17444 and control["observations_within_1d"] == 90
    near(control["finite_no_blend"]["te_d"] * 1440, 39.1, .10)
    near(control["finite_no_blend"]["rho"], 3.21, .01)
    near(control["point_no_blend"]["chi2"]-control["finite_no_blend"]["chi2"], 42, 1)

    short = obj("peaks/KMT-2024-BLG-3326.json")["models"]
    near(short["fs"]["te_or_rise_d"], .085, .002)
    near(short["fs"]["rho"],2.61,.01)
    near(short["fs"]["chi2"]-short["bin"]["chi2"],5,.2)
    near(short["profile_long"]["chi2"]-short["fs"]["chi2"],1.6,.1)
    near(short["profile_small"]["chi2"]-short["fs"]["chi2"],1.6,.1)
    assert not short["bin"]["converged"], "3326 sampled binary optimizer did not converge"
    recurrence = {r["event"]:r for r in obj("kmt/extended/recurrence.json")["events"]}
    for name,te,rho in [("KMT-2024-BLG-0824",7.90,2.06),("KMT-2024-BLG-3278",8.16,2.58)]:
        g = recurrence[name]["fixed_geometry"]
        near(g["te_or_rise_d"],te,.02)
        near(g["rho"],rho,.02)
    og = obj("other/ogle0670-models.json")["models"]
    near(og["fs"]["te_or_rise_d"],5.13,.01)
    near(og["fs"]["rho"],1.42,.01)
    near(og["profile_short"]["chi2"]-og["fs"]["chi2"],2425,1)
    assert og["bin"]["converged"]
    near(og["fs"]["chi2"]-og["bin"]["chi2"],45.45,.05)
    b = rows("other/ogle0670_fixed_binary.csv")
    near(float(b[0]["fs_minus_fixed_bin"]),45.451,.001)
    near(float(b[-1]["fs_minus_fixed_bin"]),43.754,.001)
    win = { (r["year"],r["field"]):r for r in rows("ogle/window_by_field.csv")}
    near(float(win["2025","ALL"]["2h_PSPL"])*100,1.19,.01)
    near(float(win["2026","ALL"]["2h_PSPL"])*100,.94,.01)
    prior_replay = rows("tables/peak-replay.csv")
    assert prior_replay and all(r["status"]=="ok" and float(r["max_abs_chi2"])<=.05 for r in prior_replay)

    assert "[REMARKABLE]" not in (ROOT/"README.md").read_text()
    result = dict(status="passed",alerts=len(fits),fitted=len(successful),
                  points=sum(int(r["n"]) for r in successful),receipt_hashes=len(fits),
                  distinct_groups=9035,kmt_2026_public_pages=2248,
                  injections=len(injections),positive_control_te_min=round(control["finite_no_blend"]["te_d"]*1440,3),
                  prior_peak_replays=len(prior_replay),
                  note="Raw input receipts are compared to the published SHA-256 manifest; independently downloaded input bytes are checked by referee_replay.py.")
    print(json.dumps(result,indent=2))
    return result


if __name__ == "__main__":
    run()
