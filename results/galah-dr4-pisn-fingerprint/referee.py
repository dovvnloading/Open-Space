#!/usr/bin/env python3
"""Independent, stdlib-only referee of frozen GALAH DR4 screen receipts.

No import from screen.py, no scientific use of a missing or upper-limit value.
This is an independent arithmetic review of deposited rows, NOT an independent
download of the official 917,588-row source or a PISN-yield likelihood fit.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
RUN = HERE / "run"
CORE = ("na", "al", "mg", "co", "zn")
NEUTRON = ("ba", "sr")


def load(path: Path):
    return json.loads(path.read_text())


def must(test, message):
    if not test:
        raise AssertionError(message)


def maybe_float(s):
    if s in (None, "", "None"):
        return None
    x = float(s)
    must(math.isfinite(x), f"nonfinite value: {s}")
    return x


def pass_cuts(x, c, delta=0):
    return (
        x["na_minus_mg"] <= c["na_minus_mg_max"] + delta
        and x["al_minus_mg"] <= c["al_minus_mg_max"] + delta
        and x["co_fe"] <= c["co_fe_max"] + delta
        and x["zn_fe"] <= c["zn_fe_max"] + delta
    )


def review():
    c = load(HERE / "criteria.json")
    s = load(RUN / "scan_summary.json")
    control = load(HERE / "j1010_xing_table1.json")
    rows = list(csv.DictReader((RUN / "complete_five_rows.csv").open(newline="")))
    hits = list(csv.DictReader((RUN / "five_gate_survivors.csv").open(newline="")))
    must(s["rows"]["catalogue"] == 917588, "unexpected catalogue size")
    must(s["rows"]["metal_poor_quality_base"] == 5189, "base count differs from frozen first official run")
    must(s["rows"]["all_five_measured"] == len(rows), "five-element row count")
    must(s["rows"]["all_seven_measured"] <= len(rows), "seven-element denominator")
    must(s["rows"]["five_pattern_pass"] == len(hits), "survivors table disagrees")
    must(s["criteria_sha256"] == hashlib.sha256((HERE / "criteria.json").read_bytes()).hexdigest(),
         "criteria changed since archived computation")
    must(s["data_receipt"]["url"] == s["data_receipt"]["http"]["response_url"],
         "catalogue source redirect discrepancy")
    must(s["data_receipt"]["bytes"] > 700000000, "truncated input")
    must(len(s["data_receipt"]["sha256"]) == 64, "no pinned SHA-256")
    ids = set()
    n7, n5, n2, n5_relax_01, n5_relax_02, nmet2, nmet25, nco, nzn = 0, 0, 0, 0, 0, 0, 0, 0, 0
    min_co = float("inf")
    min_zn = float("inf")
    closest = None
    ccut = c["gate"]
    for r in rows:
        ident = r["sobject_id"]
        must(ident not in ids, f"duplicate star {ident}")
        ids.add(ident)
        must(r["gaia_dr3_source_id"].isdigit(), f"missing cross-ID {ident}")
        fe = maybe_float(r["fe_h"])
        must(c["metallicity_dex"]["min"] <= fe <= c["metallicity_dex"]["max"], f"metallicity not allowed {ident}")
        must(maybe_float(r["snr_px_ccd3"]) > c["quality"]["snr_px_ccd3_gt"], f"bad SNR {ident}")
        nmet2 += fe <= -2
        nmet25 += fe <= -2.5
        x = {}
        for el in CORE:
            x[f"{el}_fe"] = maybe_float(r[f"{el}_fe"])
            e = maybe_float(r[f"e_{el}_fe"])
            must(x[f"{el}_fe"] is not None and e is not None
                 and 0 < e <= c["quality"]["per_element_sigma_dex_max"],
                 f"bad core detection {el} for {ident}")
        x["na_minus_mg"] = x["na_fe"] - x["mg_fe"]
        x["al_minus_mg"] = x["al_fe"] - x["mg_fe"]
        must(abs(x["na_minus_mg"] - maybe_float(r["na_minus_mg"])) < 1e-8, "Na/Mg arithmetic")
        must(abs(x["al_minus_mg"] - maybe_float(r["al_minus_mg"])) < 1e-8, "Al/Mg arithmetic")
        ok5 = pass_cuts(x, ccut)
        must(ok5 == (r["five_gate"] == "True"), f"five gate incorrect for {ident}")
        n5 += ok5
        n2 += (x["na_minus_mg"] <= ccut["na_minus_mg_max"]
               and x["al_minus_mg"] <= ccut["al_minus_mg_max"])
        n5_relax_01 += pass_cuts(x, ccut, .1)
        n5_relax_02 += pass_cuts(x, ccut, .2)
        nco += x["co_fe"] <= ccut["co_fe_max"]
        nzn += x["zn_fe"] <= ccut["zn_fe_max"]
        min_co = min(min_co, x["co_fe"])
        min_zn = min(min_zn, x["zn_fe"])
        if closest is None or x["co_fe"] < closest["co_fe"]:
            closest = {"sobject_id":ident,"gaia_dr3_source_id":r["gaia_dr3_source_id"],
                       "fe_h":fe, "co_fe":x["co_fe"],"zn_fe":x["zn_fe"],
                       "na_minus_mg":x["na_minus_mg"],"al_minus_mg":x["al_minus_mg"]}
        nvalues = []
        for el in NEUTRON:
            val, error = maybe_float(r[f"{el}_fe"]), maybe_float(r[f"e_{el}_fe"])
            if val is None or error is None:
                nvalues.append(None)
            else:
                must(0 < error <= c["quality"]["per_element_sigma_dex_max"], "bad n-capture uncertainty")
                nvalues.append(val)
        seven_detected = all(v is not None for v in nvalues)
        n7 += seven_detected
        ok7 = ok5 and seven_detected and nvalues[0] <= ccut["ba_fe_max"] and nvalues[1] <= ccut["sr_fe_max"]
        must(ok7 == (r["seven_gate"] == "True"), f"seven gate incorrect for {ident}")
    must(n5 == 0 and n5 == s["rows"]["five_pattern_pass"], "five-element null not reproduced")
    must(n7 == 14 and n7 == s["rows"]["all_seven_measured"], "seven-element count not reproduced")
    must(n2 == s["rows"]["odd_even_3_pass_thresholds"], "odd-even subset inconsistent")
    must(nco == 0, "unexpected Co-depleted source")
    must(nmet25 == 0, "claimed coverage below -2.5 incorrect")
    for key, elset in (("five", CORE), ("seven", CORE + NEUTRON)):
        p = s["injection"][key]
        must(p["eligible_stars"] == (len(rows) if key == "five" else n7), "bad injection denominator")
        expected = p["eligible_stars"] * p["trials_per_star"]
        must(0 <= p["passing_trials"] <= expected, "bad injection successes")
        must(abs(p["passing_trials"]/expected - p["completeness_conditional_on_flags"]) < 1e-12,
             "injection fraction inconsistent")
    must(control["original"]["na"]["upper"] - control["original"]["mg"]["value"] < -1.35,
         "control Na/Mg upper not recovered")
    must(control["original"]["al"] is None, "control Al was not detected")
    must("upper" in control["original"]["zn"], "control Zn is a limit, not a detection")
    bound95 = 1 - math.pow(.05, 1/len(rows))
    outcome = {
        "audit_date_utc": datetime.now(timezone.utc).isoformat(),
        "status": "PASS: independent row-level threshold and detection/limit audit",
        "scope": "The 196 published complete-five detection rows and their frozen upstream aggregate receipt",
        "not_independently_verified": [
            "Re-fetching and filtering all 917588 official FITS rows",
            "Heger & Woosley 2002 and Takahashi et al. 2018 PISN versus CCSN/Ia yield likelihoods",
            "APOGEE DR17/LAMOST stellar cross-match, RV/line/blend vetting, SIMBAD/SAGA literature"
        ],
        "source_sha256": s["data_receipt"]["sha256"],
        "catalogue_rows_upstream_reported": s["rows"]["catalogue"],
        "metal_poor_quality_base_upstream_reported": s["rows"]["metal_poor_quality_base"],
        "five_measured_recounted": len(rows),
        "five_pattern_pass_recounted": n5,
        "seven_measured_recounted": n7,
        "two_odd_even_pass_recounted": n2,
        "co_fe_below_gate_recounted": nco,
        "zn_fe_below_gate_recounted": nzn,
        "fe_h_below_minus2_recounted": nmet2,
        "fe_h_below_minus2point5_recounted": nmet25,
        "five_all_thresholds_relaxed_by_0point1_dex": n5_relax_01,
        "five_all_thresholds_relaxed_by_0point2_dex": n5_relax_02,
        "minimum_co_fe": min_co,
        "minimum_zn_fe": min_zn,
        "most_co_depleted_row": closest,
        "conditional_binomial_95pct_upper_if_identical_sampling": bound95,
        "injection_5_detection_conditional": s["injection"]["five"]["completeness_conditional_on_flags"],
        "injection_7_detection_conditional": s["injection"]["seven"]["completeness_conditional_on_flags"],
        "referee_conclusion": "A real, reproducible measured-element screening null (0/196), NOT a physical PISN-free population or a discovery; the low-metallicity and Sr/Al coverage is inadequate for such an inference.",
        "science_label": "coverage-limited screen, no [REMARKABLE]"
    }
    path = HERE / "referee_report.json"
    path.write_text(json.dumps(outcome, indent=2, sort_keys=True) + "\n")
    print("REFEREE PASS", outcome["five_pattern_pass_recounted"], "/", len(rows),
          "complete-five; only", nmet2, "at [Fe/H] <= -2; zero at <= -2.5")
    return outcome


if __name__ == "__main__":
    review()
