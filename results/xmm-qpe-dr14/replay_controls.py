#!/usr/bin/env python3
"""No-network replay of previously downloaded, position-vetted QPE EPIC-pn controls.

Replay uses the same frozen finder as control_validation.py but never requests
science archives. This is calibration only, NOT a host-galaxy archive screen.
"""
import csv
import json
from pathlib import Path
import numpy as np
from control_validation import BASE, OUT, TARGETS, BINSIZE_S, checksum, find_eruptions, verify_synthetic, store

# Directly from Webbe & Young (2023), Table 1 (doi:10.1093/rasti/rzad015),
# QPE-active XMM obs IDs: GSN5, RX8 bright eruptions, eRO1 2, eRO2 9.
# From Arcodia et al. (2024), A&A 684 A64, eRO3 2 and eRO4 3.
EXPECTED = {"GSN 069": 5, "RX J1301.9+2747": 8, "eRO-QPE1": 2,
            "eRO-QPE2": 9, "eRO-QPE3": 2, "eRO-QPE4": 3}
LIT = {"WebbeYoung2023": "https://doi.org/10.1093/rasti/rzad015",
       "Arcodia2024": "https://doi.org/10.1051/0004-6361/202348881",
       "Giustini2024": "https://doi.org/10.1051/0004-6361/202450861"}
# Giustini et al. (2024), A&A 692 A15, Table B.1: 2020B also contains
# 10 fitted/partial/weak eruptions, unlike the 8 bright events in Webbe&Young.
# Mission Reference Time (s from MJDref=50814 TT).
RX_FITTED_MET = [710877090,710884700,710901674,710906274,710924568,
                 710934969,710954337,710966143,710986288,711000219]

def from_csv(path):
    with path.open() as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise ValueError("No good 300s bins in saved control " + str(path))
    times = np.array([float(row["time_s"]) for row in rows])
    start = float(times[0])
    # The saved CSV omits *bad* PN bins, including the terminal GTI closure.
    # Retain a masked 3000s end window rather than turning the last valid bin
    # into an artificial time-series boundary; mask bins remain non-photometric.
    n = int(round((times[-1] - start) / BINSIZE_S))+1+10
    t = start + np.arange(n)*BINSIZE_S
    y,e,b = (np.full(n, np.nan) for _ in range(3))
    hard = np.full(n,np.nan)
    exposure = np.zeros(n)
    good = np.zeros(n,dtype=bool)
    for row in rows:
        i = int(round((float(row["time_s"])-start)/BINSIZE_S))
        if not 0 <= i < n:
            continue
        y[i]=float(row["net_0p2_2keV_cps"])
        e[i]=float(row["net_err_cps"])
        b[i]=float(row["scaled_background_0p2_2keV_cps"])
        hard[i]=float(row["hard_2_12keV_cps"]) if row["hard_2_12keV_cps"] else np.nan
        exposure[i]=float(row["exposure_s"])
        good[i]=True
    return {"t_s": t, "rate": y, "error": e, "background": b, "hard": hard,
            "exposure": exposure,"good":good,"bin_s":BINSIZE_S,"native_bin_s":None}

def rx_individual_match(peaks, max_abs_s=1200):
    """Best one-to-one match to *all 10* fitted RX peaks, not just the 8 bright."""
    actual=[(float(x["mjd_tt"])-50814.)*86400 for x in peaks]
    unused=set(range(len(actual)))
    matches=[]
    for j, lit in enumerate(RX_FITTED_MET):
        candidates=sorted(((abs(actual[k]-lit), k) for k in unused),key=lambda v:v[0])
        if candidates and candidates[0][0] <= max_abs_s:
            residual,k=candidates[0]
            unused.remove(k)
            matches.append({"published_peak_index":j+5,"published_met_s":lit,
                            "detected_met_s":round(actual[k],1),
                            "offset_s":round(actual[k]-lit,1)})
    return {"matched_of_10_published_fitted":len(matches),
            "matching_tolerance_s":max_abs_s,
            "individual_matches":matches,
            "published_total_2020_fitted_or_partial":10,
            "published_bright_segments":8}

def main():
    tests=verify_synthetic()
    rows=[]
    total_detected=0
    recurring=0
    n_curves=0
    for name in EXPECTED:
        spec=TARGETS[name]
        filename=f"control_curve_{name.replace(' ', '_').replace('+','p')}_{spec['obsid']}.csv"
        path=OUT/filename
        if not path.is_file():
            rows.append({"source":name,"obsid":spec["obsid"],
                         "expected_published_bright_eruptions":EXPECTED[name],
                         "status":"NOT_RETRIEVED","n_peaks":None,"recurrence_range_matched":False})
            continue
        curve=from_csv(path)
        result=find_eruptions(curve,50814.0,spec["published_h"])
        found=result.get("n_peaks",0)
        matched=bool(result.get("range_match"))
        total_detected+=min(found,EXPECTED[name])
        n_curves+=1
        recurring+=int(matched)
        row={"source":name,"obsid":spec["obsid"],
             "expected_published_bright_eruptions":EXPECTED[name],
             "observed_threshold_peaks":found,
             "max_count_overlap_without_one_to_one_matching":min(found,EXPECTED[name]),
             "interpeak_hours":result.get("intervals_h",[]),
             "recurrence_range_matched":matched,
             "status":result["status"],"curve_sha256":checksum(path),
             "peaks":result.get("peaks",[]),
             "published_period_h_range":spec["published_h"],
             "net_baseline_cps":result.get("net_baseline_cps"),
             "background_baseline_cps":result.get("background_baseline_cps"),
             "binned_good_exposure_h":result.get("binned_good_exposure_h"),
             "terminal_GTI_note":"Saved CSV omits invalid bins; replay pads 10 masked bins at end; late events may be partially observed"}
        if name=="RX J1301.9+2747":
            row["rx_individual_published_peak_alignment"]=rx_individual_match(
                result.get("peaks",[]))
        if name=="eRO-QPE4" and len(result.get("intervals_h",[]))==2:
            row["published_intervals_h"]=[9.82,14.70]
            row["interval_residual_minutes"]=[
                round(60*(x-y),2) for x,y in zip(result["intervals_h"],[9.82,14.70])]
        if name=="eRO-QPE3" and result.get("intervals_h"):
            row["published_first_two_interval_h_approx"]=20.4
            row["first_interval_residual_minutes"]=round(
                60*(result["intervals_h"][0]-20.4),2)
        rows.append(row)
    n_expected=sum(EXPECTED.values())
    # Note: count overlap is not flare-level sensitivity until every published
    # peak has a one-to-one time-aligned match and background exclusions pass.
    summary={"publication":LIT,"method":"NET PPS EPIC-pn 0.2-2keV RATE1-3, 300s; 4-sigma prominence, 3.5-sigma per bin, 7-sigma integrated, 3x quiet baseline; scaled-background veto",
             "method_status":"positive-control tuned; NOT independently surveyed; terminal masked GTI padding applied in replay",
             "synthetic_fixture":tests,
             "total_positive_control_epochs":len(EXPECTED),
             "epochs_with_usable_position_verified_pn_curves":n_curves,
             "epochs_with_recurrent_peaks_in_published_range":recurring,
             "operational_recurrence_recovery_fraction_all_controls":round(recurring/len(EXPECTED),5),
             "operational_recurrence_recovery_fraction_usable_curves":round(recurring/n_curves,5) if n_curves else None,
             "published_bright_eruption_count_target":n_expected,
             "detected_peaks_capped_to_published_counts":total_detected,
             "eruption_count_overlap_upper_bound":round(total_detected/n_expected,5),
             "one_to_one_published_eruption_recovery_fraction":"UNDEFINED - individual published peak times not all independently matched",
             "archive_screen_permitted":False,
             "new_candidates_examined_after_validation":0,
             "validated_new_qpes":0,
             "scientific_blockers":[
                 "Exact individual published peak times not tested for all 29 known eruptions",
                 "PPS current source extraction position verified but independent full field background FBKTSR light curves untested",
                 "Particle flare GTIs, EPIC pile-up, and 0.2-2keV spectral-temperature/hardness evolution not independently validated",
                 "No SDSS/Legacy/Gaia galaxy-nucleus and stellar sample matched; no measured selection completeness or false-positive rate",
                 "No valid host-nuclear archive screen performed, therefore no candidate-specific TNS/ATel/Swift/Chandra novelty audit",
                 "Tuning on known controls does not constitute an independent validation set"
             ]}
    store("control_replay.json",{"summary":summary,"controls":rows})
    store("control_gate_replay.json",summary)
    print(json.dumps({"summary":summary,
                      "controls":[{"name":r["source"],"n":r.get("observed_threshold_peaks"),
                                   "recurrence":r["recurrence_range_matched"]}
                                  for r in rows]},indent=2))
if __name__=="__main__":
    main()
