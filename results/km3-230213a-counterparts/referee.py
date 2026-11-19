#!/usr/bin/env python3
"""Independent readback: verify generated evidence, source counts, and null limitations.

Separate source implementation from analysis.py. Does not certify raw survey
completeness or native FRB localizations. Exits nonzero on contradictory claims.
"""
import csv
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parent
with (ROOT/"known_17.csv").open() as f: known=list(csv.DictReader(f))
with (ROOT/"known_replay.csv").open() as f: replay=list(csv.DictReader(f))
with (ROOT/"new_associations.csv").open() as f: new=list(csv.DictReader(f))
stat=json.loads((ROOT/"control_statistics.json").read_text())
receipts=json.loads((ROOT/"retrieval_receipts.json").read_text())
assert len(known)==17 and len(replay)==17
assert {r["known_id"] for r in replay}=={r["id"] for r in known}
assert all(r["status"]=="KNOWN_ALREADY_PUBLISHED" for r in replay)
assert all(float(r["distance_deg"])<=3.0 for r in replay)
assert all(r["p_chance"]=="NOT_ESTIMABLE_SELECTED_SAMPLE" for r in replay)
assert all(r["status"]=="UNVERIFIED_CATALOGUE_COINCIDENCE" for r in new)
assert all(r["p_physical_association"]=="NOT_ESTIMABLE_NO_EXPOSURE_FUNCTION" for r in new)
assert all(0<float(r["p_conditional_skyshift"])<=1. for r in new)
assert all(float(r["p_upper_bound_16_trials"])>=float(r["p_conditional_skyshift"]) for r in new)
assert all(float(r["sep_deg"])<=3.0 for r in new)
for family, s in stat["families"].items():
    assert s["n_shifts"]==999
    assert s["observed_in_3deg"]==sum(r["family"]==family for r in new)
    assert s["temporal_parent_count"]>=s["observed_in_3deg"]
    assert 0<s["p_shift_count_ge_observed"]<=1.
for name, val in stat["public_catalogue_status"].items():
    if val["status"] not in ("DOWNLOADED_FULL_PUBLIC_SNAPSHOT","DOWNLOADED_CATALOG_2"):
        assert name not in ("TNS","CHIME") or all(
            r["family"] not in (("TNS_OPTICAL","TNS_TDE_LONG") if name=="TNS" else ("CHIME_FRB",))
            for r in new)
# Independently recompute each published count from the frozen, time-filtered
# source positions with a vector-dot-product sphere formula (not the discovery
# implementation's haversine function). Recompute the entire 999-sky null.
import random
with (ROOT/"temporal_parent_sample.csv").open() as f:
    parents=list(csv.DictReader(f))
assert len(parents)==sum(d["temporal_parent_count"] for d in stat["families"].values())
def dot_dist(ref_ra, ref_dec, row):
    a,b,c,d=map(math.radians,(ref_ra,ref_dec,float(row["ra"]),float(row["dec"])))
    x=math.sin(b)*math.sin(d)+math.cos(b)*math.cos(d)*math.cos(a-c)
    return math.degrees(math.acos(max(-1.,min(1.,x))))
shifter=random.Random(230213)
shifts=[shifter.uniform(6.5,353.5) for _ in range(999)]
for fam, fstat in stat["families"].items():
    entries=[r for r in parents if r["family"]==fam]
    assert len(entries)==fstat["temporal_parent_count"]
    obs=sum(dot_dist(94.3,-7.8,r)<=3 for r in entries)
    assert obs==fstat["observed_in_3deg"]
    counts=[sum(dot_dist((94.3+delta)%360,-7.8,r)<=3 for r in entries) for delta in shifts]
    assert min(counts)==fstat["skyshift_min_count"]
    assert max(counts)==fstat["skyshift_max_count"]
    assert counts.count(0)==fstat["skyshift_n_zero"]
    assert abs((1+sum(x<=obs for x in counts))/1000-fstat["p_shift_count_le_observed"]) < 5.1e-7
    assert abs((1+sum(x>=obs for x in counts))/1000-fstat["p_shift_count_ge_observed"]) < 5.1e-7
    assert abs(sum(counts)/len(counts)-fstat["empirical_mean_count"])<0.00001

assert receipts["controls"]["valid_shifts"]==999
out={"referee_status":"PASS_WITH_STRICT_SCOPE_LIMITATIONS","known_entries_replayed":len(replay),
     "catalogue_matches":len(new),"family_counts":{k:v["observed_in_3deg"] for k,v in stat["families"].items()},
     "files_checked":["known_17.csv","known_replay.csv","new_associations.csv","temporal_parent_sample.csv",
                      "control_statistics.json","retrieval_receipts.json"],
     "unvalidated":["catalogue completeness","survey sky-exposure/footprints",
                    "original optical onset/peak and forced photometry",
                    "CHIME native localization likelihood","GRB/X-ray/LAT/radio/MIR full searches",
                    "physical-model neutrino production","all-family trials calibration"],
     "discovery_validated":False,
     "interpretation":"No source identity or discovery-level p is established by this replay."}
(ROOT/"referee_results.json").write_text(json.dumps(out,indent=2,sort_keys=True)+"\n")
print(json.dumps(out,sort_keys=True))
