#!/usr/bin/env python3
"""Recover published Xing et al. J1010 Table 1 without misusing limits."""
import json
from pathlib import Path
here = Path(__file__).resolve().parent
data = json.loads((here / "j1010_xing_table1.json").read_text())
a = data["original"]
required = ("na", "mg", "al", "co", "zn", "ba", "sr")
measured = [x for x in required if isinstance(a.get(x), dict) and "value" in a[x]]
limited = [x for x in required if isinstance(a.get(x), dict) and "upper" in a[x]]
missing = [x for x in required if a.get(x) is None]
na_mg_upper = round(a["na"]["upper"] - a["mg"]["value"], 2)
assert na_mg_upper == -1.36
assert a["co"]["value"] == -0.72
assert a["zn"]["upper"] == -0.25
assert a["ba"]["upper"] == -1.37
assert a["sr"]["upper"] == -2.25
assert a["mg"]["sigma"] == 0.12
assert "al" in missing and "zn" in limited
result = {
    "control": data["object"], "FeH_dex": data["FeH"],
    "original_data": "Xing et al. 2023 Table 1; spectroscopic literature values, not GALAH-derived",
    "na_minus_mg_central_upper_dex": na_mg_upper,
    "measured_required": measured, "censored_upper_required": limited,
    "not_measured_required": missing,
    "eligible_for_GALAH_five_element_detection_only_gate": False,
    "limitation": "An upper limit is not a detected depletion; later high-resolution analyses challenged the pure-PISN interpretation.",
    "validated_against_published_table": True,
    "independent_reduction": False,
}
out = here / "j1010_control_receipt.json"
out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
print("J1010 original Table 1 recovered: Na/Mg < -1.36, Co/Fe=-0.72; Al absent, Zn upper; NOT a qualifying five-element detection.")
