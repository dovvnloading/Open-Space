"""Verify the measurements from the released compact tables (Python 3.10+)."""
import csv
import json
import math
from pathlib import Path

DATA = Path(__file__).resolve().parent / "data"


def rows(name):
    with (DATA / name).open() as f:
        return list(csv.DictReader(f))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    all_events = rows("events.csv")
    matched = rows("matched_controls.csv")
    for obs, expected, slow, removed, fast_controls in [
        ("1020240141", 13, 7, 5, 131),
        ("5505010502", 34, 32, 0, 555),
    ]:
        burst = [r for r in all_events if r["obsid"] == obs and abs(float(r["delta_us"])) < 10]
        cleaned = [r for r in burst if r["cleaned"]=="1" and 50<=int(r["pi"])<=800]
        rejected = [r for r in burst if r["cleaned"]=="0"]
        require(len(burst)==expected, "Unexpected unfiltered multiplicity")
        require(len(cleaned)==slow and all(int(r["flags"])==16 for r in cleaned), "Clean event flags differ")
        require(len(rejected)==removed and all(int(r["flags"])==24 and
            float(r["pi_ratio"])>1.1+120/int(r["pi"]) for r in rejected), "Rejected pulse-height ratios differ")
        require(len({r["det_id"] for r in burst})==len(burst), "Detector repeated")
        width = max(float(r["delta_us"]) for r in burst) - min(float(r["delta_us"]) for r in burst)
        require(min(float(r["deadtime_us"]) for r in burst)>width, "Deadtime check failed")
        ctrl = [r for r in matched if r["obsid"]==obs]
        require(sum(int(r["fast"]) for r in ctrl)==fast_controls and
                sum(int(r["controls"]) for r in ctrl)==20*slow, "Detector-/energy-matched controls differ")
        print(f"{obs}  {len(cleaned)} cleaned 0.5–8 keV events, all slow only; "
              f"{removed} particle-filter events; {fast_controls}/{20*slow} matched fast triggers")

    attitude = rows("pointing.csv")
    for r in attitude:
        ra, dec, source_ra, source_dec = [math.radians(float(r[k])) for k in ("ra","dec","source_ra","source_dec")]
        hav = math.sin((source_dec-dec)/2)**2 + math.cos(source_dec)*math.cos(dec)*math.sin((source_ra-ra)/2)**2
        sep = 120*math.degrees(math.asin(math.sqrt(hav)))
        require(abs(sep-float(r["offset_arcmin"]))<0.001, "Coordinate separation differs")
        if r["catalog_source"].startswith("CXOU"):
            require(r["fits_target"]=="PSR_J0058-7218" and float(r["gti_min"])>12, "SMC pointing differs")
            print(f'{r["obsid"]}  {sep:.2f} arcmin from the SMC magnetar')
    checks = json.loads((DATA/"comparison.json").read_text())
    require(checks["burst_storm"]["hard_fast"]==188788 and checks["burst_storm"]["hard"]==188793,
            "Burst storm comparison differs")
    require(checks["vignetting"]["largest_ratio"]<0.00030, "Calibration map comparison differs")
    print("Packet checks passed.")


if __name__=="__main__":
    main()
