"""Join every persistent detection to full-camera measurements and record each check."""
import argparse
import gzip
import json
from pathlib import Path

import fitsio

from vetting import LAMP, POST_INSPECTION, REST, SKY, background_hypotheses


def near(table, wave, tol):
    return [k for k, v in table.items() if abs(v - wave) <= tol]


def adjudicate(scans, controls, cache):
    frames = {}
    if controls.is_file():
        with gzip.open(controls, "rt") as stream:
            records = [json.loads(line) for line in stream]
    else:
        records = [json.loads(file.read_text()) for file in controls.glob("cframe-*.json")]
    for record in records:
        name = record["source"]["url"].split("/")[-1].replace(".fits.gz", ".json")
        for m in record["measurements"]:
            key = m["targetid"], m["expid"], m["arm"], m["barycentric_A"]
            if key in frames:
                raise ValueError(f"duplicated cframe measurement {key}")
            frames[key] = m, name
    for scan in sorted(scans.glob("scan-*.jsonl")):
        tag = scan.stem.removeprefix("scan-")
        with fitsio.FITS(str(cache / f"rvtab_spectra-{tag}.fits")) as rvtab:
            tab = rvtab["RVTAB"].read(columns=["TARGETID","EXPID","VRAD_ERR"])
            errs = {(int(r["TARGETID"]),int(r["EXPID"])): float(r["VRAD_ERR"]) for r in tab}
        for line in scan.read_text().splitlines():
            s = json.loads(line)
            meta, arms = s["meta"], s["arms"]
            tid = meta["targetid"]
            for arm, a in arms.items():
                for candidate in a["candidates"]:
                    measurements = []
                    for h in candidate["hits"]:
                        exp, obs, rest = h["expid"], h["observed_A"], h["rest_A"]
                        m, file = frames[(tid, exp, arm, obs)]
                        dob = obs/rest
                        tol = (h["rest_fwhm_A"]**2 +
                               (3 * errs[tid, exp] * rest / 299792.458)**2)**0.5
                        # All 7-sigma positive shape-compatible peaks in this exposure,
                        # across the three cameras. Report any possible coincidence.
                        other = [x["observed_A"] for band in arms.values()
                                 for x in band["exposure_hits"][str(exp)] if x["shape_pass"]]
                        intervals = [(lo*dob, hi*dob) for band in arms.values()
                                     for lo,hi in band["coverage"]]
                        hypotheses = background_hypotheses(obs, m["lsf_fwhm_A"], other, intervals)
                        projected = []
                        for hypo in hypotheses:
                            companions = hypo["companions"]
                            projected.append({"anchor": hypo["anchor"], "z": hypo["z"],
                                              "accessible_companions": len(companions),
                                              "matched_companions": {name: v for name,v in companions.items()
                                                                     if v["peaks_within_one_fwhm_A"]}})
                        topo = m["topocentric_A"]
                        twidth = m["lsf_fwhm_A"] / m["topo_to_bary"]
                        sky = m["sky_controls"]
                        high_sky = [x["fiber"] for x in sky if abs(x["snr"]) >= 4]
                        weak_sky = [x["fiber"] for x in sky if abs(x["flux"])+3*x["sigma"] >= h["flux17"]]
                        traces = m["trace_controls"]
                        measurements.append(dict(
                            expid=exp, fiber=m["fiber"], controls=file,
                            barycentric_A=obs, stellar_rest_A=rest, topocentric_A=topo,
                            residual_flux17=h["flux17"], residual_snr=h["snr"], lsf_shape_p=h["shape_p"],
                            raw_target_snr=m["target_raw"].get("snr"),
                            reference_rest_hits=near(REST,rest,tol),
                            lithium_hits=near(POST_INSPECTION,rest,tol),
                            min_distance_to_li_A=min(abs(w-rest) for w in POST_INSPECTION.values()),
                            topocentric_sky_hits_A=[w for w in SKY if abs(w-topo)<=twidth],
                            topocentric_lamp_hits=near(LAMP,topo,twidth),
                            n_sky=len(sky), sky_max_abs_snr=max((abs(x["snr"]) for x in sky),default=None),
                            sky_shared_4sigma=high_sky, sky_insensitive_3sigma=weak_sky,
                            unavailable_trace_fibers=[x["fiber"] for x in traces if "sigma" not in x],
                            trace_max_abs_snr=max((abs(x["snr"]) for x in traces if "sigma" in x),default=None),
                            projected_background_redshifts=projected))
                    epochs = {x["expid"] for x in meta["exposures"]}
                    if epochs != {x["expid"] for x in measurements} or len(measurements) != len(epochs):
                        raise ValueError("missing or duplicated candidate exposure")
                    li = all(x["lithium_hits"] for x in measurements)
                    yield dict(targetid=tid, healpix=meta["healpix"], arm=arm,
                               survey=meta["survey"], rest_A=candidate["rest_A"], exposures=measurements,
                               status="excluded_stellar_LiI_window" if li else "needs_full_review",
                               natural_line_source="NIST Li I resonance doublet" if li else None,
                               unresolved_checks=[
                                   "Raw 2D CCD pixel/trace inspection and bright-neighbor transfer-function calibration",
                                   "Sky-model and telluric transmission inspection beyond line list and local sky fibres",
                                   "Full nebular/flaring/Be/WR/symbiotic classification, subthreshold companion-line fitting and imaging",
                                   "Cross-program exposure inventory and independent new observations"
                               ])


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ("scans", "controls", "cache", "output"):
        p.add_argument(key, type=Path)
    args = p.parse_args()
    rows = list(adjudicate(args.scans,args.controls,args.cache))
    args.output.write_text("".join(json.dumps(r, separators=(",",":"), allow_nan=False)+"\n" for r in rows))
    print(len(rows), "candidates;", sum(len(r["exposures"]) for r in rows), "candidate-exposures")
