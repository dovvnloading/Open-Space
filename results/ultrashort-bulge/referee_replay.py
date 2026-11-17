"""Referee replay against exact public survey bytes, not merely deposited fit tables.

The downloaded KMT 2024 pySIS tarballs and OGLE 2026 photometry are
independently hashed, then the deposited PSPL/FSPL/binary parameter vectors
are evaluated on their *original* epochs and a new PSPL/FSPL/binary fit is
optimized for two decisive cases. The 2016 positive control is refitted from
deposited original photometry. Uses only 2024/2026 public archives.
"""
import argparse
import csv
import gzip
import hashlib
import json
import urllib.request
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from peak_models import solve
from earth_control import run as recover_control

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
TARGETS = ("KMT-2024-BLG-3326", "KMT-2024-BLG-0824", "KMT-2024-BLG-3278")
OGLE = "OGLE-2026-BLG-0670"


def download(rec, destination):
    req = urllib.request.Request(rec["url"], headers={"User-Agent": "Open-Space-referee-reproducibility/1.0"})
    with urllib.request.urlopen(req, timeout=90) as resp:
        raw = resp.read()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != rec["sha256"]:
        raise AssertionError(f"Archive SHA-256 differs for {destination.name}: {actual} != {rec['sha256']}")
    if len(raw) != int(rec["size"]):
        raise AssertionError(f"Archive size differs for {destination.name}")
    destination.write_bytes(raw)
    return dict(input_sha256=actual, size=len(raw), source_url=rec["url"])


def fixture_models(event, all_peaks):
    if event == OGLE:
        return {"referee-ogle-binned": json.loads((HERE/"other/ogle0670-models.json").read_text())}
    return {name:record for name,record in all_peaks[event].items() if record}


def seed_from(record):
    m = record["models"]
    ps,fs,b = m["ps"],m["fs"],m["bin"]
    return {"ps_t0":ps["t0"],"ps_te":ps["te_or_rise_d"],"ps_u0":ps["u0_or_fall_d"],
            "fs_t0":fs["t0"],"fs_te":fs["te_or_rise_d"],"fs_u0":fs["u0_or_fall_d"],
            "fs_rho":fs["rho"],"bin_s":b["s"],"bin_q":b["q"],"bin_angle":b["alpha"],
            "bin_te":b["te_or_rise_d"]}


def run(refit=False):
    from validate import run as validate
    summary = validate(models=True)
    manifest = {r["event"]: r for r in csv.DictReader((HERE/"kmt/pySIS-2024.tsv").open(),delimiter="\t")}
    with (HERE/"ogle/photometry.csv").open() as f:
        ogle = next(r for r in csv.DictReader(f) if r["event"]=="2026-BLG-0670")
    ogle = {"url": ogle["url"],"sha256": ogle["sha256"],"size": ogle["bytes"]}
    models = json.loads(gzip.decompress((HERE/"tables/peak-comparisons.json.gz").read_bytes()))
    reports = []
    with TemporaryDirectory() as tmp:
        root = Path(tmp)
        for event in (*TARGETS, OGLE):
            if event == OGLE:
                path = root/"blg-0670.dat"
                receipt=download(ogle,path)
            else:
                path = root/(event+".tar.gz")
                receipt=download(manifest[event],path)
            from run_grid import read_curve
            job=(event,path,0.,0.)
            samples = len(read_curve(job)[0])
            comparisons = []
            for exp, rec in fixture_models(event, models).items():
                kw={"window":rec.get("window_halfwidth_d",4.),"bin_width":rec.get("bin_width_d",0.)}
                # The historical peak records without explicit fields used +/-4d
                scores=solve(job,{"fs_t0":rec["reference_jd"]},rec["cap_factor"],
                             evaluate=rec["models"],**kw)
                diffs={name:round(abs(v-rec["models"][name]["chi2"]),6) for name,v in scores.items()}
                assert max(diffs.values()) <= .05,(event,exp,diffs)
                comparison=dict(experiment=exp,models=len(diffs),
                                max_abs_chi2=max(diffs.values()),chi2_by_model=scores,
                                deposited_binary_converged=rec["models"].get("bin",{}).get("converged"))
                if refit and exp in ("initial", "referee-ogle-binned") and {"ps","fs","bin"}.issubset(rec["models"]):
                    # Preserve the exact original comparison window: solve() takes
                    # the epoch origin from fs_t0 in its seed. Moving this origin
                    # changes which epochs enter the objective and makes the
                    # absolute chi2 values scientifically non-comparable.
                    fixed_seed=seed_from(rec)
                    fixed_seed["fs_t0"]=rec["reference_jd"]
                    fresh=solve(job,fixed_seed,rec["cap_factor"],binary_trials=12,
                                profiles=True,**kw)
                    f=fresh["models"]
                    comparison["fresh_optimization"]={
                        name:{"chi2":item["chi2"],"te_d":item.get("te_or_rise_d",item.get("te_d")),
                              "rho":item.get("rho"),"converged":item.get("converged")}
                        for name,item in f.items()}
                    # Multiple likelihood minima are possible: do not silently
                    # equate a local optimizer endpoint with global evidence.
                comparisons.append(comparison)
            reports.append(dict(event=event,valid_input=receipt,native_rows=samples,experiments=comparisons))
        raw=gzip.decompress((HERE/"controls/earth-figure1.tsv.gz").read_bytes())
        stored=json.loads((HERE/"controls/earth-recovery.json").read_text())
        assert hashlib.sha256(raw).hexdigest()==stored["response_sha256"]
        replay=recover_control(raw)
        assert abs(replay["finite_no_blend"]["te_d"]-stored["finite_no_blend"]["te_d"])<.00015
        assert abs(replay["finite_no_blend"]["chi2"]-stored["finite_no_blend"]["chi2"])<.1
        control=dict(input_sha256=stored["response_sha256"],rows=replay["input_rows"],
                     te_min=replay["finite_no_blend"]["te_d"]*1440,
                     rho=replay["finite_no_blend"]["rho"],
                     finite_minus_point_chi2=replay["finite_no_blend"]["chi2"]-replay["point_no_blend"]["chi2"])
    result={"result":"PASS","offline_validation":summary,"native_input_replays":reports,
            "positive_control":control,"refit_enabled":refit}
    print(json.dumps(result,indent=2))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--refit",action="store_true",help="also repeat nonlinear PSPL/FSPL/binary optimization for all four priority events")
    args=p.parse_args()
    run(refit=args.refit)
