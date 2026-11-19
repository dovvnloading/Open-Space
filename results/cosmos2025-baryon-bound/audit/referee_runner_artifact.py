#!/usr/bin/env python3
"""Independent audit of a downloaded COSMOS2025 GitHub runner ZIP.

Usage: python referee_runner_artifact.py --zip cosmos2025-onejob-noncds-receipts.zip --out audit.json
Requires numpy and scipy; no web access, no CDS, no trust in deposited referee.json.
"""
import argparse
import hashlib
import json
import math
import zipfile
from pathlib import Path
import numpy as np
from scipy.integrate import quad

def sha(raw):
    return hashlib.sha256(raw).hexdigest()

def review(path):
    with zipfile.ZipFile(path) as z:
        files = {n: z.read(n) for n in z.namelist()}
    source = json.loads(files["mirror_receipt.json"])
    halo = json.loads(files["halo_bounds.json"])
    old = json.loads(files["independent_referee.json"])
    assert source["status"] == "BLOCKED_NO_VERIFIED_NONCDS_CATALOGUE"
    assert source["source_sha256"] is None and source["source_bytes_verified"] == 0
    assert len(source["checks"]) == 9
    assert not source["candidate_urls"] and not source["download_attempts"]
    hosts = {c["id"]: c for c in source["checks"]}
    assert hosts["zenodo_cosmos2025"]["http"]["status"] == 200
    assert hosts["zenodo_cosmos2025"]["zenodo_records_examined"] == 8
    assert hosts["zenodo_mastercatalog"]["zenodo_records_examined"] == 0
    assert hosts["irsa_tap_tables"]["mentions_cosmos2025"] is False
    assert source["derived_sha256"]["halo_bounds.json"] != sha(files["halo_bounds.json"])
    assert old["source_receipt_sha256"] != sha(files["mirror_receipt.json"])
    assert old["halo_receipt_sha256"] != sha(files["halo_bounds.json"])
    assert halo["f_baryon"] == .16 and len(halo["bins"]) == 6
    H0, OM = halo["cosmology"]["H0"], halo["cosmology"]["Om0"]
    area = halo["area_deg2_nominal"] * (math.pi / 180.)**2
    c = 299792.458
    integrand = lambda z: 1. / (H0 * math.sqrt(OM*(1+z)**3+1-OM))
    Dc = lambda z: c * quad(integrand, 0, z, epsrel=1e-10)[0]
    masses = np.asarray(halo["halo_mass_grid_msun"], float)
    volerr, tailerr = 0., 0.
    for shell in halo["bins"]:
        exact = area / 3 * (Dc(shell["z_max_inclusive"])**3 - Dc(shell["z_min_exclusive"])**3)
        volerr = max(volerr, abs(exact/shell["volume_mpc3"]-1))
        for model, massdef in (("sheth99", "fof"), ("tinker08", "200m")):
            entry = shell["models"][model]
            assert entry["mass_definition"] == massdef
            lam = np.asarray(entry["lambda_expected_count_by_halo_mass_grid"], float)
            assert np.all(lam >= 0) and np.all(np.diff(lam) <= 1e-6)
            q = entry["max_halo_quantiles_msun"]
            assert [q[k] for k in ("p500", "p950", "p990", "p999")] == sorted(q.values())
            for key, cdf in (("p500", .5), ("p950", .95), ("p990", .99), ("p999", .999)):
                lnlam = np.interp(math.log(q[key]), np.log(masses), np.log(np.maximum(lam, 1e-250)))
                tailerr = max(tailerr, abs(math.exp(lnlam)/(-math.log(cdf))-1))
                assert abs(entry["max_stellar_baryon_quantiles_log10"][key] -
                           math.log10(.16*q[key])) < 1e-8
    assert volerr < 1e-4 and tailerr < .02
    return {
        "review": "LOCAL ARTIFACT REPLAY; not external astronomical peer review",
        "run_id": source["github_run_id"],
        "zip_sha256": sha(path.read_bytes()),
        "file_sha256": {name:sha(data) for name,data in files.items()},
        "source_bytes": 0, "source_row_count": None,
        "zenodo_records_examined": 8,
        "HMF_six_shells_checked": True,
        "max_volume_relative_error": volerr,
        "max_Poisson_tail_relative_error": tailerr,
        "second_run_machine_referee": "STALE/FAILED: halo hash refers to first-run artifact",
        "scientific_status": "INCONCLUSIVE_NO_CATALOGUE_NOT_NULL",
        "ready_to_merge": False, "remarkable": False
    }

def main():
    a=argparse.ArgumentParser()
    a.add_argument("--zip", type=Path, required=True)
    a.add_argument("--out", type=Path, required=True)
    args=a.parse_args()
    report=review(args.zip)
    args.out.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps(report, indent=2))

if __name__ == "__main__":
    main()
