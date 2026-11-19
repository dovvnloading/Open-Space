#!/usr/bin/env python3
"""Independent, failure-closed review of non-CDS COSMOS2025 acquisition/HMF outputs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
from referee import independent_volume

def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def audit(run, halo):
    blockers, confirmed = [], []
    recfile = run / "mirror_receipt.json"
    if not recfile.exists():
        blockers.append("Exact mirror source receipt is absent.")
        source = {}
    else:
        source = json.loads(recfile.read_text())
        for entry in source.get("checks", []):
            assert entry["url"].startswith("https://")
            assert "http" in entry
        confirmed.append("Frozen source and alternative-host DNS/HTTP responses recorded.")
        for filename, expected in source.get("derived_sha256", {}).items():
            path = run / filename
            assert path.is_file() and sha(path) == expected, (filename, "bad output hash")
        confirmed.append("Derived files match their SHA-256 receipts.")
        if source.get("source_sha256") is None:
            assert source.get("source_bytes_verified", 0) == 0
            blockers.append("No non-CDS official COSMOS2025 784016-row FITS acquired: cannot count z>6 sources.")
        else:
            assert len(source["source_sha256"]) == 64
            assert source["source_bytes_verified"] >= 100000
            assert source.get("source_url", "").startswith("https://")
            confirmed.append("Source-byte hash and FITS extension/row-validation evidence deposited.")
            if source["status"] != "CATALOGUE_VERIFIED_BUT_NO_JOINT_POSTERIOR_OR_REFITS":
                blockers.append("Source not screened into even preliminary triage.")
            else:
                triage = json.loads((run / "fits_triage.json").read_text())
                assert triage["source_rows"] == 784016
                assert triage["all_marginal_highz_rows"] >= len(triage["top250"])
                confirmed.append("All FITS rows scanned in chunks; high-z ranking remains marginal triage only.")
    if halo.is_file():
        h = json.loads(halo.read_text())
        assert math.isclose(h["f_baryon"], .16, abs_tol=1e-12)
        assert len(h["bins"]) == 6
        for shell in h["bins"]:
            exact = independent_volume(shell["z_min_exclusive"], shell["z_max_inclusive"],
                                       h["cosmology"]["Om0"], h["cosmology"]["H0"],
                                       h["area_deg2_nominal"])
            assert abs(exact / shell["volume_mpc3"] - 1) < 1e-4
        confirmed.append("New six-shell HMF has independent comoving-volume validation.")
        blockers.append("HMF nominal unmasked field and distinct halo definitions cannot imply a strict cosmological bound.")
    else:
        blockers.append("No newly run two-HMF baseline.")
    blockers.extend([
        "No independent joint p(Mstar,z) source posterior chain; LePHARE fixed-z mass percentiles inadmissible.",
        "No two-code EAZY and Prospector/Bagpipes AGN, lines, brown-dwarf and low-z interloper fits of every candidate.",
        "No per-candidate MIRI, X-ray/radio, blends, DJA/MAST NIRSpec or literature crosschecks.",
        "No Eddington/selection, cosmic variance and trials-aware population posterior predictive analysis."
    ])
    return {"referee_type": "independent method/code-path audit, not external scientific replication",
            "disposition": "BLOCKED — NOT a measured null or violation; keep PR draft, do not merge",
            "ready_to_merge": False, "remarkable": False,
            "source_status": source.get("status", "no receipt"),
            "source_receipt_sha256": sha(recfile) if recfile.exists() else None,
            "halo_receipt_sha256": sha(halo) if halo.is_file() else None,
            "checks_confirmed": confirmed, "blocking_defects": blockers}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--halo", type=Path, required=True)
    a = p.parse_args()
    a.run.mkdir(parents=True, exist_ok=True)
    report = audit(a.run, a.halo)
    (a.run / "independent_referee.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"status": report["disposition"], "source": report["source_status"],
                      "confirmed": len(report["checks_confirmed"]), "blockers": len(report["blocking_defects"])}, indent=2))

if __name__ == "__main__":
    main()
