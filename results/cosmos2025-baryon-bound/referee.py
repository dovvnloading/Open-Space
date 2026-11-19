#!/usr/bin/env python3
"""Independent arithmetic/data-integrity referee for the COSMOS2025 audit.

Does *not* certify astrophysical candidate identification or a LCDM null.
Returns a machine-readable list of blockers, even if arithmetic passes.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad


def independent_volume(zlo, zhi, Om, H0, area):
    c = 299792.458
    h = lambda z: H0*math.sqrt(Om*(1+z)**3+(1-Om))
    d = lambda z: c*quad(lambda zp: 1/h(zp), 0,z,epsrel=1e-10)[0]
    # Exact shell formula as an independent control of the trapezoid-in-z implementation.
    return area*(math.pi/180)**2*(d(zhi)**3-d(zlo)**3)/3


def review(folder):
    blockers = []
    verified = []
    raw = folder/"run"/"receipt.json"
    hmf = folder/"run"/"halo_bounds.json"
    if not raw.exists():
        blockers.append("Catalogue scan receipt not available.")
    else:
        data = json.loads(raw.read_text())
        if data.get("status","").startswith("FAILED"):
            blockers.append("Official CDS source scan failed: " + data.get("error","unknown")[-500:])
        else:
            c=data.get("catalogue",{})
            assert c.get("row_count")==784016, "NOT a full 784016-row catalogue"
            for f in ("readme","lephare"):
                rec=data[f]
                assert len(rec["sha256"])==64 and all(x in "0123456789abcdef" for x in rec["sha256"])
                assert rec["size_bytes"]>100
                assert "cds" in rec["url"].lower()
            verified.append("Whole-source CDS record count and hashed byte receipts.")
            blockers.append("CDS mirror version not demonstrated identical to IAP v1.1 May 2026.")
            blockers.append("zPDF-u68 / marginal mass quantiles are not joint p(M*,z); cannot establish true high-z candidate denominator.")
            blockers.append("Catalogue rows lack individual source IDs/coordinates until photometry product is paired and validated.")
            full=folder/"run"/"highz_triage_all.csv.gz"
            if not full.exists():
                blockers.append("Full marginal high-z triage is not deposited; top 250 may hide candidates.")
            else:
                assert hashlib.sha256(full.read_bytes()).hexdigest()==c["triage_all_sha256"]
                with gzip.open(full,"rt",newline="") as ff:
                    allrows=list(csv.DictReader(ff))
                assert len(allrows)==c["galaxy_zfinal_or_u68_gt6"]
                verified.append("ALL marginal high-z source rows retained and compressed with an exact SHA256.")

            top=folder/"run"/"highz_triage_top250.csv"
            if top.exists():
                with top.open() as f:rows=list(csv.DictReader(f))
                assert len(rows)<=250 and len(rows)>0
                masses=[float(x["mass_u68"]) for x in rows if x["mass_u68"]]
                assert all(a>=b for a,b in zip(masses,masses[1:])), "triage not mass ranked"
                verified.append("Ranked 68%-mass-upper triage rows are monotonic.")
    phot=folder/"run"/"photometry_receipt.json"
    if not phot.exists():
        blockers.append("No hashed full 37-band photometry source and row-index/coordinate association.")
    else:
        assoc=json.loads(phot.read_text())
        assert assoc["record_count"]==784016
        assert len(assoc["sha256_raw_download"])==64
        joined=folder/"run"/"highz_candidate_photometry.csv"
        assert joined.exists()
        with joined.open(newline="") as pf:
            sources=list(csv.DictReader(pf))
        assert len(sources)==assoc["joined_shortlist_rows"]
        assert len(set(row["Id"] for row in sources))==len(sources)
        verified.append("Full photometry source hashed and top-ranked candidate IDs/coordinates unique in CDS row alignment.")
        blockers.append("CDS cross-extension row alignment is not an independent cross-release v1.1 match; no imaging or SED validation.")

    if not hmf.exists():
        blockers.append("Numerical ST/Tinker extreme-halo results not deposited.")
    else:
        a=json.loads(hmf.read_text())
        assert math.isclose(a["f_baryon"],0.16,abs_tol=1e-12)
        grid=np.array(a["halo_mass_grid_msun"])
        assert len(grid)>100 and np.all(np.diff(grid)>0)
        for b in a["bins"]:
            lo,hi=b["z_min_exclusive"],b["z_max_inclusive"]
            vol=independent_volume(lo,hi,a["cosmology"]["Om0"],a["cosmology"]["H0"],a["area_deg2_nominal"])
            assert abs(vol/b["volume_mpc3"]-1)<1e-4,(lo,hi,vol,b["volume_mpc3"])
            for model in ("sheth99","tinker08"):
                entry=b["models"][model]
                lam=np.array(entry["lambda_expected_count_by_halo_mass_grid"],float)
                assert len(lam)==len(grid) and np.all(np.isfinite(lam))
                assert np.all(lam>=0) and np.all(np.diff(lam)<=1e-5)
                assert (entry["mass_definition"]=="fof") if model=="sheth99" else (entry["mass_definition"]=="200m")
                quant=entry["max_halo_quantiles_msun"]
                vals=[quant["p500"],quant["p950"],quant["p990"],quant["p999"]]
                assert all(v is not None for v in vals)
                assert all(x<y for x,y in zip(vals,vals[1:])),(lo,hi,model,quant)
                for key in ("p500","p950","p990","p999"):
                    assert abs(entry["max_stellar_baryon_quantiles_log10"][key]-
                               math.log10(0.16*quant[key]))<1e-8
                count_at_expected=float(np.interp(math.log10(entry["M_halo_msun_expected_one"]),
                                                   np.log10(grid),lam))
                assert abs(count_at_expected-1)<0.1,(lo,hi,model,count_at_expected)
        verified.append("Sheth99/Tinker08 halo mass definitions, monotone extreme-value CDFs, baryon fractions and independently integrated shell volumes.")
        blockers.append("Halo extremes still use nominal unmasked 0.54 deg² and Poisson, no survey completeness/cosmic variance.")
    blockers += [
        "No vetted joint z–mass posterior chains from two independently fitted SED codes.",
        "No forward Eddington/selection correction, trial accounting or population posterior predictive test.",
        "No completed source-level MIRI, brown-dwarf, dusty low-z, nebular-line, AGN, blend, X-ray/radio and NIRSpec/DJA/MAST checks."
    ]
    return {"referee_disposition": "BLOCKED — do not claim detection, measured null, or [REMARKABLE]",
            "scientific_discovery_ready":False,
            "numerical_checks_passed":verified,
            "required_validation_missing":blockers}


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=Path(__file__).resolve().parent)
    p.add_argument("--out",type=Path,default=None)
    a=p.parse_args()
    report=review(a.root)
    out=a.out or a.root/"referee_report.json"
    out.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=="__main__":
    main()
