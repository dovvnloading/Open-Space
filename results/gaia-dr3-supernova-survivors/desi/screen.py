"""Screen every exact DR3 identifier in the 2026 Swan DESI DR1 spectral catalogue."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from zero_point import zpt

from fetch import CAT, HERE, digest


def screen(cache):
    receipt = json.loads((HERE/"receipt.json").read_text())
    if digest((cache/CAT).read_bytes()) != receipt["catalogue_sha256"]:
        raise ValueError("catalogue digest mismatch")
    phot = pd.read_csv(cache/CAT, dtype={"designation":str},low_memory=False)
    entries = phot[phot.designation.str.startswith("Gaia DR3 ")].copy()
    entries["source_id"] = entries.designation.str.split().str[-1]
    chunks = []
    for batch in receipt["batches"]:
        path = cache/f"gaia_{batch['batch']:03}.csv"
        if digest(path.read_bytes()) != batch["response_sha256"]:
            raise ValueError(f"response differs in {path}")
        chunks.append(pd.read_csv(path, dtype={"source_id":str}))
    frame = pd.concat(chunks).merge(entries, on="source_id", validate="one_to_one")
    if len(frame) != len(entries):
        raise ValueError("incomplete Gaia coverage")
    k = 4.74047
    pm = np.hypot(frame.pmra, frame.pmdec)
    valid = frame.astrometric_params_solved.isin([31,95]) & (frame.parallax_error>0) & pm.notna()
    zpt.load_tables()
    zp = np.full(len(frame), np.nan)
    zcols = ["phot_g_mean_mag","nu_eff_used_in_astrometry","pseudocolour","ecl_lat","astrometric_params_solved"]
    zp[valid] = zpt.get_zpt(*(frame.loc[valid,c].to_numpy() for c in zcols), _warnings=False)
    frame["zp_mas"] = zp
    par = frame.parallax - zp
    sigma = frame.parallax_error
    frame["plx_snr"] = par/sigma
    frame["vt_kms"] = k*pm/par.where(par>0)
    quality = ((frame.ruwe>0) & (frame.ruwe<1.4) & (frame.visibility_periods_used>=9)
               & frame.duplicated_source.eq(False) & (frame.ipd_frac_multi_peak<=2)
               & (frame.ipd_gof_harmonic_amplitude<0.1))
    # Keep all spectral classes; the catalogue also contains programme contaminants.
    frame["quality"] = quality
    controls = pd.read_csv(HERE.parent/"control"/"controls.csv",dtype={"source_id":str})
    fields = ["source_id","WDJname","specType","specType_confidence","parallax","parallax_error",
              "zp_mas","plx_snr","pmra","pmdec","vt_kms","quality","ruwe","ipd_frac_multi_peak",
              "nExposures","SN_B","fluxContamination_DESI","RV","RV_error"]
    frame.merge(controls,on="source_id")[["name","class",*fields]].to_csv(HERE/"control_matches.csv",index=False)
    summary = dict(catalogue_rows=len(phot),exact_dr3=len(entries),dr2_only=len(phot)-len(entries),
                   with_solution=int(valid.sum()), calibrated=int(par.notna().sum()),
                   precuts="All catalogue spectral classifications, exact Gaia DR3 designations only",
                   gates=[])
    for threshold in [3,5,10]:
        for apply_quality in [False,True]:
            keep = (frame.plx_snr>=threshold) & valid
            if apply_quality:
                keep &= quality
            s = frame.loc[keep]
            imax = s.vt_kms.idxmax()
            summary["gates"].append(dict(corrected_snr=threshold,quality=apply_quality,n=len(s),
                n_over_500=int((s.vt_kms>500).sum()),n_over_600=int((s.vt_kms>600).sum()),
                n_over_1000=int((s.vt_kms>1000).sum()),max_vt=round(s.vt_kms.max(),2),
                max_id=s.at[imax,"source_id"]))
    tail = frame[(frame.plx_snr>=3) & (frame.vt_kms>500)].sort_values("vt_kms",ascending=False)
    tail[fields].to_csv(HERE/"tail.csv",index=False)
    # Always retain the highest quality-eligible objects as an auditable envelope.
    frame[quality & (frame.plx_snr>=5)].nlargest(10,"vt_kms")[fields].to_csv(HERE/"top10.csv",index=False)
    (HERE/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path)
    screen(parser.parse_args().cache)
