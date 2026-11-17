"""Replay a composition-first screen on the hash-verified DESI/Gaia parent."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from zero_point import zpt

HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser()
parser.add_argument("cache", type=Path)
cache = parser.parse_args().cache
receipt = json.loads((HERE.parent/"desi"/"receipt.json").read_text())

def verify(p, expected):
    if hashlib.sha256(p.read_bytes()).hexdigest() != expected:
        raise ValueError(f"Fingerprint differs for {p.name}")

catfile = cache/"DESI_DR1_catalogue_1.0.csv.gz"
verify(catfile, receipt["catalogue_sha256"])
cat = pd.read_csv(catfile, dtype={"designation":str},low_memory=False)
isdr3 = cat.designation.str.startswith("Gaia DR3 ",na=False)
pattern = r"^D[A-Z]*[QZO][A-Z]*$"
unmatched = cat.loc[~isdr3 & cat.specType.str.match(pattern,na=False),["WDJname","designation","specType"]]
unmatched.to_csv(HERE/"dr2_only.csv",index=False)
cat = cat.loc[isdr3].copy()
cat["source_id"] = cat.designation.str.split().str[-1]
chunks = []
for batch in receipt["batches"]:
    f = cache/f"gaia_{batch['batch']:03}.csv"
    verify(f,batch["response_sha256"])
    chunks.append(pd.read_csv(f, dtype={"source_id":str}))
data = pd.concat(chunks).merge(cat, on="source_id", validate="one_to_one")
assert len(data) == len(cat) == receipt["dr3_ids"]
# Retain *any* Q, Z or O suffix on the Swan WD type, including mixed DAZ and DAQ.
selected = data.specType.str.match(pattern,na=False)
d = data.loc[selected].copy()
valid = d.astrometric_params_solved.isin([31,95]) & d.parallax.notna() & (d.parallax_error>0)
pm = np.hypot(d.pmra,d.pmdec)
quality = (d.ruwe.gt(0) & d.ruwe.lt(1.4) & d.visibility_periods_used.ge(9)
           & d.duplicated_source.eq(False) & d.ipd_frac_multi_peak.le(2)
           & d.ipd_gof_harmonic_amplitude.lt(0.1))
zpt.load_tables()
d["zp_mas"] = np.nan
cols = ["phot_g_mean_mag","nu_eff_used_in_astrometry","pseudocolour","ecl_lat","astrometric_params_solved"]
d.loc[valid,"zp_mas"] = zpt.get_zpt(*(d.loc[valid,col].to_numpy() for col in cols), _warnings=False)
d["pm_masyr"] = pm
common = valid & (pm>35) & quality
summary = {"dr2_only_class_objects":len(unmatched),"spectral_catalogue":receipt["rows"],"exact_dr3":len(data),"class_counts":d.specType.value_counts().sort_index().to_dict(),
           "class_gate":len(d),"astrometric_solution":int(valid.sum()),"pm35":int((valid & (pm>35)).sum()),
           "pm35_quality":int(common.sum()),"correction_available_after_quality":int((common & d.zp_mas.notna()).sum())}
for key, plx in [("raw",d.parallax),("corrected",d.parallax-d.zp_mas)]:
    den = plx+2*d.parallax_error
    d[key+"_stat_kms"] = 4.74047*pm/den.where(den>0)
    eligible = common & den.gt(0)
    d[key+"_pass"] = eligible & d[key+"_stat_kms"].gt(350)
    summary[key] = {"finite_parent":int(eligible.sum()),"positive_denom_failures":int((common & den.notna() & den.le(0)).sum()),
                    "over350":int(d[key+"_pass"].sum())}
zmissing = d[common & d.zp_mas.isna()]
summary["uncalibrated_max_raw_stat_kms"] = float(zmissing.raw_stat_kms.max())
rows = d[d.raw_pass | d.corrected_pass].copy()
cols = ["source_id","WDJname","specType","specType_confidence","pm_masyr","parallax","parallax_error","zp_mas","ruwe",
        "raw_stat_kms","corrected_stat_kms","raw_pass","corrected_pass","nExposures","SN_B","RV","RV_error"]
rows.sort_values("corrected_stat_kms",ascending=False)[cols].to_csv(HERE/"targets.csv",index=False)
ctrl = pd.read_csv(HERE.parent/"control"/"controls.csv",dtype={"source_id":str})
matches=ctrl.merge(data[["source_id","specType"]],on="source_id",how="left").merge(
    d[["source_id","raw_pass","corrected_pass"]],on="source_id",how="left")
summary["control_recovery"] = {group:{"reference":len(g),"in_parent":int(g.specType.notna().sum()),
   "in_class":int(g.raw_pass.notna().sum()),"raw_pass":int(g.raw_pass.eq(True).sum()),
   "corrected_pass":int(g.corrected_pass.eq(True).sum())} for group,g in matches.groupby("class")}
(HERE/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
print(json.dumps(summary,indent=2))
