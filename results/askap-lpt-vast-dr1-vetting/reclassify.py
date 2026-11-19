#!/usr/bin/env python3
"""Independent science-class re-evaluation of already downloaded published evidence.
Never promote mere positional coincidence to a flare-star or confirmed AGN/LPT.
"""
import csv,json,math
from collections import Counter,defaultdict
from pathlib import Path
D=Path("results/askap-lpt-vast-dr1-vetting/output")
def number(v):
 try:
  x=float(v)
  return x if math.isfinite(x) else None
 except (TypeError,ValueError):return None
def read(name):
 with (D/name).open(newline="") as f:return list(csv.DictReader(f))
def write(name,rows):
 keys=list(rows[0])
 with (D/name).open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=keys,extrasaction="ignore");w.writeheader();w.writerows(rows)
def rebuild():
 old=read("classifications.csv")
 c=json.loads((D/"counterpart_receipts.json").read_text())
 epo=read("epoch_image_audit.csv")
 assert len(old)==60 and len({x["source"] for x in old})==60 and len(epo)==693
 flagged={x["source"] for x in epo if x.get("flags")}
 counts=Counter();cat_counts=defaultdict(Counter);gaia_stars=[];ambig=[]
 for row in old:
  sid=row["source"];match=c[sid]
  gaia=match.get("gaia",{})
  sep=number(gaia.get("arcsec"))
  rec=gaia.get("record",{})
  plx=number(rec.get("Plx"));err=number(rec.get("e_Plx"));mag=number(rec.get("Gmag"))
  significance=plx/err if plx is not None and err is not None and err>0 else None
  close=gaia.get("status")=="match" and sep is not None and sep<2
  strong_star=bool(close and ((significance is not None and significance>=5) or (mag is not None and mag<8)))
  row["gaia_parallax_significance"]=round(significance,2) if significance is not None else ""
  row["gaia_arcsec"]=sep if sep is not None else ""
  ps=match.get("panstarrs",{})
  prec=ps.get("record",{})
  nd=number(prec.get("nDetections"))
  gmag=number(prec.get("gMeanPSFMag"))
  rmag=number(prec.get("rMeanPSFMag"))
  reliable_ps=bool(ps.get("status")=="match" and (nd or 0)>=2 and ((gmag is not None and 0<gmag<35) or (rmag is not None and 0<rmag<35)))
  row["panstarrs_photometry_valid"]=reliable_ps
  for name in ("legacy_dr10","unwise","panstarrs","swift"):
   row[name]=match.get(name,{}).get("status","unavailable")
  for name,v in match.items():cat_counts[name][v.get("status","unavailable")]+=1
  if row["psr"]:
   row["class"]="known_pulsar"
   row["class_basis"]="ATNF prior match within 2 arcsec (from frozen PR 200)"
  elif strong_star:
   row["class"]="plausible_stellar_counterpart_no_confirmed_flare"
   row["class_basis"]="Gaia within 2 arcsec and either parallax >=5 sigma or G<8; no flare proven"
   gaia_stars.append(sid)
  elif close:
   row["class"]="gaia_position_match_source_nature_undetermined"
   row["class_basis"]="Gaia position match with no secure stellar parallax; not proven star or AGN"
   ambig.append(sid)
  else:
   dr=match.get("legacy_dr10",{})
   kind=str(dr.get("record",{}).get("type","")).upper()
   if dr.get("status")=="match" and kind in ("EXP","DEV","SER","REX"):
    row["class"]="possible_extragalactic_host_agn_or_scintillation"
    row["class_basis"]="positional association with extended Legacy object only; AGN nature unproven"
   else:
    row["class"]="unresolved_polarisation_and_host_incomplete"
    row["class_basis"]="no confirmed physical identification; missing V and some catalogue coverage"
  counts[row["class"]]+=1
 assert sum(counts.values())==60
 assert len(gaia_stars)==2 and len(ambig)==6,(gaia_stars,ambig)
 write("classifications.csv",old)
 s=json.loads((D/"summary.json").read_text())
 s.update({"classes":dict(counts),"onoff_classes":{x["source"]:x["class"] for x in old if x["onoff"]=="True"},
   "racs_i_noncontemporaneous_spectra":sum(x.get("racs_alpha","")!="" for x in old),
   "native_automated_flagged_epochs":sum(bool(x.get("flags")) for x in epo),
   "native_automated_flagged_sources":len(flagged),
   "catalogue_status_counts":{k:dict(v) for k,v in cat_counts.items()},
   "gaia_stellar_candidates_with_support":gaia_stars,
   "gaia_ambiguous_position_matches":ambig,
   "scope":"native Stokes I cutout inspection and public catalogue matches only; no calibrated Stokes V/short-cadence time series",
   "image_claim_limit":"119px centered cutouts cannot reject full ASKAP field-edge and distant sidelobe contamination",
   "stokes_v_non_detection_claimed":False,
   "discovery_validated":False})
 (D/"summary.json").write_text(json.dumps(s,indent=2,sort_keys=True)+"\n")
 ledger={"method":"independent referee reclassification of committed source/counterpart receipts",
   "inputs":["classifications.csv pre-revision","counterpart_receipts.json","epoch_image_audit.csv"],
   "Gaia_rule":"within 2 arcsec AND parallax significance >=5 OR G<8; only plausible star, not flare",
   "gaia_plausible_stars":gaia_stars,"gaia_ambiguous_matches":ambig,
   "catalogue_nonmatches":"not upper limits; require footprint, exposure, optical depth and reliability",
   "panstarrs_requirement":"nDetections>=2 and valid 0<g/r PSF magnitude<35; many zero-detection catalogue rows cannot prove optical counterparts",
   "classification_counts":dict(counts),
   "no_new_LPT":True}
 (D/"reclassification_receipts.json").write_text(json.dumps(ledger,indent=2,sort_keys=True)+"\n")
 print(json.dumps(ledger,indent=2))
if __name__=="__main__":rebuild()
