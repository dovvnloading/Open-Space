"""Summarize frozen ground-census and force-law shards; verify the data and roster."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parent
GROUND_FIELDS="spkid name status reason reported_arc_days bulk_lines available_station_nights selected supported supported_arc_days used used_arc_days used_count_sites rms_coordinate_arcsec reduced_chi2 inverse_a_au_inv formal_sigma_inverse_a radius_au inbound_radial_au_d v_inf_km_s selected_sha256 bulk_astrometry_sha256 elapsed_seconds".split()
MODEL_FIELDS="spkid name model status reason selected supported arc_days preperihelion_count used used_arc_days rms_coordinate_arcsec reduced_chi2 inverse_a_au_inv formal_sigma_inverse_a radius_au inbound_radial_au_d v_inf_km_s parabolic_profile_delta_chi2 achieved_profile_inverse_a_at_50au a1 a2 a3 selected_sha256".split()
def write(path,cols,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=cols,extrasaction="ignore",lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

def checked(archive):
    data=json.loads(gzip.decompress(archive.read_bytes()))
    for r in data["objects"]:
        selected=r.get("selected_mpc80")
        if selected and hashlib.sha256(("\n".join(selected)+"\n").encode("ascii")).hexdigest()!=r["selected_sha256"]:
            raise ValueError(f"{archive}: astrometry hash mismatch for {r['spkid']}")
        yield r

def main():
    seen=set()
    targets=set()
    for name,count in [("census",578),("census_short",1615)]:
        evidence=ROOT/name/"evidence"
        parsed=[]
        for archive in sorted(evidence.glob("shard-*.json.gz")):
            for r in checked(archive):
                spk=r["spkid"]
                if spk in seen: raise ValueError("duplicate census ID "+spk)
                seen.add(spk)
                b=r.get("boundaries",{}).get("250au",{})
                if r["status"]=="inbound 250 au":
                    if abs(b["radius_au"]-250)>1e-5 or b["inbound_radial_au_d"]>=0: raise ValueError("invalid inbound crossing")
                    alpha,sigma=b["inverse_a_au_inv"],b["formal_sigma_inverse_a"]
                    if not sigma>0: raise ValueError("invalid formal error")
                    if alpha < -3*sigma or ((b.get("v_inf_km_s") or 0)>=1 and alpha< -1.5*sigma):
                        targets.add(spk)
                else: b={}
                parsed.append({**r,**b,"selected":len(r.get("selected_mpc80",[]))})
        if len(parsed)!=count: raise ValueError(f"expected {count} rows in {name}")
        write(evidence/"summary.csv",GROUND_FIELDS,sorted(parsed,key=lambda v:int(v["spkid"])))
    # The fragment has a separate pre-split parent calculation.
    targets.remove("1000158")
    root=ROOT/"robust"
    listed={r["spkid"] for f in root.glob("targets-*.json") for r in json.loads(f.read_text())}
    if targets!=listed: raise ValueError("force-law roster differs from census trigger")
    models=[]
    for f in sorted((root/"evidence").glob("*-part-*.json.gz")):
        for row in checked(f):
            for k,m in row["models"].items():
                b=m.get("boundaries",{}).get("250au",{})
                success=bool(b) and b["inbound_radial_au_d"]<0
                typ="inbound 250 au" if success else ("fit failed" if m.get("failed") else "no inbound 250 au")
                q=m.get("parabolic_profile",{})
                a=m.get("state_ssb_icrf",[])[6:]
                models.append({**row,**m,**(b if success else {}),"model":k,"status":typ,
                               "reason":m.get("failed",row.get("error")),"parabolic_profile_delta_chi2":q.get("delta_chi2_sky_only"),
                               "achieved_profile_inverse_a_at_50au":q.get("achieved"),**dict(zip(("a1","a2","a3"),a))})
    if {r["spkid"] for r in models}!=targets: raise ValueError("incomplete model outputs")
    write(root/"evidence/summary.csv",MODEL_FIELDS,sorted(models,key=lambda r:(int(r["spkid"]),r["model"])))

if __name__=="__main__": main()
