"""Calibrated *screening* coincidences for frozen Cao 2024 source units.

The null shifts all sources by the same longitude offsets at fixed Galactic latitude,
excluding |delta-l|<5 degrees. This is a conditional latitude null, NOT a physical
identification probability or a control for spatially varying survey exposure.
"""
from __future__ import annotations
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.coordinates import SkyCoord
import astropy.units as u

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
CAT=OUT/"catalogs"
N_SHIFTS=512
SEED=20261119
CLASSES={
    "GeV_4FGL_DR4":("fermi_dr4.csv","all_sky"),
    "pulsar_ATNF":("atnf_pulsars.csv","all_sky"),
    "SNR_Green2024":("green_snrs.csv","all_sky"),
    "SMGPS_new_SNR_candidate":("smgps_candidates.csv","smgps"),
    "SMGPS_known_SNRcat_subset":("smgps_snrcat.csv","smgps"),
    "eRASS1_hard_Xray":("erass1_hard.csv","erass_approx"),
    "Gaia_DR3_young_OC":("gaia_young_clusters.csv","all_sky"),
    "TeV_HESS_HGPS_legacy":("hess_hgps.csv","hess_approx"),
    "TeV_HAWC4_2026":("hawc4_sources.csv","hawc4_approx"),
}
PC_KPC=3.0856775814913673e21
TEV_ERG=1.602176634

def js(path,obj):
    Path(path).write_text(json.dumps(obj,indent=2,sort_keys=True,default=str)+"\n")

def nfloat(v):
    try:
        z=float(v)
        return z if np.isfinite(z) else float("nan")
    except Exception:return float("nan")

def vec(l,b):
    l=np.deg2rad(np.asarray(l,dtype=float));b=np.deg2rad(np.asarray(b,dtype=float))
    return np.column_stack((np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)))

def angular(chord):
    return np.rad2deg(2*np.arcsin(np.minimum(1,np.asarray(chord)/2)))

def footprint(name,l,b):
    l=np.asarray(l,dtype=float)%360; b=np.asarray(b,dtype=float)
    if name=="all_sky":return np.ones(l.shape,dtype=bool)
    if name=="erass_approx":return (l>=180)&(l<360)
    if name=="smgps":
        return (((l>=252)&(l<=358))|((l>=2)&(l<=61))) & (np.abs(b)<=1.5)
    if name=="hess_approx":
        return ((l>=250)|(l<=65))&(np.abs(b)<=3.0)
    if name=="hawc4_approx":
        dec=SkyCoord(l=l*u.deg,b=b*u.deg,frame="galactic").icrs.dec.deg
        return (dec>=-25.0)&(dec<=65.0)
    raise ValueError(name)

def search_tree(tree,cat,ls,bs):
    """Nearest angular separation minus counterpart radius, and nearest row ID.

    The per-match and null queries use *identical* radius treatment.
    64 neighbours is safe for the screened sparse source classes; a dense
    catalogue would need a dedicated radius query and separate validation.
    """
    k=min(64,len(cat))
    if k<1:raise ValueError("empty catalogue")
    ds,ids=tree.query(vec(ls,bs),k=k)
    if k==1:
        ds=ds[:,None];ids=ids[:,None]
    sep=angular(ds)
    radius=cat["radius_deg"].fillna(0).to_numpy(float)
    mar=sep-radius[ids]
    best=np.argmin(mar,axis=1)
    rr=np.arange(len(best))
    return mar[rr,best],ids[rr,best],sep[rr,best]

def flux_25_100(n0,index):
    """KM2A fitted power-law energy flux, 25-100 TeV only.

    Cao 2024 KM2A N0 unit: 1e-16 TeV^-1 cm^-2 s^-1 at E0=50 TeV.
    Power-law extrapolation is NOT an observed calorimetric luminosity.
    """
    n0=nfloat(n0);g=nfloat(index)
    if not (0<n0<100 and 1<g<6):return float("nan")
    a,b=25.0,100.0
    integ=np.log(b/a) if abs(g-2)<1e-10 else (b**(2-g)-a**(2-g))/(2-g)
    return n0*1e-16*(50**g)*integ*TEV_ERG

def energetics(source,counterpart,sep):
    d=nfloat(counterpart.get("distance_kpc"))
    edot=nfloat(counterpart.get("edot_erg_s"))
    age=nfloat(counterpart.get("age_yr"))
    flux=flux_25_100(source.get("km2a_n0"),source.get("km2a_index"))
    radius=nfloat(source.get("r39_deg"))
    out={"distance_kpc":d,"edot_erg_s":edot,"age_yr":age,
         "flux_25_100_erg_cm2_s":flux,
         "projected_offset_pc":np.deg2rad(sep)*d*1000 if np.isfinite(d) else np.nan,
         "r39_pc":np.deg2rad(radius)*d*1000 if np.isfinite(radius) and np.isfinite(d) else np.nan,
         "luminosity_25_100_erg_s":4*math.pi*(d*PC_KPC)**2*flux if np.isfinite(d) and d>0 and np.isfinite(flux) else np.nan}
    out["gamma_to_edot_25_100"]=out["luminosity_25_100_erg_s"]/edot if np.isfinite(edot) and edot>0 and np.isfinite(out["luminosity_25_100_erg_s"]) else np.nan
    out["energy_screen"] = ("unknown_distance_or_edot" if not(np.isfinite(d) and 0<d<30 and np.isfinite(edot) and edot>0 and np.isfinite(out["gamma_to_edot_25_100"]))
                            else "implausible_eta_gt_one" if out["gamma_to_edot_25_100"]>1
                            else "energetically_possible_not_proven")
    return out

def main():
    p=CAT/"1lhaaso_sources.csv"
    if not p.exists():raise SystemExit("Missing validated source table; cannot perform a match.")
    sources=pd.read_csv(p)
    assert len(sources)==90 and int(sources.uhe.sum())==43
    rec=json.loads((OUT/"receipts.json").read_text())
    primary=sources[~sources.known_tev].copy().reset_index(drop=True)
    # All unassociated 1LHAASO sources are screened, not just 43 UHE.
    shift=np.random.default_rng(SEED).uniform(5.0,355.0,size=N_SHIFTS)
    rows=[];pairs=[];global_rows=[]
    allsky_vectors={}; candidate_counts={}
    for classname,(filename,mask) in CLASSES.items():
        file=CAT/filename
        receipt_name={
            "GeV_4FGL_DR4":"fermi_4fgl_dr4_v35",
            "pulsar_ATNF":"atnf_psrcat_v2_latest",
            "SNR_Green2024":"green_2024",
            "SMGPS_new_SNR_candidate":"smgps_new_snr_candidates",
            "SMGPS_known_SNRcat_subset":"snrcat_smgps_overlap",
            "eRASS1_hard_Xray":"erass1_dr1_hard",
            "Gaia_DR3_young_OC":"gaia_dr3_hunt_ref",
            "TeV_HESS_HGPS_legacy":"hess_hgps_legacy",
            "TeV_HAWC4_2026":"hawc4_arxiv_2602_00263"}[classname]
        if not file.exists() or rec.get(receipt_name,{}).get("status")!="loaded":
            for _,s in primary.iterrows():
                rows.append({"source":s["name"],"uhe":bool(s.uhe),"class":classname,
                             "status":"catalogue_pending","covered":False})
            global_rows.append({"class":classname,"status":"catalogue_pending","observed_matches":None})
            continue
        catalog=pd.read_csv(file)
        catalog=catalog[np.isfinite(pd.to_numeric(catalog.ra,errors="coerce"))
                        & np.isfinite(pd.to_numeric(catalog.dec,errors="coerce"))].reset_index(drop=True)
        if len(catalog)==0:
            global_rows.append({"class":classname,"status":"empty_catalogue_unusable"})
            continue
        tree=cKDTree(vec(catalog.l,catalog.b))
        matched_by_shift=np.zeros(N_SHIFTS,dtype=int)
        observed_total=0
        included=0
        for _,s in primary.iterrows():
            name=str(s["name"])
            l=float(s["l"]);b=float(s["b"])
            # A fixed screening aperture: source 95%-stat uncertainty plus
            # r39 39%-containment radius and 0.2 deg systematic allowance.
            radius=np.clip(nfloat(s["r95_deg"])+nfloat(s["r39_deg"])+0.20,0.35,1.50)
            covered=bool(footprint(mask,[l],[b])[0])
            control_l=(l+shift)%360
            use=footprint(mask,control_l,np.full(N_SHIFTS,b))
            # Survey footprint is not a real exposure map; status flagged.
            status="screened" if mask=="all_sky" else "screened_approx_footprint"
            if not covered:
                rows.append({"source":name,"uhe":bool(s.uhe),"class":classname,
                             "status":"outside_survey_footprint","covered":False,
                             "control_n":0})
                continue
            if int(use.sum())<64:
                rows.append({"source":name,"uhe":bool(s.uhe),"class":classname,
                             "status":"too_few_covered_controls","covered":True,
                             "control_n":int(use.sum())})
                continue
            included+=1
            obs,ids,seps=search_tree(tree,catalog,[l],[b])
            margin=float(obs[0])
            matches=margin<=radius
            ctrl,_,_=search_tree(tree,catalog,control_l[use],np.full(use.sum(),b))
            ctrl=np.asarray(ctrl,dtype=float)
            nearest_p=(1+int(np.sum(ctrl<=margin)))/(1+len(ctrl))
            chance=(1+int(np.sum(ctrl<=radius)))/(1+len(ctrl))
            if matches:observed_total+=1
            if mask=="all_sky":
                matched_by_shift+=(ctrl<=radius).astype(int)
            rows.append({"source":name,"uhe":bool(s.uhe),"class":classname,
                         "status":status,"covered":True,"screen_radius_deg":round(float(radius),5),
                         "nearest_name":str(catalog.iloc[int(ids[0])]["name"]),
                         "nearest_deg":round(float(seps[0]),6),
                         "nearest_margin_deg":round(margin,6),
                         "n_matched":int(matches),
                         "chance_any_in_radius":chance,
                         "p_nearest_margin":nearest_p,
                         "control_hits":int(np.sum(ctrl<=radius)),
                         "control_n":len(ctrl),
                         "null_lon_shifts_seed":SEED})
            if matches:
                # Save every candidate within the aperture, not just nearest.
                neighbour_ids=tree.query_ball_point(vec([l],[b])[0],2*np.sin(np.deg2rad(radius+max(catalog.radius_deg.max(),0))/2))
                for index in neighbour_ids:
                    other=catalog.iloc[int(index)]
                    cossep=np.clip(np.dot(vec([l],[b])[0],vec([other["l"]],[other["b"]])[0]),-1,1)
                    sep=np.rad2deg(np.arccos(cossep))
                    shape=nfloat(other.get("radius_deg",0))
                    if sep>radius+shape:continue
                    pr=(1+int(np.sum(ctrl<=sep-shape)))/(1+len(ctrl))
                    item={"source":name,"class":classname,
                          "counterpart":str(other["name"]),"separation_deg":round(sep,6),
                          "counterpart_radius_deg":round(shape,5),
                          "control_p_at_pair_margin":pr,
                          "source_screen_radius_deg":radius,
                          "association_status":"geometrical_only"}
                    if classname=="pulsar_ATNF":item.update(energetics(s,other,sep))
                    if classname=="Gaia_DR3_young_OC":
                        item["cluster_log_age_yr"]=nfloat(other.get("log_age"))
                        item["cluster_distance_pc"]=nfloat(other.get("distance_pc"))
                    pairs.append(item)
        if mask=="all_sky" and included>0:
            p_global=(1+int(np.sum(matched_by_shift>=observed_total)))/(1+N_SHIFTS)
            null_mean=float(np.mean(matched_by_shift))
            null_sd=float(np.std(matched_by_shift,ddof=1))
        else:
            p_global=np.nan;null_mean=np.nan;null_sd=np.nan
        global_rows.append({"class":classname,"status":"screened" if mask=="all_sky" else "approximate_footprint_only",
                            "eligible_sources":included,"observed_matches":observed_total,
                            "null_mean_matches":null_mean,
                            "null_sd_matches":null_sd,
                            "p_global_ge_observed":p_global,
                            "n_common_shifts":N_SHIFTS if mask=="all_sky" else 0,
                            "footprint":mask})
    pairs=pd.DataFrame(pairs)
    screen=pd.DataFrame(rows)
    global_report=pd.DataFrame(global_rows)
    # The primary published TeV-unassociated set is NOT automatically the
    # strict unassociated set. Missing Fermi/ATNF/Green => unknown.
    core=["GeV_4FGL_DR4","pulsar_ATNF","SNR_Green2024"]
    strict=[]
    for _,s in primary.iterrows():
        sel=screen[(screen.source==s["name"]) & screen["class"].isin(core)]
        established= len(sel)==3 and (sel.status=="screened").all()
        strict.append({"source":s["name"],"uhe":bool(s.uhe),
                       "published_tev_dark":True,
                       "status":"screen_complete" if established else "incomplete_class_coverage",
                       "no_geometrical_GeV_PSR_SNR":bool(established and (sel.n_matched.fillna(0)==0).all()) if established else None})
    pd.DataFrame(strict).to_csv(OUT/"strict_dark.csv",index=False)
    # Bonferroni correction across all pairs tested, without pretending
    # independent positions or turning a small p into a physical ID.
    if len(pairs):
        pvals=pd.to_numeric(pairs.control_p_at_pair_margin,errors="coerce")
        pairs["bonferroni_all_pairs"]=np.minimum(1,pvals*len(pairs))
        pairs=pairs.sort_values(["bonferroni_all_pairs","control_p_at_pair_margin","separation_deg"])
    pairs.to_csv(OUT/"candidate_pairs.csv",index=False,float_format="%.9g")
    screen.to_csv(OUT/"source_class_controls.csv",index=False,float_format="%.9g")
    global_report.to_csv(OUT/"class_census.csv",index=False,float_format="%.9g")
    summary={
      "source_count":len(sources),"uhe_source_count":int(sources.uhe.sum()),
      "published_tev_unassociated":len(primary),
      "published_tev_unassociated_uhe":int(primary.uhe.sum()),
      "strict_no_GeV_PSR_SNR_if_all_loaded":int(sum(x["no_geometrical_GeV_PSR_SNR"] is True for x in strict)),
      "core_three_GeV_ATNF_Green_screen_complete":bool(all(x["status"]=="screen_complete" for x in strict)),
      "strict_census_complete":False,
      "strict_census_scope_note":"The 5 no-match cases are only for 4FGL-DR4, ATNF and Green; not global non-detections across missing radio/X-ray/cluster catalogues.",
      "all_requested_catalogues_complete":False,
      "controls_per_class_expected":N_SHIFTS,
      "shift_seed":SEED,
      "chance_p_definition":"(1+number of latitude-matched longitude controls with a match)/(1+eligible control count)",
      "minimum_resolvable_p":1/(N_SHIFTS+1),
      "radius_rule_deg":"clip(1LHAASO ePos95stat + r39 + 0.2, 0.35, 1.5); counterpart catalog radius added",
      "physical_identifications":0,
      "discovery_claim":False,
      "caveats":["longitude shifts do not reproduce longitude-dependent sensitivity/crowding",
                 "surveys with approximate sky footprints cannot certify absence or global significance",
                 "radio SNR candidates, X-ray points and GeV coincidences are not proof of same accelerator",
                 "KM2A power-law integration 25-100 TeV requires independent spectral consistency checks",
                 "pair p values are positional screening only; Bonferroni shown for multiplicity",
                 "known literature association proposals must be checked before any novelty label"]}
    js(OUT/"summary.json",summary)
    print("PUBLISHED DARK:",len(primary),"UHE:",int(primary.uhe.sum()))
    print("STRICT FULL COVERAGE:",summary["strict_census_complete"])
    print(global_report.to_string(index=False))

if __name__=="__main__":main()
