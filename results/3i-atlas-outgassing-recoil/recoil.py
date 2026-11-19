#!/usr/bin/env python3
"""Measured production-rate to RTN momentum budget; scenarios NOT inferred jet axes."""
from __future__ import annotations
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

HERE=Path(__file__).resolve().parent
OUT=HERE/"generated"
AU_M=149597870700.
DAY_S=86400.
AMU=1.66053906660e-27
MOLECULE_KG={"H2O":18.01528*AMU,"CO2":44.0095*AMU,"CO":28.0101*AMU}
RADII_KM=[.22,.42,.60,1.10,1.30,1.50,2.80]
DENSITIES=[300,500,1000]
GAS_SPEEDS=[200,500,1000]
EFFICIENCIES=[.1,.3,.5,1.]

def dump(path,obj):
    (OUT/path).write_text(json.dumps(obj,indent=2,allow_nan=False)+"\n")

def get_jpl(sbdb):
    parameters={}
    for q in sbdb.get("orbit",{}).get("model_pars",[]):
        name=str(q.get("name","")).upper()
        if name in {"A1","A2","A3","DT"}:
            try:
                val=float(q["value"]); sigma=float(q["sigma"]) if q.get("sigma") else None
            except (TypeError,ValueError): continue
            parameters[name]={"value":val,"sigma":sigma,"units":q.get("units"),
                              "kind":q.get("kind")}
    return parameters

def main():
    orbit=json.loads((OUT/"orbit.json").read_text())
    sbdb=json.loads((OUT/"sbdb.json").read_text())
    fit=orbit["fit_models"]["r2_rtn"]
    coeff=fit["radial_transverse_normal_A_1e8_au_d2"]
    sigma=fit["A_formal_1sigma_1e8_au_d2"]
    if not all(math.isfinite(x) for x in coeff):
        raise RuntimeError("Unusable independent orbit coefficients")
    norm=math.sqrt(sum(x*x for x in coeff))
    pars=get_jpl(sbdb)
    jpl=[pars.get(k,{}).get("value") for k in ["A1","A2","A3"]]
    if all(x is not None for x in jpl):
        # JPL model parameters are in au/d^2, not 10^-8 au/d^2.
        jpl=[x*1e8 for x in jpl]
    else: jpl=None
    comparisons={"independent_inverse_square_A_1e8_au_d2":coeff,
        "independent_formal_sigma_1e8_au_d2":sigma,
        "independent_norm_1e8_au_d2":norm,
        "independent_direction_RTN_unit":[x/norm for x in coeff],
        "jpl_model_parameters":pars,
        "jpl_A_1e8_au_d2":jpl,
        "jpl_minus_independent_A_1e8_au_d2":[jpl[i]-coeff[i] for i in range(3)] if jpl else None,
        "warning":"JPL delayed peak DT (when present) and different weighting/systematics are not fitted by our symmetric r^-2 model. Formal sigmas are not a systematic error budget."}
    dump("jpl_comparison.json",comparisons)
    with (HERE/"production.csv").open(newline="") as handle:
        rows=list(csv.DictReader(handle))
    groups=defaultdict(list)
    for row in rows: groups[row["epoch"]].append(row)
    epoch_results=[]
    grid=[]
    for epoch,members in sorted(groups.items(),key=lambda x:x[1][0]["date"]):
        r=float(members[0]["r_au"])
        if any(abs(float(row["r_au"])-r)>.03 for row in members):
            raise ValueError("Group has inconsistent heliocentric distances: "+epoch)
        # Gas mass rate from directly measured H2O, CO2, CO and a separate
        # explicit OH -> water PROXY. Upper limits are included as though
        # their gas is present: a favorable scenario, not a total-gas bound.
        kg_s=0.
        gas_members=[]
        for row in members:
            species=row["species"]
            q=float(row["q_s"])
            if species in MOLECULE_KG:
                value=q*MOLECULE_KG[species]
                kg_s+=value
                gas_members.append({"species":species,"q_molecules_s":q,
                    "kg_s":value,"kind":row["kind"],"reference":row["reference"],
                    "note":row["note"]})
            elif species=="OH":
                qwater=q*1.1
                value=qwater*MOLECULE_KG["H2O"]
                kg_s+=value
                gas_members.append({"species":"H2O_from_OH_proxy","q_OH_s":q,
                    "assumed_H2O_to_OH_factor":1.1,"kg_s":value,
                    "reference":row["reference"],"kind":row["kind"]})
        a_RTN=[x*1e-8*AU_M/(DAY_S*DAY_S)/r**2 for x in coeff]
        a=math.sqrt(sum(x*x for x in a_RTN))
        entry={"epoch":epoch,"date":members[0]["date"],"r_au":r,
               "measured_species":[{"species":m["species"],"kind":m["kind"],
                   "Q_s":float(m["q_s"]),"Q_sigma_s":float(m["sigma_s"]) if m["sigma_s"] else None,
                   "facility":m["facility"],"reference":m["reference"]} for m in members],
               "included_gas_mass_loss_kg_s":kg_s,"gases_in_recoil":gas_members,
               "fitted_RTN_accel_m_s2":a_RTN,"fitted_norm_m_s2":a,
               "missing_species_are_not_zero":True,
               "jet_geometry_unmeasured":True}
        epoch_results.append(entry)
        if kg_s>0:
            for radius in RADII_KM:
                for density in DENSITIES:
                    mass=4/3*math.pi*(radius*1000)**3*density
                    force_need=mass*a
                    for vel in GAS_SPEEDS:
                        force_available=kg_s*vel
                        eta_required=force_need/force_available
                        grid.append({"epoch":epoch,"r_au":r,"radius_km":radius,"density_kg_m3":density,
                            "gas_speed_m_s":vel,"nucleus_mass_kg":mass,
                            "measured_species_mass_flow_kg_s":kg_s,
                            "required_force_N":force_need,
                            "available_force_N_at_eta_1":force_available,
                            "minimum_efficiency_needed":eta_required,
                            "magnitude_feasible_at_eta_le_1":eta_required<=1})
    dump("recoil_epochs.json",epoch_results)
    dump("recoil_grid.json",grid)
    # Pure-input physical sanity invariants: kg/s * m/s / kg -> m/s^2.
    assert all(x["minimum_efficiency_needed"]>=0 for x in grid)
    assert all(abs(x["required_force_N"]/x["available_force_N_at_eta_1"]-
                 x["minimum_efficiency_needed"])<1e-8 for x in grid)
    fid=lambda e,R,rho=500,v=500: next((q for q in grid if q["epoch"]==e and q["radius_km"]==R and q["density_kg_m3"]==rho and q["gas_speed_m_s"]==v),None)
    summary={"fit":comparisons,"n_production_entries":len(rows),"n_epoch_groups":len(groups),
        "n_grid_rows":len(grid),"physical_limits":{
            "radii_km":RADII_KM,"rho_kg_m3":DENSITIES,
            "exhaust_speed_m_s":GAS_SPEEDS,"momentum_efficiency_sample":EFFICIENCIES,
            "mass_basis":"spherical uniform density",
            "required_RTN_direction":comparisons["independent_direction_RTN_unit"],
            "sum_species_limit":"Only listed contemporaneous gas fluxes; unmeasured species and distributed ice sources are not constrained.",
            "direction_limit":"Q data provide no jet direction; eta applies to the vector sum and can be near zero even if gas is abundant.",
            "model_limit":"Symmetric r^-2 RTN independent fit; JPL has delayed response; outburst timing is sparse."},
        "example_eta_at_500m_s_density500":{k:fid(k,R)["minimum_efficiency_needed"] if fid(k,R) else None
                 for k,R in [("jwst_aug06",.42),("jwst_aug06",1.3),
                    ("juice_nov02",.42),("juice_nov02",1.3),
                    ("jwst_dec23",.42),("jwst_dec23",1.3)]}}
    dump("summary.json",summary)
    report=[
        "# 3I/ATLAS independent fit and outgassing recoil: computed data",
        "",
        "**Claim scope:** this is a physical sufficiency grid, **not** an identification of a jet direction or a proof of anomalous forces. All numbers are reproduced from the frozen input rules and public sources; the independent fit omits JPL's peak-delay DT and astrometric debiasing.",
        "",
        f"Retained full-arc sky positions: **{fit['fit_count']}** from {orbit['n_supported']} supported; independent weighted RMS **{fit['rms_fit_per_coordinate_arcsec']:.3f} arcsec per coordinate**.",
        "",
        "The independently fitted symmetric inverse-square RTN A coefficients, in 10^-8 au/day^2: "+", ".join(f"{v:+.4f} +/- {e:.4f}" for v,e in zip(coeff,sigma))+".",
        "",
        f"Fitted unit RTN direction: {comparisons['independent_direction_RTN_unit']}. The nongravitational acceleration at heliocentric distance r is fitted as A/r^2.",
        "",
        "### Flux-to-recoil scenario (density 500 kg/m^3; 500 m/s; fully collimated efficiency eta <= 1)",
        "",
        "| Gas measurement epoch | Radius (km) | Required eta | Magnitude possible at eta <= 1? |",
        "| --- | ---: | ---: | --- |",
    ]
    for k in ["jwst_aug06","juice_nov02","jwst_dec23"]:
        for R in [.42,1.3]:
            g=fid(k,R)
            if g:
                report.append(f"| {k} | {R:.2f} | {g['minimum_efficiency_needed']:.3g} | {'Yes' if g['magnitude_feasible_at_eta_le_1'] else 'No'} |")
    report += ["","The numerical grid also covers 0.22–2.8 km radii, 300–1000 kg/m^3 densities and 200–1000 m/s exhaust speeds. A 'No' in one row is **not** proof of excess acceleration; missing gas species, uncertain active fraction, a smaller nucleus or alternative gas speeds remain possible.","",
      "Scalar measured Q cannot establish the radial, transverse and normal momentum vector: observed gas jets and nucleus spin/orientation would be needed. Large-aperture OH can arise from distributed icy grains. A fit with fixed r^-2 cannot resolve perihelion asymmetry, and fitting uncertainty is not exclusively statistical.",
      "","Detailed values: [orbit.json](orbit.json), [jpl_comparison.json](jpl_comparison.json), [recoil_epochs.json](recoil_epochs.json), [recoil_grid.json](recoil_grid.json) and [retrieval receipts](retrieval_receipts.json).",
      ""]
    (OUT/"REPORT.md").write_text("\n".join(report))

if __name__=="__main__": main()
