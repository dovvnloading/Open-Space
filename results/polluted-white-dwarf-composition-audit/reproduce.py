"""Re-run both abundance checks. Usage: python reproduce.py [--cache directory]."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
from itertools import combinations
from urllib.request import urlopen

import numpy as np

HERE = Path(__file__).resolve().parent
ELEMENTS = ("Na", "Mg", "Ca", "Cr", "Fe")
NUMBERS = (11, 12, 20, 24, 26)
ROWS = ("WD J2236-5548", "O'Brien 2023")
CHO = {"CC", "CH", "CL", "CLL", "CE", "R", "CHONDRITES"}


def inputs(cache):
    cache.mkdir(parents=True, exist_ok=True)
    for name, info in json.loads((HERE / "sources.json").read_text()).items():
        path = cache / name
        if not path.exists():
            with urlopen(info["url"], timeout=60) as response:
                path.write_bytes(response.read())
        if hashlib.sha256(path.read_bytes()).hexdigest() != info["sha256"]:
            raise ValueError(f"Source checksum mismatch for {name}")


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def measurement(row, element):
    key = f"log({element}/H(e))"
    return float(row[key]), float(row[key+"e"])


def cool_grid(cache, overshoot, temperature, gravity, calcium, numbers=NUMBERS):
    def interpolate_file(gname):
        by_temp = {}
        with (cache / f"timescales_He_g{gname}_ov{overshoot}.csv").open() as f:
            for key, *values in csv.reader(f):
                if key == "Teff":
                    temperature_key = float(values[0])
                    by_temp[temperature_key] = {}
                elif key == "CaHe":
                    abundances = np.array(values[::-1], dtype=float)
                elif key in {str(z) for z in numbers}:
                    by_temp[temperature_key][int(key)] = np.interp(
                        calcium, abundances, np.array(values[::-1], dtype=float))
        temperatures = sorted(by_temp)
        return np.array([np.interp(temperature, temperatures,
                        [by_temp[t][z] for t in temperatures]) for z in numbers])
    lower, upper = [interpolate_file(gname) for gname in (800, 825)]
    return lower + (upper-lower) * (gravity-8.0)/0.25


def warm_grid(cache, temperature, gravity, filenames=("SV_He_80.csv", "SV_He_85.csv"),
              numbers=(20, 22, 26, 28)):
    def interpolate_file(gname):
        with (cache / gname).open() as f:
            rows = dict((row[0], np.array(row[1:], dtype=float)) for row in csv.reader(f))
        return np.array([np.interp(temperature, rows["T:"], rows[str(z)])
                         for z in numbers])
    lower, upper = [interpolate_file(gname) for gname in filenames]
    return lower+(upper-lower)*(gravity-8.0)/0.5


def sample_metric(cache):
    # Buchan et al. 2025, Table 2 and Appendix A. The order follows
    # chemistry_info.usual_elements in the authors' public sample selector.
    row = next(r for r in read_csv(cache/"PEWDD_2024_raw.csv") if r[""]=="1828")
    assert row["Star"]=="WD 1232+563" and not row["log(Ti/H(e))"]
    solar = {"Ti":4.97,"Ca":6.30,"Ni":6.20,"Fe":7.46,"Mg":7.55,"Si":7.51,"O":8.69}
    atomic = {"Ti":22,"Ca":20,"Ni":28,"Fe":26,"Mg":12,"Si":14,"O":8}
    t,g,ca = float(row["T_eff"]),float(row["logg"]),measurement(row,"Ca")[0]
    files = ("tdif_He_no_overshoot_abedop1_800.csv", "tdif_He_no_overshoot_abedop1_850.csv")
    results = {}
    for label, order in (("catalog",("Ca","Ni","Fe","Mg","Si","O")),
                         ("corrected",("Ti","Ca","Fe","Mg","Si","O"))):
        numbers = [atomic[z] for z in order]
        tau = 10**np.array([warm_grid(cache,t,g,files,numbers),
                            cool_grid(cache,0,t,g,ca,numbers)])
        measured = [measurement(row, "Ni" if z=="Ti" else z) for z in order]
        values, errors = np.array(measured).T
        scaled = values - np.array([solar[z] for z in order])
        corrs = [np.corrcoef(times, scaled)[0,1] for times in tau]
        mg = tau[:,order.index("Mg")]
        depth = 5*np.mean(5*mg/(100000+5*mg))*max(0,np.mean(corrs))
        pairs=[]
        for i,j in combinations(range(len(order)),2):
            a,b=tau[:,i]/tau[:,j]
            score=abs(np.log10(a/b) + depth/np.log(10)*(a-b))/np.hypot(errors[i],errors[j])
            pairs.append((score,order[i]+"/"+order[j]))
        score,pair=max(pairs)
        results[label]={"metric":score, "limiting_pair":pair, "decline_depth":depth,
                        "test5_pass":bool(score>=2)}
    return results


def wd1232(cache, stars, meteorites):
    row = next(r for r in stars if r["Star"]=="WD 1232+563" and r["Paper"]=="Badenas-Agusti 2024")
    xu = next(r for r in stars if r["Star"]=="WD 1232+563" and r["Paper"]=="Xu 2019")
    nickel, nickel_err = measurement(row, "Ni")
    assert not row["log(Ti/H(e))"] and nickel==-8.92 and nickel_err==0.15
    assert measurement(xu, "Ni")[1]==-1  # PEWDD encodes an upper limit with -1.
    # Primary paper, Table 5: Ti=-8.92; Ni was not detected in this analysis.
    log_ca, log_ti, log_fe, log_ni = warm_grid(cache, float(row["T_eff"]), float(row["logg"]))
    reference = next(r for r in meteorites if r["Names"]=="Orgueil" and
                     all(r[f"[{z}/Mg]"] for z in ("Ni", "Fe", "Ti", "Ca")))
    def ratio(a,b): return float(reference[f"[{a}/Mg]"])-float(reference[f"[{b}/Mg]"])
    photocatalog = nickel-measurement(row,"Fe")[0]
    photopaper = nickel-measurement(row,"Ca")[0]
    return {
        "test5_published_selection_metric":sample_metric(cache),
        "catalog_apparent_log_Ni_Fe":photocatalog,
        "Orgueil_log_Ni_Fe":ratio("Ni","Fe"),
        "Ni_deficit_factor_build_up":10**(ratio("Ni","Fe")-photocatalog),
        "Ni_deficit_factor_steady":10**(ratio("Ni","Fe")-(photocatalog-(log_ni-log_fe))),
        "tau_Ni_over_tau_Fe":10**(log_ni-log_fe),
        "primary_log_Ti_Ca":photopaper,
        "primary_steady_log_Ti_Ca":photopaper-(log_ti-log_ca),
        "Orgueil_log_Ti_Ca":ratio("Ti","Ca"),
        "independent_Xu_log_Ti_Ca":measurement(xu,"Ti")[0]-measurement(xu,"Ca")[0],
        "Xu_2019_Ni_upper_limit":measurement(xu,"Ni")[0],
    }


def phase_fit(observed, log_taus, compositions, sigma):
    # tau in units of the Mg sinking time. Profile over dilution offset and
    # nonnegative time since switch-off analytically at each duration of supply.
    tau = 10**(log_taus-log_taus[1])
    durations = np.r_[0.0, np.geomspace(1e-5, 120, 321)]
    buildup = np.zeros((len(durations),len(tau)))
    buildup[1:] = np.log10(-np.expm1(-durations[1:,None]/tau)*tau)
    weight = 1/sigma**2
    centered = lambda x: x-(x*weight).sum(axis=-1,keepdims=True)/weight.sum()
    slope = centered(np.log10(np.e)/tau)
    chi = np.full(len(compositions), np.inf)
    endpoint = np.zeros((len(compositions),2))
    for duration, phase in zip(durations,buildup):
        residual = centered(observed-compositions-phase)
        since = np.maximum(0,-(residual*(slope*weight)).sum(axis=1)/(slope**2*weight).sum())
        score = ((residual+since[:,None]*slope)**2*weight).sum(axis=1)
        better = score<chi
        chi[better]=score[better]
        endpoint[better,0]=duration
        endpoint[better,1]=since[better]
    return chi,endpoint


def wd2236(cache, stars, meteorites):
    row = next(r for r in stars if (r["Star"],r["Paper"])==ROWS)
    observed, quoted = np.array([measurement(row,z) for z in ELEMENTS]).T
    rows = [r for r in meteorites if all(r[f"[{z}/Mg]"] for z in ELEMENTS)]
    compositions = np.array([[float(r[f"[{z}/Mg]"]) for z in ELEMENTS] for r in rows])
    result = {"samples":len(rows), "names":len({r["Names"] for r in rows}),
              "log_abundances":dict(zip(ELEMENTS,observed)), "fits":{}}
    for ov in (0,1):
        log_taus = cool_grid(cache,ov,float(row["T_eff"]),float(row["logg"]),observed[2])
        trials = [cool_grid(cache,ov,t,g,c)[:2] for t in (5150,5350,5550)
                  for g in (8.10,8.17,8.24) for c in (-9.40,-9.17,-8.90)]
        ratios = [10**(a-b) for a,b in trials]
        orgueil = [float(r["[Na/Mg]"]) for r in rows if r["Names"]=="Orgueil"]
        na_mg = observed[0]-observed[1]
        entry = {"log_tau_years":dict(zip(ELEMENTS,log_taus)), "systematic_floor":{},
                 "Na_Mg_bound":{"photospheric_log_Na_Mg":na_mg,
                    "tau_Na_over_tau_Mg":10**(log_taus[0]-log_taus[1]),
                    "parameter_stress_test_min_max":[min(ratios),max(ratios)],
                    "minimum_CI_depletion_factors":[10**(o-na_mg) for o in orgueil]}}
        sel = np.argwhere([r["Class"] in CHO for r in rows]).ravel()
        for floor in (0.0,0.1,0.15,0.2):
            sigma = np.hypot(quoted,floor)
            chi,end = phase_fit(observed,log_taus,compositions,sigma)
            bestall = int(np.argmin(chi))
            bestcho = int(sel[np.argmin(chi[sel])])
            indices = [bestall,bestcho]+[i for i,r in enumerate(rows) if r["Names"] in
                    ("Orgueil","earth_bulk_mcdonough","earth_silicate","chondrites","Roda")]
            summaries = []
            for i in dict.fromkeys(indices):
                summaries.append({"name":rows[i]["Names"],"class":rows[i]["Class"],
                    "chi2":round(float(chi[i]),4),"input_row":rows[i][""],
                    "duration_tauMg":round(float(end[i,0]),4),
                    "since_tauMg":round(float(end[i,1]),4)})
            entry["systematic_floor"][str(floor)]=summaries
        # Cross-check dependence on individual elements. Keeping five columns
        # and assigning a negligible fit weight also preserves the tauMg units.
        jack={}
        for i,z in enumerate(ELEMENTS):
            sigma=np.hypot(quoted,0.15);sigma[i]=1e12
            chi,_=phase_fit(observed,log_taus,compositions,sigma)
            j=int(np.argmin(chi)); jack[z]=rows[j]["Names"]
        entry["leave_one_out_best"]=jack
        result["fits"][str(ov)]=entry
    # Independent sensitivity to a tabulation with no calcium grid dimension.
    # A mass-conservation bound works even when tau_Na < tau_Mg.
    files=("Detlev_timescales_g8_He_overshoot.csv","Detlev_timescales_g8.5_He_overshoot.csv")
    ln,lm,logq=warm_grid(cache,float(row["T_eff"]),float(row["logg"]),files,(11,12,"qcvz:"))
    k=10**(lm-ln)-1
    mg_fraction=10**(logq+observed[1])*24.305/4.0026
    minimum_retention=mg_fraction**k
    reference=min(float(r["[Na/Mg]"]) for r in rows if r["Names"]=="Orgueil")
    offset=reference-(observed[0]-observed[1])
    result["legacy_grid_mass_limit"]={
        "tau_Na_over_tau_Mg":10**(ln-lm),"log_mixed_mass_fraction":logq,
        "present_Mg_mass_over_WD_mass":mg_fraction,"decay_exponent":k,
        "min_ratio_retention_allowing_one_WD_of_Mg":minimum_retention,
        "minimum_depletion_factor":10**offset*minimum_retention,
        "log10_Mg_mass_in_WDs_required_for_CI":np.log10(mg_fraction)+offset/k,
        "single_pulse_required_age_in_tauMg":offset/(k*np.log10(np.e)),
    }
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--cache",type=Path,default=HERE / ".cache")
    args=parser.parse_args()
    inputs(args.cache)
    stars=read_csv(args.cache/"PEWDD.csv")
    meteorites=read_csv(args.cache/"meteorite_database_Mg.csv")
    findings={"WD1232+563":wd1232(args.cache,stars,meteorites),
              "WDJ2236-5548":wd2236(args.cache,stars,meteorites)}
    print(json.dumps(findings,indent=2,allow_nan=False))
