"""Post-selection stress test for supersoft and absorbed spectra. No new source selection."""
import argparse
import json
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urljoin

import numpy as np
import pandas as pd
from astropy.io import fits
from bs4 import BeautifulSoup

from catalogs import photons,pspc_response
from full_only_review import get,plain
from prepare import truth

MPE="https://erosita.mpe.mpg.de/dr1/"
BASELINE="PL2_NH3e20"
MODELS={BASELINE:("powerlaw",2.,3e20)}
for kt in (.05,.08,.12):
    for nh in (3e20,1e21):
        MODELS[f"BB{kt:g}_NH{nh:.0e}"]=("blackbody",kt,nh)
for nh in (3e21,1e22):
    MODELS[f"PL1.5_NH{nh:.0e}"]=("powerlaw",1.5,nh)
KEV_TO_ERG=1.602176634e-9
X,W=np.polynomial.legendre.leggauss(6)


def spectrum(e,model):
    family,shape,nh=model
    e=np.asarray(e)
    transmission=photons(e.reshape(-1),0.,nh).reshape(e.shape)
    if family=="powerlaw":
        return (e**(-shape))*transmission
    # Photon (not energy) Planck spectrum; the arbitrary normalization cancels.
    return (e**2/np.expm1(np.minimum(e/shape,700.)))*transmission


def flux(lo,hi,model):
    grid=np.geomspace(lo,hi,8000)
    return KEV_TO_ERG*np.trapezoid(grid*spectrum(grid,model),grid)


def rate(response,model):
    lo,hi,area=response
    steps=(hi-lo)/2
    e=(hi+lo)[:,None]/2+steps[:,None]*X
    return np.sum(steps*area*np.sum(spectrum(e,model)*W,axis=1))


def ero_response(path,arf=None):
    # The same incident-energy grid is required. No effective-area interpolation
    # or bare RMF normalization is silently substituted for a full response.
    with fits.open(path) as h:
        m=next(x for x in h if "MATRIX" in list(getattr(getattr(x,"columns",None),"names",[]) or []))
        eb=h["EBOUNDS"].data
        channels=eb["CHANNEL"][((eb["E_MIN"]+eb["E_MAX"])/2>=.2)&((eb["E_MIN"]+eb["E_MAX"])/2<2.3)]
        lo,hi=m.data["ENERG_LO"].astype(float),m.data["ENERG_HI"].astype(float)
        area=[]
        for row in m.data:
            nf=int(row["N_GRP"])
            first=np.atleast_1d(row["F_CHAN"])[:nf]
            sizes=np.atleast_1d(row["N_CHAN"])[:nf]
            cc=np.concatenate([np.arange(int(j),int(j+n)) for j,n in zip(first,sizes)]) if nf else np.array([])
            pix=np.atleast_1d(row["MATRIX"])[:len(cc)]
            area.append(pix[np.isin(cc,channels)].sum())
        area=np.asarray(area,dtype=float)
    if arf:
        if np.nanmax(area)>1.05:
            raise ValueError("Redistribution matrix already includes collecting area")
        with fits.open(arf) as h:
            t=h["SPECRESP"].data
            if not (np.allclose(t["ENERG_LO"],lo) and np.allclose(t["ENERG_HI"],hi)):
                raise ValueError("Unmatched ARF and RMF energy grids")
            area*=t["SPECRESP"]
    if len(channels)==0 or np.nanmax(area)<10 or not np.isfinite(area).all() or np.any(area<0):
        raise ValueError("Invalid soft-band area or PI channels")
    return lo,hi,area


def find_ero(work,receipts,native_ecf):
    # Follow only MPE DR1 links. A canned single telescope or pointed-source
    # response is not a calibrated response for the main all-survey catalogue.
    pages=[MPE,MPE+"AllSkySurveyData_dr1/",MPE+"eSASS4DR1/",
           MPE+"AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Main_v1.2.html"]
    scanned=set()
    candidate_groups=[]
    errors=[]
    docs=[]
    def read(url):
        local=get(url,work,receipts,maxbytes=4_000_000)
        soup=BeautifulSoup(local.read_bytes(),"html.parser")
        related=[]
        products=[]
        calibration=[]
        for a in soup.select("a[href]"):
            href=urljoin(url,a["href"]).split("#",1)[0]
            if not href.startswith(MPE):
                continue
            label=" ".join(a.stripped_strings)
            context=" ".join(a.parent.stripped_strings)
            if re.search(r"\.(?:rsp|rmf|arf)(?:\.gz)?$",href,re.I):
                products.append((href,context[:1600]))
            elif re.search("response|calibrat|conversion factor|effective area",label+" "+href,re.I) and not re.search(r"\.(?:tar|tgz|gz|fits|pdf|zip)",href,re.I):
                related.append(href)
            if re.search(r"ECF|conversion factor|response matrices",label,re.I):
                calibration.append([label,href])
        for t in soup.select("table"):
            raw=" ".join(t.stripped_strings)
            if re.search(r"energy conversion factor|count.to.flux|\bECF\b",raw,re.I) and re.search(r"black.?body|supersoft|kT|photon index",raw,re.I):
                calibration.append(["Model-dependent ECF table; retained for row/model/unit verification",raw[:850]])
        return related,products,calibration
    for _ in range(3):
        todo=list(dict.fromkeys(p for p in pages if p not in scanned))[:8]
        if not todo:
            break
        pages=[]
        with ThreadPoolExecutor(max_workers=3) as pool:
            tasks=[(url,pool.submit(read,url)) for url in todo]
            for url,future in tasks:
                scanned.add(url)
                try:
                    related,products,calibration=future.result()
                    pages.extend(related)
                    docs.append({"page":url,"published_response_links":products,"ecf_references":calibration})
                    if products:
                        candidate_groups.append((url,products))
                except Exception as exc:
                    errors.append([url,str(exc)[:240]])
    for url,links in candidate_groups:
        def declared(context):
            return bool(re.search(r"erass1|dr1",context,re.I) and
                        re.search(r"survey[- ](?:averaged|average|weighted)|main (?:source )?catalog(?:ue)?|catalog(?:ue)? flux conversion",context,re.I) and
                        re.search(r"all (?:seven|7)|combined (?:tm|module)|tm0",context,re.I))
        valid=[link for link,context in links if declared(context+" "+link)]
        for primary in valid:
            arth=None
            if re.search(r"\.rmf(?:\.gz)?$",primary,re.I):
                arfs=[s for s in valid if re.search(r"\.arf(?:\.gz)?$",s,re.I)]
                if len(arfs)!=1 or sum(bool(re.search(r"\.rmf(?:\.gz)?$",s,re.I)) for s in valid)!=1:
                    continue
                arth=arfs[0]
            elif not re.search(r"\.rsp(?:\.gz)?$",primary,re.I):
                continue
            try:
                rsp=ero_response(get(primary,work,receipts),get(arth,work,receipts) if arth else None)
                test=flux(.2,2.3,MODELS[BASELINE])/rate(rsp,MODELS[BASELINE])
                if not abs(test/native_ecf-1)<=.05:
                    raise ValueError(f"Fiducial response ECF {test:g} differs from native {native_ecf:g} by >5%")
                return rsp,{"mode":"DR1 survey response folded","response":primary,"arf":arth,
                            "documentation":url,"native_broad_ECF":native_ecf,"folded_broad_ECF":test}
            except Exception as exc:
                errors.append([primary,str(exc)[:240]])
    return None,{"mode":"fixed released eROSITA count-to-flux calibration; absorbed band ratios only",
                 "native_broad_ECF":native_ecf,"scope":"No verified DR1 survey response or evaluated alternate-model ECF adopted. Published links restrict the search to the DR1 documentation examined.",
                 "documents":docs,"unavailable":errors}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    out=a.out
    work=a.work/"spectral-stress"
    work.mkdir(parents=True,exist_ok=True)
    receipts=[]
    native=json.loads((out/"spectral_model.json").read_text())
    inputs=json.loads((out/"catalog_inputs.json").read_text())["resources"]
    ref=native["rass_response"]
    path=get(ref,work,receipts)
    if receipts[-1]["sha256"]!=next(v["sha256"] for v in inputs.values() if v["url"]==ref):
        raise ValueError("PSPC-C response checksum differs from primary analysis")
    ros=pspc_response(path)
    baseline=MODELS[BASELINE]
    ref_ross=flux(.2,2,baseline)/rate(ros,baseline)
    lo,hi,area=ros
    midpoint_ecf=flux(.2,2,baseline)/np.sum(spectrum((hi+lo)/2,baseline)*(hi-lo)*area)
    if abs(midpoint_ecf/native["rass_ecf"]-1)>.002:
        raise ValueError("Midpoint PSPC-C fiducial does not reproduce the primary analysis")
    release=Path(__file__).parent/"release/ero_main.tar.gz.json"
    sample=json.loads(release.read_text())["example"]
    ero_native=float(sample["ML_FLUX_1"])/float(sample["ML_RATE_1"])
    ero,documentation=find_ero(work,receipts,ero_native)
    calibration=[]
    factors={}
    for label,model in MODELS.items():
        common=flux(.2,2,model)
        band_factor=common/flux(.2,2.3,model)
        ros_ecf=common/rate(ros,model)
        ul_relative=band_factor/(flux(.2,2,baseline)/flux(.2,2.3,baseline))
        if ero is None:
            er_relative=ul_relative
        else:
            er_relative=(common/rate(ero,model))/(flux(.2,2,baseline)/rate(ero,baseline))
        factors[label]={"xmm":1.,"ero":er_relative,"ero_ul":ul_relative,"rass":ros_ecf/ref_ross}
        calibration.append({"label":label,"family":model[0],"gamma_or_kT_keV":model[1],"NH_cm-2":model[2],
                            "ero_band_factor":band_factor,"rass_common_ECF":ros_ecf,"relative_factors":factors[label]})
    leads=pd.read_csv(out/"input.csv",usecols=["detection_uid","name"])
    comp=pd.read_csv(out/"comparisons.csv.gz",low_memory=False)
    comp=comp.loc[comp.detection_uid.isin(leads.detection_uid)&truth(comp.photometric_lead)].copy()
    results=[]
    for row in comp.itertuples():
        record={"detection_uid":row.detection_uid,"name":row.name,"high":row.high,"low":row.low,
                "low_kind":row.low_kind,"published_guarded_ratio":row.guarded_ratio,
                "published_original_grid_floor":row.model_floor_ratio}
        weak="ero_ul" if row.low=="ero" and row.low_kind=="upper_limit" else row.low
        for label,f in factors.items():
            record[label]=row.guarded_ratio*f[row.high]/f[weak]
        record["same_model_floor"]=min(record[name] for name in MODELS)
        record["independent_epoch_floor"]=row.guarded_ratio*min(f[row.high] for f in factors.values())/max(f[weak] for f in factors.values())
        results.append(record)
    results=pd.DataFrame(results)
    if len(results):
        audited=results.groupby("detection_uid").agg(independent_epoch_floor=("independent_epoch_floor","max"),
                  tied_model_floor=("same_model_floor","max"),published_original_grid_floor=("published_original_grid_floor","max"))
    else:
        audited=pd.DataFrame(columns=["independent_epoch_floor","tied_model_floor","published_original_grid_floor"]).rename_axis("detection_uid")
    if set(audited.index)!=set(leads.detection_uid):
        raise ValueError("Every previously qualified lead must retain an assessed comparison")
    results.to_csv(out/"spectral_stress_pairs.csv.gz",index=False,compression={"method":"gzip","mtime":0})
    leads.merge(audited,left_on="detection_uid",right_index=True,validate="one_to_one").to_csv(out/"spectral_stress_leads.csv",index=False)
    summary={"primary_photometric_leads":len(leads),"assessed_pairs":len(results),"original_grid_10x":int(audited.published_original_grid_floor.ge(10).sum()),
             "tied_stress_10x":int(audited.tied_model_floor.ge(10).sum()),"independent_epoch_stress_10x":int(audited.independent_epoch_floor.ge(10).sum()),
             "primary_ros_ECF":native["rass_ecf"],"recomputed_primary_midpoint_ECF":midpoint_ecf,"folded_baseline_ros_ECF":ref_ross,"ero_calibration":documentation,"models":calibration,
             "scope":"Post-selection, finite illustrative shape grid. Observed energy flux with Morrison & McCammon absorption. Blackbody photon density E^2/(exp(E/kT)-1); PL photon density E^(-Gamma). Neutral column held per model. Instrument band/source morphology, the XMM catalogue ECF, Poisson errors, and source/detection selection stay fixed. eROSITA uses only the declared calibration mode. No event classification, arbitrary-temperature bound, or corrected catalogue flux is inferred.",
             "ratios":"Multiply each existing count-based guarded amplitude by the ratio of the high and low conversion factors relative to its common fiducial. For an eROSITA map limit, vary only the band ratio at fixed released map ECF. RASS map limits use the refolded PSPC-C ECF. Both a common shape and independently allowed shapes across epochs are recorded. The lead statistic is the strongest surviving qualified pair.",
             "inputs":sorted(receipts,key=lambda item:item["url"])}
    (out/"spectral_stress.json").write_text(json.dumps(plain(summary),allow_nan=False,indent=2)+"\n")
    print("Supersoft/hard follow-up",len(leads),len(results),summary["independent_epoch_stress_10x"],flush=True)


if __name__=="__main__":
    main()
