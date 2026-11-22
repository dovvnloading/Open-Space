"""Audit RASS-high source colours against the survey response and original 1RXS BSC.

Reads the completed comparison and host registers. All original sky positions,
identifiers, table units and fit metadata accompany the diagnostic model grid.
"""
import argparse
import hashlib
import io
import itertools
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd
import requests
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.votable import parse

from spectral_stress import flux, rate

TAP="https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync"
BANDS={"A":(11,41),"B":(52,201),"C":(52,90),"D":(91,201),"broad":(11,235)}
INPUTS=[]


def get(url,params=None):
    rec={"url":url,"params":params}
    try:
        r=requests.get(url,params=params,timeout=(12,100),headers={"User-Agent":"Open-Space RASS hardness provenance"})
        rec.update(status=r.status_code,final_url=r.url,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest())
        r.raise_for_status()
        return r.content
    except Exception as exc:
        rec["error"]=str(exc)[:280]
        raise
    finally:
        INPUTS.append(rec)


def query(sql):
    raw=get(TAP,dict(REQUEST="doQuery",LANG="ADQL",FORMAT="votable",MAXREC=50000,QUERY=sql))
    v=parse(io.BytesIO(raw))
    for node in list(v.infos)+[n for res in v.resources for n in res.infos]:
        if node.name=="QUERY_STATUS" and str(node.value).upper() in ("ERROR","OVERFLOW"):
            raise ValueError(str(node.content or node.value))
    table=v.get_first_table()
    INPUTS[-1]["returned_fields"]=[{"name":f.name or f.ID,"unit":str(f.unit) if f.unit is not None else "","ucd":str(f.ucd or ""),
                                    "datatype":str(f.datatype)} for f in table.fields]
    return table.to_table(use_names_over_ids=True).to_pandas()


def finite(value):
    try:
        n=float(value)
        return n if np.isfinite(n) else None
    except (ValueError,TypeError):
        return None


def namekey(value):
    return re.sub(r"\s+","",str(value).strip().upper())


def literal(name):
    return "'"+str(name).replace("'","''")+"'"


def native_schema(table):
    sql="SELECT column_name,unit,description,ucd FROM TAP_SCHEMA.columns WHERE table_name="+literal(table)
    return query(sql).fillna("")


def by_names(names,table):
    rows=[]
    for i in range(0,len(names),15):
        batch=names[i:i+15]
        sql=f"SELECT * FROM {table} WHERE name IN ("+",".join(map(literal,batch))+")"
        rows.extend(query(sql).to_dict("records"))
    return rows


def source_rows(out):
    leads=pd.read_csv(out/"input.csv",low_memory=False)
    pairs=pd.read_csv(out/"comparisons.csv.gz",low_memory=False)
    passed=pairs.photometric_lead.astype(str).str.lower().isin(["true","1","1.0"])
    high=pairs.loc[passed & pairs.high.eq("rass"),["detection_uid","low","low_kind","guarded_ratio","model_floor_ratio"]].copy()
    key=["detection_uid","name","xmm_time","rass_name","rass_ra","rass_dec","xmm_xtyp"]
    use=high.merge(leads[key],on="detection_uid",validate="many_to_one")
    v=out/"vetting"/"vetting.csv"
    if v.exists():
        host=pd.read_csv(v,usecols=["id","host_id","host_type"],low_memory=False)
        use=use.merge(host,left_on="detection_uid",right_on="id",how="left",validate="many_to_one")
    return use


def first_rass_catalog(targets,raw):
    objects={namekey(r["name"]):r for r in raw}
    coords=targets[["rass_name","rass_ra","rass_dec"]].drop_duplicates()
    rows=[]
    for i in range(0,len(coords),10):
        cones=[]
        for r in coords.iloc[i:i+10].itertuples(index=False):
            cones.append(f"1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{r.rass_ra:.9f},{r.rass_dec:.9f},{60/3600:.9f}))")
        rows.extend(query("SELECT * FROM rassbsc WHERE "+" OR ".join(cones)).to_dict("records"))
    others={namekey(r["name"]):r for r in rows}
    named=[str(v.get("onerxs_name","")).strip() for v in objects.values() if isinstance(v.get("onerxs_name"),str) and str(v["onerxs_name"]).strip()]
    if named:
        others.update({namekey(r["name"]):r for r in by_names(sorted(set(named)),"rassbsc")})
    matched=[]
    for r in coords.itertuples(index=False):
        cross_id=namekey(objects.get(namekey(r.rass_name),{}).get("onerxs_name",""))
        pos=SkyCoord(r.rass_ra,r.rass_dec,unit="deg")
        for k,entry in others.items():
            sep=pos.separation(SkyCoord(entry["ra"],entry["dec"],unit="deg")).arcsec
            if sep<=60 or k==cross_id:
                matched.append(dict(target_2rxs=r.rass_name,exact_released_1rxs_id=k==cross_id,
                                    sep_from_2rxs_arcsec=round(float(sep),4),**entry))
    return matched


def channel_responses(raw):
    with fits.open(io.BytesIO(raw),memmap=False) as hdus:
        tab=next(h for h in hdus if "MATRIX" in list(getattr(getattr(h,"columns",None),"names",[]) or [])).data
        lo,hi=tab["ENERG_LO"].astype(float),tab["ENERG_HI"].astype(float)
        area={k:np.zeros(len(lo)) for k in BANDS}
        for i,row in enumerate(tab):
            n=int(row["N_GRP"])
            first=np.atleast_1d(row["F_CHAN"])[:n].astype(int)
            count=np.atleast_1d(row["N_CHAN"])[:n].astype(int)
            channels=np.concatenate([np.arange(f,f+w) for f,w in zip(first,count)]) if n else np.empty(0)
            vals=np.atleast_1d(row["MATRIX"])[:len(channels)].astype(float)
            for k,(start,end) in BANDS.items():
                area[k][i]=vals[(channels>=start)&(channels<=end)].sum()
    return {k:(lo,hi,v) for k,v in area.items()}


def model_grid(out):
    m=json.loads((out/"spectral_model.json").read_text())
    url=m["rass_response"]
    inputs=json.loads((out/"catalog_inputs.json").read_text())["resources"]
    expected={d["sha256"] for d in inputs.values() if d["url"]==url}
    data=get(url)
    if hashlib.sha256(data).hexdigest() not in expected:
        raise ValueError("RASS response checksum differs from catalogue calibration")
    resp=channel_responses(data)
    shape=[("powerlaw",float(g)) for g in np.arange(1.,5.01,.25)]+[("blackbody",k) for k in (.03,.04,.05,.06,.08,.1,.12,.16,.2,.3)]
    nhs=(5e19,1e20,2e20,3e20,5e20,1e21,3e21,1e22)
    rows=[]
    for (family,param),nh in itertools.product(shape,nhs):
        model=(family,param,nh)
        cr={k:rate(r,model) for k,r in resp.items()}
        if min(cr.values())<=0:
            continue
        ecf=flux(.2,2.,model)/cr["broad"]
        rows.append(dict(family=family,shape=param,nh_cm2=nh,hr1=(cr["B"]-cr["A"])/(cr["B"]+cr["A"]),
                         hr2=(cr["D"]-cr["C"])/(cr["D"]+cr["C"]),ecf=ecf,relative_ecf=ecf/m["rass_ecf"]))
    grid=pd.DataFrame(rows)
    base=grid.loc[grid.family.eq("powerlaw") & grid["shape"].eq(2.) & grid.nh_cm2.eq(3e20)]
    if len(base)!=1 or abs(base.relative_ecf.iloc[0]-1)>.01:
        raise ValueError("fiducial response fold disagrees with primary count-to-flux conversion")
    return grid


def diagnostic(src,native,grid,bsc):
    rec={"detection_uid":src.detection_uid,"xmmsl2_name":src.name,"rass_name":src.rass_name,
         "low_mission":src.low,"low_kind":src.low_kind,"nominal_guarded_ratio":src.guarded_ratio,
         "archived_2rxs_exact_name":native is not None,"bsc_alternatives":len(bsc)}
    if native is None:
        return dict(rec,interpretation="no exact 2RXS record; no spectral inference")
    core=re.compile(r"^(?:hardness|band_[abcd]_(?:counts|tot_counts|bkg_counts)|plaw_|bb_(?:temperature|dof|chi2)|nh_gal|source_quality_flag|extended_region_flag|nearby_src_det_flag|count_rate|exposure|onerxs_(?:hr|name)|detection_likelihood|remaining_bkg_area)")
    rec.update({k:v for k,v in native.items() if core.search(k)})
    matched=[r for r in bsc if r.get("exact_released_1rxs_id")]
    if len(matched)==1:
        rec["exact_bsc_name"]=matched[0]["name"]
        rec.update({"bsc_"+k:v for k,v in matched[0].items() if re.search(r"hardness|^hr[_12]|count_rate|exposure",k)})
    fit=finite(native.get("plaw_photon_index"));err=finite(native.get("plaw_photon_index_error"))
    dof=finite(native.get("plaw_dof"))
    rec["plaw_fit_note"]=("finite catalogue estimate; inspect uncertainty, degrees of freedom, and reduced chi-square"
                          if fit and err and err>0 and dof and dof>0 else "no constrained positive-error power-law fit identified")
    allowed=np.ones(len(grid),dtype=bool)
    used=[]
    for i in (1,2):
        val=finite(native.get(f"hardness_ratio_{i}"));error=finite(native.get(f"hardness_ratio_{i}_error"))
        # Endpoints often have nominal zero error when a band is empty.
        # A symmetric Gaussian cut on an endpoint would spuriously force a spectrum.
        if val is not None and error is not None and -1<val<1 and error>0:
            used.append(f"HR{i}")
            allowed &= np.abs(grid[f"hr{i}"]-val)<=3*error
    rec["informative_colours"]=",".join(used)
    rec["n_grid_models"]=len(grid)
    rec["n_models_within_3sigma_marginals"]=int(allowed.sum())
    if used and allowed.any():
        subset=grid.loc[allowed]
        best=subset.loc[subset.relative_ecf.idxmin()]
        rec.update(min_colour_compatible_relative_ecf=float(best.relative_ecf),
                   min_ecf_model=f"{best.family} / {best['shape']:.3g} / NH {best.nh_cm2:.3g}",
                   max_colour_compatible_relative_ecf=float(subset.relative_ecf.max()),
                   min_rass_only_guarded_ratio=float(src.guarded_ratio*best.relative_ecf))
    if not used:
        rec["interpretation"]="no usable interior, positive-error hardness ratio; spectral ambiguity persists"
    elif not allowed.any():
        rec["interpretation"]="no grid model passes both usable colours; inspect data or a broader spectral family"
    else:
        rec["interpretation"]="conditional finite-grid colour diagnostic; low-state conversion fixed. Sparse bands, correlated colour errors, complex spectra and blends can change the result"
    return rec


def swift_countercheck(dest):
    from astropy.time import Time
    source="2SXPS J125851.4+235529"
    result={"xmmsl2":"XMMSL2 J125851.2+235522","requested_name":source,"table":"swift2sxps",
            "documentation":"https://heasarc.gsfc.nasa.gov/W3Browse/swift/swift2sxps.html"}
    try:
        schema=native_schema("swift2sxps")
        if not len(schema):
            raise ValueError("no native Swift table dictionary")
        schema.to_csv(dest/"swift2sxps_column_dictionary.csv",index=False)
        native=by_names([source],"swift2sxps")
        if len(native)!=1 or str(native[0].get("name","")).strip()!=source:
            raise ValueError(f"expected one exact Swift designation, received {len(native)}")
        row=native[0]
        pd.DataFrame([row]).to_csv(dest/"swift2sxps_native.csv",index=False)
        dictionary=schema.set_index("column_name").to_dict("index")
        observations={}
        mjd_dates={}
        for k,v in row.items():
            if not re.search(r"time|date|exposure|epoch|observation|snapshot|flux|rate|hardness|variab|var_|_var|prob|flag|band|spectrum|fit_",k,re.I):
                continue
            meta=dictionary.get(k,{})
            desc=str(meta.get("description",""))
            unit=str(meta.get("unit","")).strip()
            val=None if pd.isna(v) else v
            if isinstance(val,np.generic):
                val=val.item()
            if isinstance(val,bytes):
                val=val.decode("utf8",errors="replace")
            observations[k]={"value":val,"unit":unit,"description":desc}
            t=finite(v)
            if (t is not None and 44000<t<75000 and
                ("mjd" in (desc+" "+unit).lower() or
                 (k in ("time","start_time","end_time") and unit.lower() in ("d","day","days")))):
                mjd_dates[k]=Time(t,format="mjd",scale="utc").isot
        result.update(status="one_exact_row",returned_name=str(row["name"]),ra=row.get("ra"),dec=row.get("dec"),
                      mjd_fields_as_utc=mjd_dates,measurement_fields=observations,
                      interpretation="Native Swift band fluxes, count rates and variability probabilities retain each field's published definition. Source-level mean, peak or probability fields can span multiple snapshots or observations. A catalogue start or end date cannot date a particular peak unless explicitly linked to that measurement. No Swift flux is inserted into the three-survey soft-band ratio.")
        log={k:v["value"] for k,v in observations.items() if v["value"] is not None and
             re.search(r"^(?:time|start_time|end_time|exposure|det_flag|var_flag|count_rate|flux|.*prob.*var.*)$",k)}
        print("Swift exact source",source,"fields",log,"MJD times",mjd_dates,flush=True)
    except Exception as e:
        result.update(status="incomplete",error=str(e)[:300])
    (dest/"swift2sxps_review.json").write_text(json.dumps(result,indent=2,allow_nan=False,default=str)+"\n")
    return result["status"]


def main(out):
    dest=out/"rass_hardness"
    dest.mkdir(parents=True,exist_ok=True)
    targets=source_rows(out)
    targets.to_csv(dest/"target_pairs.csv",index=False)
    target_names=sorted({str(n) for n in targets.rass_name.dropna() if str(n).strip()})
    status={"selection":"All guarded pairs with RASS as the high epoch, retaining every host class.",
            "tables":{"rass2rxs":"Boller et al. 2016; 2RXS","rassbsc":"Voges et al. 1999; 1RXS bright catalogue"},
            "bands_pi":BANDS,
            "definitions":"Published HR1=(B-A)/(B+A), HR2=(D-C)/(D+C). The HEASARC band B counts description starts at PI 51, whereas its hardness-ratio definition starts at PI 52. The forward model follows the hardness definition. The native counts are preserved unchanged.",
            "limits":"RASS-modeled observed 0.2–2 keV flux is conditioned on a finite family of absorbed power laws (shape is Gamma) and blackbodies (shape is kT in keV). Three times the catalogue marginal error is a diagnostic tolerance, not a joint confidence region. One or two colours cannot uniquely resolve a soft excess, intrinsic absorption, blends, or mixtures. The low-state mission retains its published calibration. Both ROSAT catalogues reduce the same photons. Catalogued absorption-corrected PL, blackbody, and plasma-model fluxes are not used as observed common-band fluxes.",
            "urls":{"rass2rxs":"https://heasarc.gsfc.nasa.gov/W3Browse/rosat/rass2rxs.html","rassbsc":"https://heasarc.gsfc.nasa.gov/W3Browse/rosat/rassbsc.html"}}
    status["swift_2sxps_exact_countercheck"]=swift_countercheck(dest)
    for table in ("rass2rxs","rassbsc"):
        try:
            schema=native_schema(table)
            if not len(schema):
                raise ValueError("empty schema")
            schema.to_csv(dest/(table+"_column_dictionary.csv"),index=False)
        except Exception as e:
            status[table+"_schema_error"]=str(e)[:300]
    try:
        raw=by_names(target_names,"rass2rxs")
        status["unique_2rxs_queries"]=len(target_names)
        status["exact_2rxs_rows"]=len(raw)
    except Exception as e:
        raw=[]
        status["rass2rxs_query_error"]=str(e)[:300]
    native={}
    for item in raw:
        k=namekey(item["name"])
        if k in native:
            raise ValueError(f"non-unique 2RXS designation {k}")
        native[k]=item
    pd.DataFrame(raw if raw else [],columns=None if raw else ["name"]).to_csv(dest/"rass2rxs_native.csv",index=False)
    try:
        bsc=first_rass_catalog(targets,raw)
        status["bsc_rows_by_target"]=len(bsc)
    except Exception as e:
        bsc=[]
        status["rassbsc_query_error"]=str(e)[:300]
    pd.DataFrame(bsc if bsc else [],columns=None if bsc else ["target_2rxs","name"]).to_csv(dest/"rassbsc_native.csv",index=False)
    try:
        grid=model_grid(out)
        grid.to_csv(dest/"pspc_colour_grid.csv",index=False)
        check=[]
        for s in targets.itertuples(index=False):
            alt=[r for r in bsc if namekey(r["target_2rxs"])==namekey(s.rass_name)]
            check.append(diagnostic(s,native.get(namekey(s.rass_name)),grid,alt))
        pd.DataFrame(check).to_csv(dest/"source_diagnostics.csv",index=False)
        status["pair_diagnostics"]=len(check)
        report=[r for r in check if r["xmmsl2_name"]=="XMMSL2 J125851.2+235522"]
        compact=("rass_name","nominal_guarded_ratio","hardness_ratio_1","hardness_ratio_1_error","hardness_ratio_2","hardness_ratio_2_error",
                 "plaw_photon_index","plaw_photon_index_error","plaw_dof","informative_colours","n_models_within_3sigma_marginals",
                 "min_colour_compatible_relative_ecf","min_rass_only_guarded_ratio","interpretation")
        print("J125851.2+235522",[{k:r.get(k) for k in compact} for r in report],flush=True)
    except Exception as e:
        status["fold_error"]=str(e)[:300]
    status["receipts"]=INPUTS
    (dest/"status.json").write_text(json.dumps(status,indent=2)+"\n")
    print("2RXS colour audit",len(target_names),"guarded ROSAT-high sources;",len(raw),"exact 2RXS records;",len(bsc),"BSC counterparts",flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("out",type=Path)
    main(parser.parse_args().out)
