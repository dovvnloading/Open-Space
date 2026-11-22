"""Audit 1RXS BSC and FSC as earlier reductions of the same RASS observations."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord, search_around_sky

from catalogs import Archive, NASA, get_tap_table, get_tdat, numeric

TABLES = {"rassbsc":(60.,18000,19000), "rassfsc":(90.,105000,106000)}
COLUMNS = ["detection_uid","xmm_name","xmmsl2_legacy_hint","rass_ul_pair","rass_ul_high_missions",
           "catalogue","legacy_name","legacy_ra","legacy_dec","legacy_poserr_arcsec",
           "xmm_separation_arcsec","ero_separation_arcsec","search_radius_arcsec",
           "plausible","n_alternatives","multiple","ambiguous","needs_review","status",
           "rate_cps","rate_error_cps","detection_likelihood","flux_0p2_2_cgs","flux_error_cgs",
           "high_mission_for_ratio","diagnostic_high_flux","diagnostic_high_flux_error",
           "diagnostic_ratio_high_over_legacy","diagnostic_guarded_ratio_high_over_legacy"]


class Audit(Archive):
    def __init__(self, work, out):
        super().__init__(work,out)
        self.summary = {"catalogues":{}, "meaning":"1RXS and 2RXS reduce the same ROSAT survey epoch. A historical counterpart is never counted as a fourth epoch."}

    def save(self):
        result = dict(self.summary,inputs=self.records,unsuccessful_routes=self.failures)
        (self.out/"rass_legacy_status.json").write_text(json.dumps(result,indent=2)+"\n")


def read_legacy(arc,table):
    for reader,args in [(get_tap_table,(NASA,table,"rass")),(get_tdat,(table,"rass"))]:
        try:
            data,mapping,units,desc = reader(arc,*args)
            # The TAP reader separately checks its result against SELECT COUNT(*).
            # The release-scale bounds catch empty/truncated historical mirrors.
            _,minimum,maximum = TABLES[table]
            if not minimum<=len(data)<=maximum:
                raise ValueError(f"Unexpected historical release size {len(data)} for {table}")
            data["ra"],data["dec"] = numeric(data,"ra"),numeric(data,"dec")
            data["rate"],data["errrate"] = numeric(data,"rate"),numeric(data,"errrate")
            data["ml"],data["pos"] = numeric(data,"ml"),numeric(data,"pos")
            if "pos" in mapping:
                # HEASARC 1RXS error_radius is released in arcsec; normalize if
                # the TAP metadata declares a different angular unit.
                unit = str(units.get(mapping["pos"],"")).lower()
                if "arcmin" in unit:
                    data["pos"] *= 60.
                elif "deg" in unit:
                    data["pos"] *= 3600.
            valid = data.ra.between(0,360,inclusive="left") & data.dec.between(-90,90)
            if not valid.all():
                raise ValueError("invalid 1RXS coordinates")
            arc.summary["catalogues"][table] = {"status":"ok","rows":len(data),
                    "reader":reader.__name__, "location":args[:-1],"mapping":mapping,
                    "units":{k:str(units.get(v,"")) for k,v in mapping.items()},
                    "legacy_schema":{str(k):str(v) for k,v in desc.items()}}
            arc.save()
            return data
        except Exception as exc:
            arc.note(table+" "+reader.__name__,exc)
    arc.summary["catalogues"][table] = {"status":"unavailable"}
    arc.save()
    return None


def bool_series(s):
    return s.astype(str).str.lower().isin(["true","1","1.0","t","yes"])


def extract_leads(out):
    sources = pd.read_csv(out/"input.csv",low_memory=False)
    if sources.detection_uid.duplicated().any():
        raise ValueError("duplicate legacy-audit detection IDs")
    pairs = pd.read_csv(out/"comparisons.csv.gz",usecols=["detection_uid","low","high","low_kind","photometric_lead"])
    candidates = pairs.loc[bool_series(pairs.photometric_lead) & pairs.low.eq("rass") &
                           pairs.low_kind.eq("upper_limit")]
    high = candidates.groupby("detection_uid").high.agg(lambda vals:"|".join(sorted(set(vals))))
    sources["rass_ul_high_missions"] = sources.detection_uid.map(high).fillna("")
    sources["rass_ul_pair"] = sources.rass_ul_high_missions.ne("")
    return sources


def find_matches(sources,cat,table):
    limit = TABLES[table][0]
    x = SkyCoord(sources.xmm_ra.to_numpy()*u.deg,sources.xmm_dec.to_numpy()*u.deg)
    c = SkyCoord(cat.ra.to_numpy()*u.deg,cat.dec.to_numpy()*u.deg)
    ids,targets,_,_ = search_around_sky(x,c,limit*u.arcsec)
    pairs = set(zip(ids,targets))
    valid = (bool_series(sources.ero_quality)&~bool_series(sources.ero_ambiguous)&
             sources.ero_ra.between(0,360)&sources.ero_dec.between(-90,90))
    er = np.flatnonzero(valid)
    ecoords = SkyCoord(sources.ero_ra.iloc[er].to_numpy()*u.deg,sources.ero_dec.iloc[er].to_numpy()*u.deg)
    if len(er):
        ii,jj,_,_ = search_around_sky(ecoords,c,limit*u.arcsec)
        pairs.update(zip(er[ii],jj))
    lookup = {i:j for j,i in enumerate(er)}
    for i,j in sorted(pairs):
        legacy = cat.iloc[j]
        # A catalogue source anywhere within either historical search radius is
        # kept as a plausible alternative. Wider ROSAT positions remain reviewable.
        yield i, {"catalogue":table,"legacy_name":str(legacy["name"]).strip(),
                  "legacy_ra":legacy.ra,"legacy_dec":legacy.dec,"legacy_poserr_arcsec":legacy.pos,
                  "xmm_separation_arcsec":x[i].separation(c[j]).arcsec,
                  "ero_separation_arcsec":ecoords[lookup[i]].separation(c[j]).arcsec if i in lookup else np.nan,
                  "search_radius_arcsec":limit,"rate_cps":legacy.rate,"rate_error_cps":legacy.errrate,
                  "detection_likelihood":legacy.ml,"plausible":True}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--work",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    arc=Audit(a.work,a.out)
    sources=extract_leads(a.out)
    model=json.loads((a.out/"spectral_model.json").read_text())
    ecf=model["rass_ecf"]
    matches = {i:[] for i in range(len(sources))}
    if len(sources):
        for table in TABLES:
            data=read_legacy(arc,table)
            if data is not None:
                for i,hit in find_matches(sources,data,table):
                    matches[i].append(hit)
    complete = all(arc.summary["catalogues"].get(t,{}).get("status")=="ok" for t in TABLES) if len(sources) else True
    output=[]
    for i,row in sources.iterrows():
        assoc = matches[i]
        alternatives = len(assoc)
        # Count separate catalogue entries, even if the nominal centroids overlap.
        # This is deliberately conservative for ambiguous historical assignments.
        blank = {"catalogue":"","legacy_name":"","plausible":False}
        for hit in assoc or [blank]:
            status = "plausible historical detection" if alternatives else ("clear within radii" if complete else "coverage incomplete")
            ambiguous = alternatives>1 or not complete
            entry = dict(hit,detection_uid=row.detection_uid,xmm_name=row["name"],
                         xmmsl2_legacy_hint=row.get("xmm_oldrass"),rass_ul_pair=row.rass_ul_pair,
                         rass_ul_high_missions=row.rass_ul_high_missions,n_alternatives=alternatives,
                         multiple=alternatives>1,ambiguous=ambiguous,
                         needs_review=bool(row.rass_ul_pair and (alternatives or not complete)),status=status)
            if alternatives and ecf and hit["rate_cps"]>0 and hit["rate_error_cps"]>0:
                f,err=hit["rate_cps"]*ecf,hit["rate_error_cps"]*ecf
                entry.update(flux_0p2_2_cgs=f,flux_error_cgs=err)
                for high in row.rass_ul_high_missions.split("|"):
                    if high not in ("xmm","ero"):
                        continue
                    hf,he=row[high+"_flux"],row[high+"_err"]
                    record = dict(entry,high_mission_for_ratio=high,diagnostic_high_flux=hf,
                                  diagnostic_high_flux_error=he,diagnostic_ratio_high_over_legacy=hf/f,
                                  diagnostic_guarded_ratio_high_over_legacy=max(0.,hf-1.96*he)/(f+1.96*err))
                    output.append(record)
                if row.rass_ul_pair:
                    continue
            output.append(entry)
    audit = pd.DataFrame(output,columns=COLUMNS)
    audit.to_csv(a.out/"rass_legacy.csv",index=False,float_format="%.10g")
    arc.summary.update(total_leads=len(sources),completed_both=complete,
                       n_leads_with_plausible_matches=sum(bool(v) for v in matches.values()),
                       n_affected_rass_ul_leads=int(audit.loc[audit.needs_review.eq(True),"detection_uid"].nunique()),
                       n_leads_with_multiple_matches=sum(len(v)>1 for v in matches.values()),
                       radii_arcsec={k:v[0] for k,v in TABLES.items()},rass_ecf_cgs_per_cps=ecf,
                       definition="All candidates in 60 arcsec (BSC) or 90 arcsec (FSC) about a slew position or a reliable eROSITA centroid are retained. These radii are the historical XMMSL2 User Guide §6 RASS match radii. Fluxes use the same PSPC-C response as the primary analysis. Ratios are diagnostic and do not substitute automatically for the map or the 2RXS entries.")
    arc.save()
    print(json.dumps({k:v for k,v in arc.summary.items() if k.startswith("n_") or k=="completed_both"}),flush=True)


if __name__=="__main__":
    main()
