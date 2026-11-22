"""Count unique configurations, known contaminants and the reasons a lead stops."""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from prior_pairs import vectors
from audit import rate


def summarize(evidence, output):
    decisions=pd.read_csv(evidence/"system_dispositions.csv",dtype={"dr3_images":str})
    images={r["id"]:r for r in (json.loads(l) for l in (evidence/"imaging"/"imaging.jsonl").read_text().splitlines())}
    targets=decisions[decisions.cohort=="science"].copy()
    targets["wide"] = targets.system_id.str.startswith("FW")
    complete=set(targets.loc[~targets.wide,"dr3_images"])
    wide=set(targets.loc[targets.wide,"dr3_images"])
    system_stats=[]
    for _,g in targets.groupby("dr3_images"):
        rows=g.sort_values("system_id")
        row=rows.iloc[0]
        record=images[row.representative]
        bands=record.get("bands",{})
        system_stats.append(dict(representative=row.representative,dr3_images=row.dr3_images,
            field_ids=" ".join(rows.system_id),passed_gaia=True,prior=bool(g.prior.any()),
            imaging_status=row.imaging_status,geometry=row.geometry,qso_anchor=bool(g.qso_anchor.any()),
            stellar_or_redshift_conflict=bool(g.stellar_or_redshift_conflict.any()),
            any_extended_band=any(b["extended_light"] for b in bands.values()),
            at_least_two_usable_filters=len({k.split("/")[-1] for k,v in bands.items() if v["usable_morphology"]})>=2,
            independent_colours=bool(record.get("consistent_two_band_photometry")),
            two_band_light=row.imaging_status=="extended_light",
            two_band_light_and_geometry=bool(row.full_morphology),
            favours_free_centre_over_each_single_host=bool(record.get("free_centroid_preferred_to_each_single_host_in_both_bands")),
            anchored_light_and_geometry=bool(g.plausible_before_novelty.any())))
    union=pd.DataFrame(system_stats)
    union.to_csv(output/"unique_science_configurations.csv",index=False)
    def count(f):
        return {k:int(f[k].sum()) for k in ["passed_gaia","qso_anchor","prior","stellar_or_redshift_conflict",
            "at_least_two_usable_filters","independent_colours","any_extended_band","two_band_light","geometry",
            "two_band_light_and_geometry","favours_free_centre_over_each_single_host","anchored_light_and_geometry"]}
    available=union[~union.prior]
    positive=available[available.anchored_light_and_geometry]
    positive.to_csv(output/"unlisted_two_band_shortlist.csv",index=False)
    result=dict(primary_only=len(complete-wide),wide_only=len(wide-complete),in_both=len(complete&wide),
                unique_science=count(union),unique_unlisted=count(available),
                unlisted_two_band_ids=positive.representative.tolist(),
                imaging_outcomes={str(k):int(v) for k,v in union.imaging_status.value_counts().items()})
    references=pd.read_csv(evidence/"known_lenses.csv")
    for arm,path in [("primary",evidence),("wider",evidence/"wide")]:
        fields=pd.read_csv(path/"gaia_fields.csv.gz",dtype={"field_id":str})
        science=fields[fields.cohort=="science"].reset_index(drop=True)
        tree=cKDTree(vectors(science.ra,science.dec))
        radius=2*np.sin(np.deg2rad(6/3600)/2)
        checks=[]
        for family in ["GLQ binaries","GLQ NIQs","GraL VI asterisms","GLQ doubles","GLQ quads"]:
            refs=references.loc[references.catalog==family].drop_duplicates(["name","ra","dec"])
            queries=tree.query_ball_point(vectors(refs.ra,refs.dec),radius)
            for ref,ix in zip(refs.itertuples(),queries):
                nearby=science.iloc[ix]
                passing=targets[targets.system_id.isin(nearby.system_id)]
                checks.append(dict(family=family,name=ref.name,ra=ref.ra,dec=ref.dec,
                    eligible=bool(len(ix)),gaia=bool(len(passing)),
                    morphology=bool(passing.full_morphology.any()),
                    anchored=bool(passing.plausible_before_novelty.any()),
                    fields=" ".join(nearby.system_id)))
        frame=pd.DataFrame(checks)
        frame.to_csv(output/(arm+"_published_comparisons.csv"),index=False)
        result[arm+"_published_comparisons"]={name:dict(catalogue_entries=len(g),eligible=int(g.eligible.sum()),
               gaia=int(g.gaia.sum()),morphology=int(g.morphology.sum()),anchored=int(g.anchored.sum()))
               for name,g in frame.groupby("family")}
    # Check zero triggering even in random fields for which quasar anchoring has no meaning.
    for arm,wide in [("primary",False),("wider",True)]:
        subset=decisions[(decisions.cohort=="control") & (decisions.system_id.str.startswith("RW")==wide)]
        fields=json.loads(((evidence/"wide" if wide else evidence)/"summary.json").read_text())["cohorts"]["control"]["fields"]
        result[arm+"_controls"]=dict(after_gaia=rate(len(subset),fields),
                light_and_model=rate(int(subset.full_morphology.sum()),fields))
    arms=[pd.read_csv(path/"known_recovery.csv") for path in (evidence,evidence/"wide")]
    benchmark=pd.concat(arms,ignore_index=True).groupby("name",as_index=False)[
        ["in_eligible_parent","passed_gaia","passed_light_and_model","passed_with_quasar_anchor"]].any()
    benchmark.to_csv(output/"known_recovery_either_aperture.csv",index=False)
    result["benchmark_either_aperture"]={"total":len(benchmark),**{k:int(benchmark[k].sum())
        for k in benchmark if k!="name"}}
    result["recovered_confirmed_lenses"]=benchmark.loc[benchmark.passed_with_quasar_anchor,"name"].tolist()
    (output/"combined_summary.json").write_text(json.dumps(result,indent=2)+"\n")


if __name__=="__main__":
    ev,out=map(Path,sys.argv[1:3])
    out.mkdir(parents=True,exist_ok=True)
    summarize(ev,out)
