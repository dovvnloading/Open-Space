"""Execute study stages on a runner; no inference depends on a missing archive."""
import argparse
import hashlib
import json
import os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from prior_pairs import vectors


def load_components(out):
    return pd.read_csv(out/"gaia_components.csv", dtype={"system_id":str,"source_id":str,"field_id":str})


def catalogue_stage(out):
    import catalogues
    from prior_pairs import extra_catalogues
    with ThreadPoolExecutor(2) as pool:
        jobs = [pool.submit(catalogues.load_lenses,out),pool.submit(extra_catalogues,out)]
        for job in jobs:
            job.result()


def flag(value):
    return str(value).lower() == "true"


def system_groups(components):
    systems = []
    seen = {}
    aliases = {}
    for ident, rows in components.groupby("system_id", sort=True):
        # Identical sets of DR3 images have identical pixels and lens geometry.
        # Keep every field in the denominator, including overlapping targeted apertures.
        key = tuple(sorted(rows.source_id.astype(str)))
        if key not in seen:
            seen[key] = ident
            by_source = rows.sort_values("source_id")
            systems.append(dict(id=ident, field_id=str(by_source.field_id.iloc[0]),
                                components=by_source.to_dict("records")))
        aliases[ident] = seen[key]
    return systems,aliases


def inspect(out):
    import catalogues
    import imaging
    from prior_pairs import matches
    components = load_components(out)
    systems, aliases = system_groups(components)
    with ThreadPoolExecutor(2) as pool:
        remote = pool.submit(catalogues.annotate,components,out)
        images = pool.submit(imaging.assess,systems,out/"imaging",workers=6)
        extra = matches(components,out)
        flags = remote.result()
        flags.to_csv(out/"system_catalogue_flags.csv",index=False)
        fits = {x["id"]:x for x in images.result()}
    prior_inputs = set(components.loc[components.input_lens_catalogue.map(flag),"system_id"])
    extra_ok = all(any(r.get("catalog")==tag and r.get("complete") for r in
                       json.loads((out/"supplemental_receipts.json").read_text()))
                   for tag in ["Chen et al. 2026 MGQPC", "SuGOHI VI"])
    flags["supplemental_prior"] = flags.system_id.isin(extra)
    flags["prior_input"] = flags.system_id.isin(prior_inputs)
    by_id = flags.set_index("system_id")
    rows = []
    for ident, representative in aliases.items():
        result, c = fits[representative], by_id.loc[ident]
        cohort = "science" if ident.startswith("F") else "control"
        prior = bool(c.prior_input or c.prior_match or c.supplemental_prior)
        rejected = bool(c.star_evidence or c.discordant_redshifts)
        plausible = bool(c.qso_evidence and not rejected and result["imaging_model_pass"])
        if prior:
            disposition = "previously_listed"
        elif rejected:
            disposition = "stellar_or_redshift_conflict"
        elif not c.qso_evidence:
            disposition = "no_independent_quasar_anchor"
        elif result["imaging_status"] in ("unavailable","inconclusive"):
            disposition = "inconclusive_imaging"
        elif not result["imaging_model_pass"]:
            disposition = "fails_light_or_geometry"
        elif not (c.lens_catalogues_complete and c.public_catalogues_complete and extra_ok):
            disposition = "publication_or_catalogue_coverage_incomplete"
        else:
            disposition = "individual_review_required"
        rows.append(dict(system_id=ident,cohort=cohort,representative=representative,
                         prior=prior,imaging_status=result["imaging_status"],
                         geometry=result["lens_model"]["physical_geometry"],
                         full_morphology=result["imaging_model_pass"],qso_anchor=c.qso_evidence,
                         stellar_or_redshift_conflict=rejected, plausible_before_novelty=plausible,
                         disposition=disposition))
    columns = ["system_id","cohort","representative","prior","imaging_status","geometry","full_morphology",
               "qso_anchor","stellar_or_redshift_conflict","plausible_before_novelty","disposition"]
    decisions = pd.DataFrame(rows,columns=columns)
    decisions.to_csv(out/"system_dispositions.csv",index=False)
    shortlist = decisions[(decisions.cohort=="science") & decisions.plausible_before_novelty.eq(True)
                          & decisions.prior.eq(False)]
    short = set(shortlist.representative)
    if short:
        import epochs
        epochs.assess([s for s in systems if s["id"] in short],out/"epochs")
    print("Prior to individual publication checks:",shortlist.to_string(index=False),flush=True)


def summarize(out):
    counts = json.loads((out/"gaia-counts.json").read_text())
    fields = pd.read_csv(out/"gaia_fields.csv.gz",dtype={"field_id":str})
    decisions = pd.read_csv(out/"system_dispositions.csv")
    known = pd.read_csv(out/"known_lenses.csv")
    bench = known[known.benchmark.map(flag)].copy()
    science = fields[fields.cohort=="science"].reset_index(drop=True)
    pos = vectors(science.ra,science.dec)
    tree = cKDTree(pos)
    rad = 2*np.sin(np.deg2rad(6/3600)/2)
    decision = decisions.set_index("system_id")
    benchmark = []
    for lens, nearby in zip(bench.itertuples(),tree.query_ball_point(vectors(bench.ra,bench.dec),rad)):
        subset = science.iloc[nearby]
        accepted = subset[subset.stage=="passed_gaia"]
        ids = list(accepted.system_id)
        detail = decision.reindex(ids).dropna(subset=["cohort"])
        benchmark.append(dict(name=lens.name,ra=lens.ra,dec=lens.dec,fields=" ".join(subset.system_id),
                              in_eligible_parent=bool(len(subset)),passed_gaia=bool(ids),
                              passed_light_and_model=bool(detail.full_morphology.map(flag).any()),
                              passed_with_quasar_anchor=bool(detail.plausible_before_novelty.map(flag).any())))
    bench = pd.DataFrame(benchmark,columns=["name","ra","dec","fields","in_eligible_parent","passed_gaia",
                                            "passed_light_and_model","passed_with_quasar_anchor"])
    bench.to_csv(out/"known_recovery.csv",index=False)
    count = lambda col:int(bench[col].map(flag).sum())
    result = dict(gaia=counts, benchmark=dict(total=len(bench),eligible=count("in_eligible_parent"),
                  gaia=count("passed_gaia"),morphology=count("passed_light_and_model"),
                  anchored=count("passed_with_quasar_anchor")), cohorts={})
    for cohort in ["science","control"]:
        subset = decisions[decisions.cohort==cohort]
        result["cohorts"][cohort] = dict(fields=int((fields.cohort==cohort).sum()),gaia=len(subset),
                unique_image_sets=subset.representative.nunique(), prior=int(subset.prior.map(flag).sum()),
                lens_geometry=int(subset.geometry.map(flag).sum()),light_and_geometry=int(subset.full_morphology.map(flag).sum()),
                anchored_with_light_and_geometry=int(subset.plausible_before_novelty.map(flag).sum()),
                dispositions={str(k):int(v) for k,v in subset.disposition.value_counts().items()},
                imaging={str(k):int(v) for k,v in subset.imaging_status.value_counts().items()})
    result["git_commit"] = os.getenv("GITHUB_SHA")
    result["run_id"] = os.getenv("GITHUB_RUN_ID")
    (out/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("stage",choices=["catalogues","inspect","summarize"])
    args.add_argument("output",type=Path)
    p = args.parse_args()
    p.output.mkdir(parents=True,exist_ok=True)
    {"catalogues":catalogue_stage,"inspect":inspect,"summarize":summarize}[p.stage](p.output)
