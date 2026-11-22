"""Execute study stages on a runner; no inference depends on a missing archive."""
import argparse
import hashlib
import json
import os
import shutil
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

from prior_pairs import vectors


def load_components(out):
    paths = [out/"gaia_components.csv", out/"wide"/"gaia_components.csv"]
    return pd.concat([pd.read_csv(p, dtype={"system_id":str,"source_id":str,"field_id":str})
                      for p in paths if p.is_file()],ignore_index=True)


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


def pixel_evidence(systems, out):
    """Reuse only the same image coordinates, magnitudes and pixel-analysis revision."""
    import imaging
    digest = hashlib.sha256(Path(imaging.__file__).read_bytes()).hexdigest()
    previous = Path(__file__).parent/"evidence"
    historical = previous/"summary.json"
    cache, todo = {}, []
    def inputs(obj):
        return tuple(sorted((str(c["source_id"]), float(c["ra"]), float(c["dec"]), float(c["g"]))
                            for c in obj["components"]))
    if historical.exists() and os.environ.get("FRESH_PIXELS","").lower() != "true":
        metadata = json.loads(historical.read_text())
        old_digest = metadata.get("imaging_code_sha256")
        if old_digest == digest and (previous/"imaging"/"imaging.jsonl").exists():
            old_systems, _ = system_groups(load_components(previous))
            old_keys = {str(s["id"]):inputs(s) for s in old_systems}
            for line in (previous/"imaging"/"imaging.jsonl").read_text().splitlines():
                measurement = json.loads(line)
                key = old_keys.get(str(measurement["id"]))
                if key and "imaging_status" in measurement and "lens_model" in measurement:
                    cache[key] = measurement
    dest = out/"imaging"
    dest.mkdir(exist_ok=True)
    reused = []
    for obj in systems:
        prior = cache.get(inputs(obj))
        if prior is None:
            todo.append(obj)
            continue
        prior = dict(prior, reused_from_id=prior["id"], id=obj["id"])
        if prior.get("diagnostic_png"):
            shutil.copy2(previous/"imaging"/prior["diagnostic_png"], dest/prior["diagnostic_png"])
        reused.append(prior)
    print("Independent image sets:",len(systems),"reused:",len(reused),"requested:",len(todo),flush=True)
    new = imaging.assess(todo,dest,workers=6)
    complete = sorted(reused + new,key=lambda r:r["id"])
    with (dest/"imaging.jsonl").open("w") as f:
        for row in complete:
            f.write(json.dumps(row,allow_nan=False)+"\n")
    return complete

def inspect(out):
    import catalogues
    from prior_pairs import matches
    components = load_components(out)
    systems, aliases = system_groups(components)
    with ThreadPoolExecutor(2) as pool:
        remote = pool.submit(catalogues.annotate,components,out)
        images = pool.submit(pixel_evidence,systems,out)
        extra = matches(components,out)
        flags = remote.result()
        flags.to_csv(out/"system_catalogue_flags.csv",index=False)
        fits = {x["id"]:x for x in images.result()}
    prior_inputs = set(components.loc[components.input_lens_catalogue.map(flag),"system_id"])
    extra_ok = all(any(r.get("catalog")==tag and r.get("complete") for r in
                       json.loads((out/"supplemental_receipts.json").read_text()))
                   for tag in ["Chen et al. 2026 MGQPC", "SuGOHI VI", "SuGOHI current compilation"])
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
        from deflector import assess as check_deflector
        followup = []
        for obj in systems:
            if obj["id"] in short:
                followup.append(dict(obj,center_radec_deg=fits[obj["id"]]["lens_model"]["center_radec_deg"]))
        epochs.assess(followup, out/"epochs")
        check_deflector(followup, out/"deflectors")
    print("Prior to individual publication checks:",shortlist.to_string(index=False),flush=True)


def summarize_arm(out, label):
    folder = out if label == "primary" else out/"wide"
    counts = json.loads((folder/"gaia-counts.json").read_text())
    fields = pd.read_csv(folder/"gaia_fields.csv.gz",dtype={"field_id":str})
    all_decisions = pd.read_csv(out/"system_dispositions.csv")
    decisions = all_decisions[all_decisions.system_id.isin(fields.system_id)]
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
    bench.to_csv(folder/"known_recovery.csv",index=False)
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
    result["imaging_code_sha256"] = hashlib.sha256((Path(__file__).parent/"imaging.py").read_bytes()).hexdigest()
    (folder/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(label, json.dumps(result, indent=2))
    return result


def summarize(out):
    from audit import audit
    for arm in ("primary","wide"):
        folder = out if arm=="primary" else out/"wide"
        if not (folder/"gaia-counts.json").exists():
            continue
        summarize_arm(out,arm)
        measured = audit(folder,out)
        (folder/"audit.json").write_text(json.dumps(measured,sort_keys=True,indent=2)+"\n")
    packet(out)


def packet(out):
    """Publish self-contained records for plausible survivors and near misses."""
    components = load_components(out)
    flags = pd.read_csv(out/"system_catalogue_flags.csv").set_index("system_id")
    table = pd.read_csv(out/"system_dispositions.csv")
    spectra = pd.read_csv(out/"external_matches.csv",dtype=str)
    literature = pd.read_csv(out/"literature_matches.csv",dtype=str)
    extra = pd.read_csv(out/"supplemental_matches.csv",dtype=str)
    images = {r["id"]:r for r in (json.loads(line) for line in (out/"imaging"/"imaging.jsonl").read_text().splitlines())}
    for_review = table[table.cohort.eq("science") & table.full_morphology.map(flag)]
    result = []
    for r in for_review.itertuples():
        ident = r.system_id
        result.append(dict(system_id=ident,disposition=r.disposition,prior=flag(r.prior),
            components=components[components.system_id==ident].to_dict("records"), imaging=images[r.representative],
            catalogue_flags=flags.loc[ident].to_dict(), spectra=spectra[spectra.system_id==ident].to_dict("records"),
            lens_matches=literature[literature.system_id==ident].to_dict("records"),
            additional_matches=extra[extra.system_id==ident].to_dict("records")))
    with (out/"review_packets.jsonl").open("w") as handle:
        for row in result:
            # pandas JSON maps missing measurements to null, preserving text identifiers.
            handle.write(pd.Series(row).to_json(double_precision=14,force_ascii=False)+"\n")
    comp = components.groupby("system_id").source_id.apply(lambda x:" ".join(sorted(x.astype(str))))
    table["dr3_images"] = table.system_id.map(comp)
    table.to_csv(out/"system_dispositions.csv",index=False)
    unique_sets = table.groupby("dr3_images")
    counts = dict(gaia_unique_image_sets=len(unique_sets),records=len(table),
                  both_apertures=sum(any(g.system_id.str.match(r"^[FR]W")) and
                                     any(~g.system_id.str.match(r"^[FR]W")) for _,g in unique_sets),
                  unique_science_sets=table[table.cohort.eq("science")].dr3_images.nunique(),
                  unique_control_sets=table[table.cohort.eq("control")].dr3_images.nunique())
    (out/"image_set_overlap.json").write_text(json.dumps(counts,indent=2)+"\n")


if __name__ == "__main__":
    args = argparse.ArgumentParser()
    args.add_argument("stage",choices=["catalogues","inspect","summarize"])
    args.add_argument("output",type=Path)
    p = args.parse_args()
    p.output.mkdir(parents=True,exist_ok=True)
    {"catalogues":catalogue_stage,"inspect":inspect,"summarize":summarize}[p.stage](p.output)
