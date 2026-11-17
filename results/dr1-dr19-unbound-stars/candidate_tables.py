"""Publish a consistent reference distance for comparisons, plus every alternative and orbit."""

from make_queue import rows, save, number
from infer import load_inputs
import numpy as np

NOTES = {
    "1717410430709663616": (
        "primary",
        "One DESI exposure; five lines and two spectral fits",
        "Disk crossing; central origins disfavoured; systemic RV unmeasured",
    ),
    "4444095982221534720": (
        "secondary",
        "One DESI exposure; three usable lines",
        "Escape depends on MW mass and distance; central origin weakly sampled",
    ),
    "2752545725694798336": (
        "distance-sensitive",
        "Four DESI nights and one older SDSS epoch",
        "Disk-compatible; short dwarf distances yield bound orbits",
    ),
    "6340918677332110720": (
        "distance-sensitive",
        "Two DESI nights; SPDist outside training",
        "Remote plane crossing on long distance models; dwarf distance permits a bound halo orbit",
    ),
    "4435527217525195264": (
        "distance-sensitive",
        "Two DESI exposures on one night",
        "Dwarf photometry permits bound halo motion",
    ),
    "4550108286326235520": (
        "distance-sensitive",
        "One DESI exposure",
        "Dwarf photometry permits bound halo motion",
    ),
    "2294414169958935296": (
        "nearby-screen survivor",
        "One DESI exposure",
        "Escape changes with distance scale; ejection undetermined",
    ),
    "4406622293782568192": (
        "previously published positive energy",
        "Two BOSS nights; discrepant surface gravities",
        "MINESweeper already tabulated positive energy; shorter distance gives bound motion",
    ),
}

if __name__ == "__main__":
    stars, gs = load_inputs()
    ss = {s["source_id"]: s for s in stars}
    p = rows("inference.csv")
    o = rows("orbits.csv")
    detail = []
    for sid, (level, epochs, origin) in NOTES.items():
        s = ss[sid]
        g = gs[sid]
        r = next(x for x in p if x["source_id"] == sid and x["model"] == "bj_photogeo")
        orb = next(
            x for x in o if x["source_id"] == sid and x["model"] == "bj_photogeo"
        )
        detail.append(
            dict(
                source_id=sid,
                survey_id=s["targetid"],
                classification=level,
                rv_evidence=epochs,
                origin_assessment=origin,
                ruwe=g["ruwe"],
                pmra=g["pmra"],
                pmdec=g["pmdec"],
                pmra_error=g["pmra_error"],
                pmdec_error=g["pmdec_error"],
                parallax=g["parallax"],
                parallax_error=g["parallax_error"],
                parallax_zp=s["zp"],
                rv_sigma_with_floor=np.hypot(
                    number(s["e_rv"]), number(s.get("rv_floor") or 1)
                ),
                **{k: v for k, v in r.items() if k != "source_id"},
                **{k: v for k, v in orb.items() if k not in ("source_id", "model")},
            )
        )
    save("candidates.csv", detail)
    save("candidate_models.csv", [r for r in p if r["source_id"] in NOTES])
    save("candidate_origins.csv", [r for r in o if r["source_id"] in NOTES])
