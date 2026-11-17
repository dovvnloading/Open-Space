"""Forced trailed photometry of one named, deliberately unmasked moon."""
import argparse
import json
from pathlib import Path
import numpy as np
from search import Cuts, file_hash, load_npz, prepare, sample_paths, score, snr, validate_cube


def run(path, body, dest, support):
    cube=load_npz(path)
    if str(cube["unmasked_control_body"])!=str(body):
        raise ValueError("control cube must have just the selected body unmasked")
    cut=Cuts(psf_support=support)
    maps=prepare(cube,validate_cube(cube),cut)
    xy=cube["known_path"][list(cube["known_id"]).index(body)]
    # Translate the entire predicted track. These local controls assess the
    # same region of the images; overlapping positions are correlated.
    offsets=[(0.,0.)]
    for radius in [8.,12.,16.]:
        for a in np.arange(0,360,15):
            offsets.append((np.cos(np.deg2rad(a))*radius,np.sin(np.deg2rad(a))*radius))
    dx,dy=np.array(offsets).T
    num,info=sample_paths(maps,xy[None,:,:,0]+dx[:,None,None],xy[None,:,:,1]+dy[:,None,None],cut)
    fit=score(num,info,cube["visit"],cut)
    match=np.all(fit["n"][1:]==fit["n"][0],axis=1)
    report=dict(moon=body,cube_sha256=file_hash(path),psf_support=support,
                exposure_id=cube["exposure_id"].tolist(),visits=np.unique(cube["visit"]).tolist(),
                supported=bool(fit["eligible"][0]),n_exposures=fit["n"][0].tolist(),
                formal_snr_per_frame=snr(num[0],info[0]).tolist(),visit_scores=fit["visit_snr"][0].tolist(),
                dropout_scores=fit["dropout"][0].tolist(),pooled_score=float(fit["snr"][0]),
                fitted_jy=float(fit["flux"][0]) if fit["eligible"][0] else None,
                two_visit_cuts_passed=bool(fit["rank"][0]>=1),
                matched_coverage_local_controls=int(match.sum()),
                local_controls_at_least_as_high=int(np.sum(fit["snr"][1:][match]>=fit["snr"][0])))
    Path(dest).write_text(json.dumps(report,indent=2)+"\n")
    print(body,report["n_exposures"],report["visit_scores"],report["pooled_score"])


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("cube")
    p.add_argument("body",type=int)
    p.add_argument("output")
    p.add_argument("--support",type=float,default=.90)
    a=p.parse_args()
    run(a.cube,a.body,a.output,a.support)
