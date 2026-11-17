"""Estimate track-bank spacing at the actual six-detector search seasons."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.special import ndtri
from scipy.stats import qmc


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--orbit-dir',type=Path,required=True)
    args=parser.parse_args()
    spec=importlib.util.spec_from_file_location('orbit',args.orbit_dir/'propagate.py')
    orbit=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(orbit)
    inputs=json.loads((args.orbit_dir/'orbit_input.json').read_text())
    reference=orbit.family(inputs)
    sample=qmc.Sobol(7,scramble=True,seed=19832006).random_base2(15)
    pa=np.deg2rad(79.)
    def stress_offsets(sample):
        major,minor=ndtri(sample[:,3:5]).T
        iras=np.stack([18*major*np.sin(pa)+7*minor*np.cos(pa),
                       18*major*np.cos(pa)-7*minor*np.sin(pa)],axis=1)
        rr=40*np.sqrt(sample[:,5])
        aa=2*np.pi*sample[:,6]
        akari=np.stack([rr*np.cos(aa),rr*np.sin(aa)],axis=1)
        return [iras,akari]
    epochs=['2025-06-01','2025-06-20','2025-07-15','2025-11-25','2025-12-10','2025-12-24',
            '2026-06-01','2026-06-20','2026-07-15']
    def coordinates(states):
        pairs=[]
        for epoch in epochs:
            ra,de=orbit.sky(states,orbit.epoch(epoch))
            pairs.extend((ra*np.cos(np.deg2rad(-49.7)),de))
        return np.stack(pairs,axis=1)*3600
    bank=coordinates(reference)
    tree=cKDTree(bank)
    out={'epochs_utc':epochs,'grid':len(bank),'cases':{}}
    other=qmc.Sobol(7,scramble=True,seed=19832007).random_base2(15)
    stress=orbit.family(inputs,parameters=sample[:,:3],angular_offsets=stress_offsets(sample))
    cases=[
        ("nominal_offgrid", orbit.family(inputs,parameters=sample[:,:3]), bank, tree),
        ("historical_centroid_stress", stress, bank, tree),
    ]
    stress_bank=coordinates(stress)
    independent_stress=orbit.family(inputs,parameters=other[:,:3],angular_offsets=stress_offsets(other))
    cases.append(("independent_stress_to_searched_stress_bank",independent_stress,stress_bank,cKDTree(stress_bank)))
    for name,states,bank,tree in cases:
        projected=coordinates(states)
        rms,indices=tree.query(projected)
        delta=(projected-bank[indices]).reshape(len(projected),-1,2)
        offsets=np.linalg.norm(delta,axis=-1).max(axis=1)
        out['cases'][name]=dict(tracks=len(projected),max_epoch_offset_arcsec_percentiles=
            dict(zip(['p50','p90','p95','p99','max'],np.percentile(offsets,[50,90,95,99,100]))),
            fraction_within_half_native_pixel=float((offsets<3.075).mean()),
            fraction_within_one_native_pixel=float((offsets<6.15).mean()))
    (Path(__file__).resolve().parent/'orbit_spacing.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
