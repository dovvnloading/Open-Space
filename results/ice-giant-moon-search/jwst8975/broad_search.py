"""Screen the full inner-system grid in bounded-memory radial batches."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geometry
import search


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube', type=Path)
    p.add_argument('out', type=Path)
    p.add_argument('--radii', nargs=3, type=float, default=[52000,140000,500])
    p.add_argument('--phase-step', type=float, default=1)
    args = p.parse_args()
    cube = search.load_npz(args.cube)
    kernels = search.validate_cube(cube)
    cuts = search.Cuts(psf_support=.8)
    maps = search.prepare(cube, kernels, cuts)
    r0, r1, dr = args.radii
    if dr <= 0 or r0 > r1:
        raise ValueError('invalid radius grid')
    args.out.mkdir(parents=True, exist_ok=False)
    with open(args.out/'radial_bins.csv','x',newline='') as summary, open(args.out/'passing.csv','x',newline='') as hits:
        bins, peaks = csv.writer(summary), csv.writer(hits)
        bins.writerow(['a_km','templates','eligible','positive_passing','negative_passing','max_positive_rank','max_negative_rank'])
        peaks.writerow(['a_km','phase_deg','rank','total_snr','flux_Jy','visit_snr','dropout_snr','n'])
        with TemporaryDirectory() as tmp:
            fname = Path(tmp)/'tracks.npz'
            radii = np.arange(r0,r1+dr*1e-9,dr)
            for r in range(0,len(radii),4):
                group = radii[r:r+4]
                geometry.generate(args.cube, fname, grid=[group[0], group[-1], dr, args.phase_step])
                bank = search.load_npz(fname)
                search.validate_tracks(bank,cube)
                cache = search.TemplateCache(cube,bank,kernels)
                num, info = cache.measure(maps,len(bank['x']),cuts)
                obs = search.score(num,info,cube['visit'],cuts)
                neg = search.score(-num,info,cube['visit'],cuts)
                for a in group:
                    sel = bank['a_km']==a
                    pos, inv = obs['rank'][sel], neg['rank'][sel]
                    bins.writerow([a,sel.sum(),obs['eligible'][sel].sum(),(pos>=1).sum(),(inv>=1).sum(),pos.max(),inv.max()])
                for k in np.flatnonzero(obs['rank']>=1):
                    peaks.writerow([bank['a_km'][k],bank['phase_deg'][k],obs['rank'][k],
                                    obs['snr'][k],obs['flux'][k],*[json.dumps(obs[f][k].tolist())
                                    for f in ['visit_snr','dropout','n']]])
                summary.flush()
                hits.flush()
        (args.out/'parameters.json').write_text(json.dumps(dict(cube_sha256=search.file_hash(args.cube),
            search_sha256=search.file_hash(Path(search.__file__)), geometry_sha256=search.file_hash(Path(geometry.__file__)),
            script_sha256=search.file_hash(Path(__file__)), cuts=vars(cuts), radii=args.radii,
            phase_step=args.phase_step, epoch_tdb_jd=float(bank['epoch_jd_tdb']),
            exposure_id=cube['exposure_id'].tolist()),indent=2)+'\n')


if __name__=='__main__':
    main()
