"""Measure a held-out known moon at its spacecraft ephemeris and displaced controls."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import search


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('control_cube', type=Path)
    p.add_argument('science_cube', type=Path, help='reference PSF and enclosed-flux fraction')
    p.add_argument('output', type=Path)
    args = p.parse_args()
    cube = search.load_npz(args.control_cube)
    reference = search.load_npz(args.science_cube)
    if not np.array_equal(cube['exposure_id'], reference['exposure_id']):
        raise ValueError('reference PSF does not match the exposures')
    cube.update(psf=reference['psf'], psf_fraction=reference['psf_fraction'])
    body = int(cube['unmasked_control_body'])
    track = cube['known_path'][np.flatnonzero(cube['known_id']==body).item()]
    cuts = search.Cuts(psf_support=.8)
    maps = search.prepare(cube, search.validate_cube(cube), cuts)
    offsets = [(0, 0)]+[(rad*np.cos(a), rad*np.sin(a)) for rad in [7, 11, 15]
                        for a in np.arange(12)*np.pi/6]
    px = np.stack([track[:, :, 0]+x for x, y in offsets])
    py = np.stack([track[:, :, 1]+y for x, y in offsets])
    num, info = search.sample_paths(maps, px, py, cuts)
    scores = search.score(num, info, cube['visit'], cuts)
    report = dict(body=body, exposure_id=cube['exposure_id'].tolist(), visits=np.unique(cube['visit']).tolist(),
                  per_exposure_snr=search.snr(num, info)[0].tolist(), visit_snr=scores['visit_snr'][0].tolist(),
                  exposures_per_visit=scores['n'][0].tolist(), dropout=scores['dropout'][0].tolist(),
                  full_rank=float(scores['rank'][0]), offsets_px=offsets, displaced_rank=scores['rank'][1:].tolist(),
                  cube_sha256=search.file_hash(args.control_cube), science_psf_cube_sha256=search.file_hash(args.science_cube))
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(body,report['per_exposure_snr'],report['full_rank'])


if __name__ == '__main__':
    main()
