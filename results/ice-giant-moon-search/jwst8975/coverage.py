"""Select a field-coverage stratum using PSF support, independent of detected flux."""
import argparse
import csv
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import search


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube', type=Path)
    p.add_argument('bank', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('coverage_csv', type=Path)
    args = p.parse_args()
    cube, bank = search.load_npz(args.cube), search.load_npz(args.bank)
    kernels = search.validate_cube(cube)
    cuts = search.Cuts(psf_support=.8)
    maps = search.prepare(cube, kernels, cuts)
    cached = search.TemplateCache(cube, bank, kernels)
    _, info = cached.measure(maps, len(bank['x']), cuts)
    score = search.score(np.zeros_like(info), info, cube['visit'], cuts)
    use = score['eligible']
    keys = {'x','y','x_path','y_path','orbit_id','elements','a_km','phase_deg'}
    np.savez_compressed(args.output, **{k: (v[use] if k in keys else v) for k,v in bank.items()})
    with open(args.coverage_csv, 'w', newline='') as out:
        writer = csv.writer(out)
        writer.writerow(['orbit_id','a_km','phase_deg','eligible','exposures_visit1','exposures_visit2'])
        for i, flag in enumerate(use):
            writer.writerow([bank['orbit_id'][i],bank['a_km'][i],bank['phase_deg'][i],
                             bool(flag), *score['n'][i]])
    print('eligible', use.sum(), 'of', len(use))


if __name__ == '__main__':
    main()
