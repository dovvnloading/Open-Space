"""Test a radial family of peaks against the azimuthal median projected ring halo."""
import argparse
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import numpy as np
from scipy.stats import binned_statistic

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geometry
import search


def subtract_ring_profile(cube):
    images = cube['image'].copy()
    yy, xx = np.indices(images.shape[1:])
    east_pix = -(xx-cube['center_xy'][0])
    north_pix = yy-cube['center_xy'][1]
    for k in range(len(images)):
        m = len(cube['path_ra'][k])//2
        los, east, north = geometry.tangent_basis(cube['path_ra'][k,m],cube['path_dec'][k,m])
        pole = geometry.pole(str(cube['planet']),cube['mjd'][k]+2400000.5)
        depth = -(east_pix*(pole@east)+north_pix*(pole@north))/(pole@los)
        radius = np.sqrt(east_pix**2+north_pix**2+depth**2)
        valid = cube['valid'][k] & np.isfinite(images[k])
        edges = np.arange(0,radius.max()+1,.65)
        med = binned_statistic(radius[valid], images[k][valid], 'median', bins=edges)[0]
        count = np.histogram(radius[valid],bins=edges)[0]
        good = np.isfinite(med) & (count>=20)
        mids = (edges[1:]+edges[:-1])/2
        model = np.interp(radius,mids[good],med[good])
        images[k] -= model
    return {**cube,'image':images}


def summarize(cube, bank, cuts, psfs):
    maps = search.prepare(cube, psfs, cuts)
    cache = search.TemplateCache(cube,bank,psfs)
    num, info = cache.measure(maps,len(bank['x']),cuts)
    stats = search.score(num,info,cube['visit'],cuts)
    k = np.argmax(stats['rank'])
    return dict(best_phase_deg=float(bank['phase_deg'][k]), templates_passing=int((stats['rank']>=1).sum()),
                best_rank=float(stats['rank'][k]), best_visit_snr=stats['visit_snr'][k].tolist(),
                nominal_phase_visit_snr=stats['visit_snr'][int(np.argmin(abs(bank['phase_deg']-102)))].tolist())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('cube', type=Path)
    p.add_argument('belinda', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    cube = search.load_npz(args.cube)
    psfs = search.validate_cube(cube)
    cuts = search.Cuts(psf_support=.8)
    with TemporaryDirectory() as tmp:
        f = Path(tmp)/'grid.npz'
        geometry.generate(args.cube,f,grid=[55000,55000,500,.5])
        bank = search.load_npz(f)
    report = dict(a_km=55000, cube_sha256=search.file_hash(args.cube),
                  original=summarize(cube,bank,cuts,psfs),
                  ring_subtracted=summarize(subtract_ring_profile(cube),bank,cuts,psfs))
    belinda = search.load_npz(args.belinda)
    j = np.flatnonzero(belinda['known_id']==714).item()
    coords = belinda['known_path'][j]
    body = dict(x=coords[None,:,coords.shape[1]//2,0], x_path=coords[None,:,:,0],
                y_path=coords[None,:,:,1], phase_deg=np.array([0]))
    report['belinda_ring_subtracted'] = summarize(subtract_ring_profile(belinda),body,cuts,psfs)
    args.output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    main()
