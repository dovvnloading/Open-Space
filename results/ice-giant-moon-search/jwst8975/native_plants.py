"""Plant in native detector coordinates, resample, and search the full real-image bank."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import warnings

from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
import numpy as np
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import geometry
import prepare
import search


class Detector:
    def __init__(self, path, record, cube, ee_path):
        if search.file_hash(path) != record['sha256']:
            raise ValueError(f'native FITS hash mismatch {path}')
        with fits.open(path, memmap=False) as h:
            sci = h['SCI']
            self.w = WCS(sci.header,naxis=2)
            self.shape = sci.data.shape
            self.native_scale = np.sqrt(sci.header['PIXAR_SR'])*geometry.ARCSEC
        self.ra, self.dec = self.w.all_pix2world(record['center_native_x'],record['center_native_y'],0)
        self.ref = cube['center_xy']
        self.scale = float(cube['scale_arcsec'])
        y,x = np.indices(cube['image'].shape[1:])
        dx,dy = self.native(x,y)
        self.coords = np.array([dy,dx])
        table = np.genfromtxt(ee_path, names=True)
        band = str(cube['band']).split('/')
        band = next(f for f in reversed(band) if f in table.dtype.names)
        edges = np.r_[0,table['aper_radius']]
        brightness = np.diff(np.r_[0,table[band]])/(np.pi*np.diff(edges**2))
        self.half = 7
        os = 8
        length = self.half*2+1
        xy = (np.arange(-self.half,self.half+1)[:,None]+(np.arange(os)+.5)/os-.5).ravel()
        self.phase_grid = np.zeros((8,8,length,length))
        for j in range(8):
            for i in range(8):
                rho = np.hypot(xy[:,None]-j/8,xy[None,:]-i/8)*self.native_scale
                z = brightness[np.searchsorted(edges[1:],rho)]*(self.native_scale/os)**2
                self.phase_grid[j,i] = z.reshape(length,os,length,os).sum((1,3))

    def native(self,x,y):
        e = -(x-self.ref[0])*self.scale/geometry.ARCSEC
        n = (y-self.ref[1])*self.scale/geometry.ARCSEC
        ra,dec = prepare.inverse_sky(self.ra,self.dec,e,n)
        return self.w.all_world2pix(ra,dec,0,quiet=True)

    def resampled_source(self,x,y,flux,jitter):
        nx,ny = self.native(x,y)
        # Registration challenge, independent zero-mean displacement per exposure.
        nx, ny = nx+jitter[0], ny+jitter[1]
        tx,ty = np.rint(nx*8).astype(int),np.rint(ny*8).astype(int)
        image = np.zeros(self.shape)
        weights = np.ones(len(tx))/(len(tx)-1)
        weights[[0,-1]] *= .5
        for x8,y8,w in zip(tx,ty,weights,strict=True):
            ix,iy = x8//8,y8//8
            stamp = self.phase_grid[y8%8,x8%8]
            patch = search.clipped_stamp(image.shape,iy-self.half,ix-self.half,stamp)
            if patch is not None:
                cut, part = patch
                image[cut] += part*(flux*w)
        # Input is Jy/native pixel; native-to-output area factor converts through MJy/sr.
        return ndimage.map_coordinates(image,self.coords,order=1,mode='constant')*(self.scale/self.native_scale)**2


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['cube','grid','covered','manifest','fits_dir','ee','out']:
        p.add_argument(name,type=Path)
    p.add_argument('--n',type=int,default=40)
    p.add_argument('--flux',type=float,nargs='+',default=[3e-6])
    p.add_argument('--seed',type=int,default=30517)
    p.add_argument('--jitter',type=float,default=.3,help='native pixel RMS per coordinate')
    a = p.parse_args()
    cube, grid, bank = [search.load_npz(x) for x in [a.cube,a.grid,a.covered]]
    if a.n > len(bank['x']):
        raise ValueError('too few geometrically eligible off-grid sources')
    kernels = search.validate_cube(cube)
    search.validate_tracks(grid,cube)
    search.validate_tracks(bank,cube)
    cuts = search.Cuts(psf_support=.8)
    records = json.loads(a.manifest.read_text())
    if [r['exposure_id'] for r in records] != cube['exposure_id'].tolist():
        raise ValueError('native manifest order differs from cube')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', FITSFixedWarning)
        detectors = [Detector(a.fits_dir/r['filename'],r,cube,a.ee) for r in records]
    cache = search.TemplateCache(cube,grid,kernels)
    base = search.scan(search.prepare(cube,kernels,cuts),grid,cube['visit'],cuts,cache=cache)
    basehits = np.flatnonzero(base['rank']>=1)
    rng = np.random.default_rng(a.seed)
    a.out.mkdir(parents=True,exist_ok=False)
    outrows = []
    visits = cube['visit']
    with open(a.out/'native_plants.csv','x',newline='') as f:
        fields = ['flux','orbit_id','a_km','phase_deg','recovered','status','recovered_orbit_id']
        rows = csv.DictWriter(f,fieldnames=fields)
        rows.writeheader()
        for flux in a.flux:
            for i in rng.choice(len(bank['x']),a.n,replace=False):
                row = dict(flux=flux,orbit_id=bank['orbit_id'][i],a_km=bank['a_km'][i],phase_deg=bank['phase_deg'][i],
                           recovered=False,status='missed',recovered_orbit_id='')
                hit = search.matched(basehits,bank['x'][i],bank['y'][i],grid,base,visits,cuts,1.5)
                if hit is not None:
                    row['status']='preexisting_peak'
                else:
                    image = cube['image'].copy()
                    jitter = rng.normal(0,a.jitter,(len(image),2))
                    for j,d in enumerate(detectors):
                        image[j] += d.resampled_source(bank['x_path'][i,j],bank['y_path'][i,j],flux,jitter[j])
                    maps = search.prepare({**cube,'image':image},kernels,cuts)
                    result = search.scan(maps,grid,visits,cuts,cache=cache)
                    peaks = np.flatnonzero(result['rank']>=1)
                    peaks = peaks[np.argsort(-result['rank'][peaks])]
                    match = search.matched(peaks,bank['x'][i],bank['y'][i],grid,result,visits,cuts,1.5)
                    if match is not None:
                        row.update(recovered=True,status='recovered',recovered_orbit_id=grid['orbit_id'][match])
                rows.writerow(row)
                f.flush()
                outrows.append(row)
    metadata = dict(seed=a.seed,jitter_sigma_px_per_axis=a.jitter,cuts=vars(cuts),
                    only_fixed_calibrated_noise=True,summary=search.summarize_plants(outrows),
                    sha256={str(k):search.file_hash(getattr(a,k)) for k in ['cube','grid','covered','manifest','ee']},
                    script_sha256=search.file_hash(Path(__file__)))
    (a.out/'report.json').write_text(json.dumps(metadata,indent=2)+'\n')


if __name__=='__main__':
    main()
