"""Blind resolved-population search and family-wise permutation calibration."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import fft, ndimage

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'measurement'
STEP = 1.5 / 60
RADII = np.array([3., 6., 12.])
DM = np.arange(16., 23.01, .5)
FIELDS = {'bootes': (209, 212, 13.5, 15.5), 'coma': (185.75, 188.75, 22.9, 24.9)}
# Reference centres and moduli from McConnachie 2012, Table 1.
CONTROLS = {'bootes': (210.025, 14.5, 19.11), 'coma': (186.7458, 23.904, 18.2)}
SEED = 20261117
NNULL = 99


def truth(series):
    return series.astype(str).str.lower().isin(['true', 't', '1', '1.0', 'yes'])


def catalogue(field):
    files = list((ROOT / 'catalog').glob(f'{field}.csv.gz'))
    if len(files) != 1:
        raise ValueError(f'Expected catalog/{field}.csv.gz')
    data = pd.read_csv(files[0])
    keep = truth(data.star) & truth(data.clean) & data.g.between(17,23.5) & data.r.between(16,23.5)
    for col in ('gerr','rerr'):
        keep &= (data[col] > 0) & (data[col] < .2)
    return data.loc[keep].copy()


def weights(g, r, gerr, rerr, moduli=DM):
    """Gaussian colour likelihood along the pre-He-flash isochrone."""
    track = np.genfromtxt(ROOT / 'model/parsec.csv', names=True, delimiter=',')
    seq = track[track['stage'] == 0]
    seq = seq[np.argsort(seq['g'])]
    colour = seq['g'] - seq['r']
    sigma = np.sqrt(.05**2 + gerr**2 + rerr**2)
    locus = np.interp(g - moduli[:,None], seq['g'], colour, left=np.nan, right=np.nan)
    return np.nan_to_num(np.exp(-.5*((g-r-locus)/sigma)**2), nan=0.)


class FieldGrid:
    def __init__(self, field):
        self.ra0, ra1, self.dec0, dec1 = FIELDS[field]
        self.cos = np.cos(np.deg2rad((self.dec0+dec1)/2))
        self.nx = int(np.floor((ra1-self.ra0)*self.cos/STEP))
        self.ny = int(np.floor((dec1-self.dec0)/STEP + 1e-9))
        self.shape = (self.ny, self.nx)
        self.y,self.x = np.indices(self.shape)
        dx = np.minimum(self.x, self.nx-self.x)
        dy = np.minimum(self.y, self.ny-self.y)
        d = np.hypot(dx, dy)*STEP*60
        self.kernels, self.valid = [], []
        for rad in RADII:
            on = d <= rad + 1e-8
            off = (d > 2*rad + 1e-8) & (d <= 3*rad + 1e-8)
            alpha = on.sum()/off.sum()
            diff = fft.rfft2(on.astype(float)-alpha*off)
            var = fft.rfft2(on.astype(float)+alpha**2*off)
            self.kernels.append((diff, var))
            p = int(np.ceil(3*rad/(60*STEP)))
            self.valid.append((self.x>=p)&(self.x<self.nx-p)&(self.y>=p)&(self.y<self.ny-p))
        self.valid = np.array(self.valid)

    def pixels(self, ra, dec):
        ix = np.floor((ra-self.ra0)*self.cos/STEP).astype(int)
        iy = np.floor((dec-self.dec0)/STEP).astype(int)
        good = (ix>=0)&(ix<self.nx)&(iy>=0)&(iy<self.ny)
        return (iy*self.nx+ix)[good], good

    def coords(self, y, x):
        return self.ra0+(x+.5)*STEP/self.cos, self.dec0+(y+.5)*STEP

    def dist(self, ra, dec, ra0, dec0):
        return np.hypot((ra-ra0)*self.cos,dec-dec0)*60

    def moments(self, pix, weight):
        def spatial(w):
            binned = np.stack([np.bincount(pix, weights=row, minlength=self.nx*self.ny)
                               for row in w]).reshape((-1,*self.shape))
            return fft.rfft2(binned)
        wfft, vfft = spatial(weight), spatial(weight**2)
        diff = np.array([fft.irfft2(wfft*k[0],s=self.shape) for k in self.kernels])
        var = np.array([fft.irfft2(vfft*k[1],s=self.shape) for k in self.kernels])
        return diff, var

    def scores(self, diff, var):
        result = diff/np.sqrt(1+np.maximum(0,var))
        return np.where(self.valid[:,None],result,-np.inf)


def peaks(grid, scores, floor=-np.inf):
    best = scores.max(axis=(0,1))
    yy,xx = np.where((best==ndimage.maximum_filter(best,size=3))&np.isfinite(best)&(best>floor))
    ranked = np.argsort(-best[yy,xx],kind='stable')
    ra,dec = grid.coords(yy,xx)
    rows = []
    for i in ranked:
        k,j = np.unravel_index(np.argmax(scores[:,:,yy[i],xx[i]]),scores.shape[:2])
        if any(grid.dist(ra[i],dec[i],r['ra'],r['dec']) < 2*max(RADII[k],r['radius_arcmin']) for r in rows):
            continue
        rows.append(dict(rank=len(rows)+1,ra=ra[i],dec=dec[i],dm=DM[j],radius_arcmin=RADII[k],score=best[yy[i],xx[i]]))
    return rows


def run():
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    nulls, all_peaks, geometry = {}, [], []
    for field in FIELDS:
        grid = FieldGrid(field)
        stars = catalogue(field)
        pix, inside = grid.pixels(stars.ra.to_numpy(), stars.dec.to_numpy())
        stars = stars.loc[inside]
        weight = weights(*(stars[c].to_numpy() for c in ['g','r','gerr','rerr']))
        if not len(stars):
            raise ValueError(f'Empty selection in {field}')
        score = grid.scores(*grid.moments(pix,weight))
        all_peaks.extend(dict(field=field,**p) for p in peaks(grid,score))
        for k,rad in enumerate(RADII):
            # Area applies to detection centres; masks/clipping of individual sightlines enter separately.
            geometry.append(dict(field=field,radius_arcmin=rad,n_stars=len(stars),n_centres=grid.valid[k].sum(),
                                 grid_area_sqdeg=grid.valid[k].sum()*STEP**2))
        cells = (pix//grid.nx//10)*100 + pix%grid.nx//10
        groups = [np.where(cells==key)[0] for key in np.unique(cells)]
        ids = np.arange(len(pix))
        maxima = []
        for _ in range(NNULL):
            for group in groups:
                ids[group] = rng.permutation(group)
            maxima.append(float(grid.scores(*grid.moments(pix,weight[:,ids])).max()))
        nulls[field] = maxima
    nulls['global'] = np.max([nulls[k] for k in FIELDS],axis=0)
    frame = pd.DataFrame(all_peaks)
    frame['global_p'] = [(1+(nulls['global']>=value).sum())/(NNULL+1) for value in frame.score]
    frame['retained'] = frame.global_p <= .05
    controls = []
    for field,(ra,dec,dm) in CONTROLS.items():
        sub = frame[frame.field==field].copy()
        sub['offset_arcmin'] = FieldGrid(field).dist(sub.ra,sub.dec,ra,dec)
        sub['dm_offset'] = sub.dm-dm
        local = sub[sub.offset_arcmin<=np.maximum(6,sub.radius_arcmin)].sort_values('rank')
        best = local.iloc[0] if len(local) else sub.nsmallest(1,'offset_arcmin').iloc[0]
        controls.append(dict(field=field,ref_ra=ra,ref_dec=dec,ref_dm=dm,**{c:best[c] for c in [
            'rank','ra','dec','dm','radius_arcmin','score','global_p','offset_arcmin','dm_offset']},
            spatial_recovery=bool(len(local) and local.retained.any()),
            joint_recovery=bool((local.retained & (local.dm_offset.abs()<=1)).any())))
    frame.to_csv(OUT/'peaks.csv',index=False,float_format='%.8f')
    pd.DataFrame(controls).to_csv(OUT/'controls.csv',index=False,float_format='%.8f')
    pd.DataFrame(nulls).to_csv(OUT/'field_maxima_null.csv',index=False,float_format='%.12g')
    pd.DataFrame(geometry).to_csv(OUT/'geometry.csv',index=False,float_format='%.8f')
    threshold = np.sort(nulls['global'])[-5]
    (OUT/'calibration.json').write_text(json.dumps(dict(seed=SEED,permutations=NNULL,
        n_retained=int(frame.retained.sum()),strict_threshold=float(threshold),
        p_method='(1 + null maxima >= observed) / 100; maximum over both fields and all search cells'),indent=2)+'\n')


if __name__ == '__main__':
    run()
