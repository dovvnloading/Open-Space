"""Catalogue-level Plummer injections with a Kroupa IMF and the frozen detector."""
import json
from concurrent.futures import ProcessPoolExecutor
from itertools import product
import os
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from search import (ROOT, OUT, FIELDS, CONTROLS, DM, STEP, SEED, FieldGrid,
                    catalogue, weights, peaks)

MAGNITUDES = [-2., -4., -6.]
SIZES_PC = [5., 30., 100.]
MODULI = [17., 19., 21., 23.]
PER_QUADRANT = 6
BAD_MASK = sum(1 << bit for bit in [1,2,3,5,6,10,11,12])


def luminosity_grid():
    a = np.genfromtxt(ROOT/'model/parsec.csv', names=True, delimiter=',')
    a = a[a['stage'] != 9]  # omit terminal, extinguished isochrone placeholder
    edges = np.unique(np.concatenate([np.linspace(x,y,17) for x,y in zip(a['initial_mass'][:-1],a['initial_mass'][1:])]))
    m = (edges[:-1]+edges[1:])/2
    # Kroupa dN/dm, continuous at 0.5 solar masses; quadrature over the track.
    counts = np.diff(edges)*np.where(m<.5,m**-1.3,.5*m**-2.3)
    g = np.interp(m,a['initial_mass'],a['g']); r = np.interp(m,a['initial_mass'],a['r'])
    v = g-.5784*(g-r)-.0038  # Lupton 2005 SDSS relation applied approximately to DECam AB.
    counts /= np.sum(counts*10**(-.4*v))
    return g,r,counts


class Campaign:
    def __init__(self, field):
        self.field = field
        self.grid = FieldGrid(field)
        self.threshold = json.loads((OUT/'calibration.json').read_text())['strict_threshold']
        raw = pd.read_csv(ROOT/f'catalog/{field}.csv.gz',usecols=['ra','dec','gdepth','rdepth','maskbits'])
        self.depth_tree = cKDTree(np.column_stack([raw.ra*self.grid.cos,raw.dec]))
        self.depth = raw[['gdepth','rdepth']].to_numpy()
        self.mask = raw.maskbits.to_numpy().astype(np.int64)
        stars = catalogue(field)
        pix, inside = self.grid.pixels(stars.ra.to_numpy(),stars.dec.to_numpy())
        w = weights(*(stars[c].to_numpy()[inside] for c in ['g','r','gerr','rerr']))
        self.base = self.grid.moments(pix,w)
        self.base_peaks = peaks(self.grid,self.grid.scores(*self.base),self.threshold)
        self.model = luminosity_grid()
        # Uniformly sample continuous positions inside valid 3 arcmin-centre cells.
        ra,dec = self.grid.coords(self.grid.y,self.grid.x)
        cr,cd,_ = CONTROLS[field]
        valid = self.grid.valid[0] & (self.grid.dist(ra,dec,cr,cd)>25.5)
        self.cells = {}
        for cy,cx in product(range(2),repeat=2):
            self.cells[(cy,cx)] = np.where(valid & ((self.grid.y>=self.grid.ny/2)==cy)
                                                & ((self.grid.x>=self.grid.nx/2)==cx))

    def observe(self, ra,dec,g,r,rng):
        sep, idx = self.depth_tree.query(np.column_stack([ra*self.grid.cos,dec]))
        d = self.depth[idx]
        eligible = (sep<=STEP) & ((self.mask[idx]&BAD_MASK)==0) & np.isfinite(d).all(axis=1)
        mags = []; errs = []
        for j,source_mag in enumerate([g,r]):
            f = 10**((22.5-source_mag)/2.5)
            # Add a 0.01-mag multiplicative noise floor to local 5-sigma PSF noise.
            sigma = np.hypot(10**((22.5-d[:,j])/2.5)/5,f*.01*np.log(10)/2.5)
            flux = f+rng.normal(size=len(f))*sigma
            with np.errstate(divide='ignore',invalid='ignore'):
                mags.append(22.5-2.5*np.log10(flux)); errs.append(2.5/np.log(10)*sigma/flux)
            eligible &= (flux>0) & (errs[-1]<.2)
        eligible &= (mags[0]>=17)&(mags[0]<=23.5)&(mags[1]>=16)&(mags[1]<=23.5)
        pix, inside = self.grid.pixels(ra[eligible],dec[eligible])
        w = weights(*(v[eligible][inside] for v in [*mags,*errs]))
        return pix,w

    def one(self,mv,rh,dm,cy,cx,rng):
        y,x = self.cells[(cy,cx)]
        i = rng.integers(len(y))
        ra0,dec0 = self.grid.coords(y[i]+rng.uniform(-.5,.5),x[i]+rng.uniform(-.5,.5))
        g,r,n = self.model
        # Restrict the draw to stars within 2 mag of the faint selection limit;
        # total expected light includes the full stellar mass interval.
        eligible = (g+dm <= 25.5)&(r+dm <= 25.5)
        counts = n[eligible]*10**(-.4*mv)
        chosen = rng.choice(np.where(eligible)[0],size=rng.poisson(counts.sum()),p=counts/counts.sum())
        dist_pc = 10**((dm+5)/5)
        u = rng.random(len(chosen)); angle = rng.uniform(0,2*np.pi,len(chosen))
        radius_deg = np.rad2deg(rh/dist_pc)*np.sqrt(u/(1-u))
        ra = ra0+radius_deg*np.cos(angle)/self.grid.cos; dec = dec0+radius_deg*np.sin(angle)
        pix,w = self.observe(ra,dec,g[chosen]+dm,r[chosen]+dm,rng)
        if len(pix):
            diff,var = self.grid.moments(pix,w)
            result = peaks(self.grid,self.grid.scores(self.base[0]+diff,self.base[1]+var),self.threshold)
        else:
            result = self.base_peaks
        def match(row):
            return abs(row['dm']-dm)<=1 and self.grid.dist(row['ra'],row['dec'],ra0,dec0)<=max(6,row['radius_arcmin'])
        existing = any(match(p) for p in self.base_peaks)
        recovered = next((p for p in result if match(p)),None)
        return dict(field=self.field,mv_expected=mv,rhalf_pc=rh,dm=dm,quadrant_y=cy,quadrant_x=cx,
                    ra=ra0,dec=dec0,n_selected=len(pix),already_above=existing,
                    recovered=bool(recovered and not existing),recovered_score=recovered['score'] if recovered else np.nan)


def run_field(field):
    campaign = Campaign(field)
    rng = np.random.default_rng(SEED+100+list(FIELDS).index(field))
    rows = []
    for mv,rh,dm,cy,cx in product(MAGNITUDES,SIZES_PC,MODULI,range(2),range(2)):
        for _ in range(PER_QUADRANT):
            rows.append(campaign.one(mv,rh,dm,cy,cx,rng))
    return rows


def summarize(df, keys):
    out = df.groupby(keys,as_index=False).agg(trials=('recovered','size'),successes=('recovered','sum'),
                                               mean_stars=('n_selected','mean'),preexisting=('already_above','sum'))
    p = out.successes/out.trials; z = 1.959963984540054
    denom = 1+z**2/out.trials
    center = (p+z**2/(2*out.trials))/denom
    err = z*np.sqrt(p*(1-p)/out.trials+z**2/(4*out.trials**2))/denom
    out['fraction']=p; out['wilson_lo']=np.maximum(0,center-err); out['wilson_hi']=np.minimum(1,center+err)
    return out


def run():
    with ProcessPoolExecutor(max_workers=min(2,os.cpu_count() or 1)) as exe:
        rows = [row for field_rows in exe.map(run_field,FIELDS) for row in field_rows]
    trials = pd.DataFrame(rows)
    trials.to_csv(OUT/'injection_trials.csv.gz',index=False,float_format='%.7f',compression={'method':'gzip','mtime':0})
    keys = ['field','mv_expected','rhalf_pc','dm']
    summarize(trials,keys).to_csv(OUT/'completeness.csv',index=False,float_format='%.6f')
    summarize(trials,keys+['quadrant_y','quadrant_x']).to_csv(OUT/'spatial_completeness.csv',index=False,float_format='%.6f')


if __name__ == '__main__':
    run()
