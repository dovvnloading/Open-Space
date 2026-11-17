"""Moving-track and source-injection tests on native SPHEREx DR1 and QR3 cutouts."""
import argparse
import csv
import gzip
import io
import hashlib
import importlib.util
import json
import warnings
from pathlib import Path
from urllib.request import urlopen

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from scipy.special import ndtri
from scipy.stats import qmc
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.constants import c, h, k as boltzmann, astronomical_unit

ROOT = Path(__file__).resolve().parent
N_CONTROL = 1500
N_INJECT = 320
SEED = 20261117
FLUX_LEVELS_AB = np.arange(16., 22.51, .5)


def filtered(pixels, mask):
    """Local sky subtraction and smoothing, with flagged samples omitted."""
    valid = mask.astype(float)
    pixel_values = np.where(mask,pixels,0.)
    small_weight = gaussian_filter(valid, .8, mode="constant")
    wide_weight = gaussian_filter(valid, 3., mode="constant")
    small = gaussian_filter(pixel_values,.8,mode="constant") / np.maximum(small_weight,1e-8)
    wide = gaussian_filter(pixel_values,3.,mode="constant") / np.maximum(wide_weight,1e-8)
    signal = small-wide
    signal[(small_weight < .7) | (wide_weight < .7)] = np.nan
    return signal


def sample(image,x,y):
    return map_coordinates(image,[np.asarray(y),np.asarray(x)],order=1,
                           mode="constant",cval=np.nan,prefilter=False)


def astro_xy(wcs, ra, dec):
    return np.asarray(wcs.all_world2pix(ra,dec,0,quiet=True))


class Frame:
    def __init__(self, row, cache):
        path = cache / row["name"]
        if not path.exists():
            with urlopen(row["url"],timeout=120) as response:
                path.write_bytes(response.read())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=row["sha256"]:
            raise ValueError(f"Checksum mismatch for {path}")
        self.name=row["name"]
        self.t=(row["t_min"]+row["t_max"])/2
        self.group = 0 if self.t<60950 else 1 if self.t<61100 else 2
        with fits.open(path) as hdu:
            ih=hdu["IMAGE"].header
            self.wcs=WCS(ih,naxis=2)
            self.aux=np.array([ih["CRPIX1A"]-1,ih["CRPIX2A"]-1])
            # MJy/sr to microJy per pixel at the TAN projection origin.
            factor=1e12*proj_plane_pixel_area(self.wcs)*(np.pi/180)**2
            self.pixels=(hdu["IMAGE"].data-hdu["ZODI"].data)*factor
            self.good=(hdu["FLAGS"].data==0) & np.isfinite(self.pixels) & (hdu["VARIANCE"].data>0)
            self.response=filtered(self.pixels,self.good).astype(np.float32)
            self.epsfs=hdu["EPSF"].data.copy()
            self.sampling=int(hdu["EPSF"].header["OVSMPX"])
            assert self.sampling==hdu["EPSF"].header["OVSMPY"]==5

    def zone(self,x,y):
        # Auxiliary A coordinates are zero-based original detector pixels.
        xd,yd=np.array([x,y])-self.aux
        return int(np.argmin((xd-self.epsfs["XCENTER"])**2+(yd-self.epsfs["YCENTER"])**2))

    def injected(self,x,y,read_xy):
        """Apply a 100-microJy point source to actual science pixels and remeasure."""
        n=22
        ix,iy=round(x),round(y)
        out=np.zeros(len(read_xy))
        ind=np.where((abs(read_xy[:,0]-ix)<n-1)&(abs(read_xy[:,1]-iy)<n-1))[0]
        if not len(ind):return out
        ny,nx=self.pixels.shape
        xmin,xmax=max(0,ix-n),min(nx,ix+n+1)
        ymin,ymax=max(0,iy-n),min(ny,iy+n+1)
        if xmax<=xmin or ymax<=ymin:return out
        cut=self.pixels[ymin:ymax,xmin:xmax]
        mask=self.good[ymin:ymax,xmin:xmax]
        psf=self.epsfs["EPSF"][self.zone(x,y)]
        sy,sx=(np.array(psf.shape)-1)/2
        xx,yy=np.meshgrid(np.arange(xmin,xmax),np.arange(ymin,ymax))
        model=map_coordinates(psf,[(yy-y)*self.sampling+sy,(xx-x)*self.sampling+sx],
                              order=1,mode="constant",cval=0.,prefilter=False)
        # Normalize on the full PSF stamp; at detector edges the coverage cut
        # discards incomplete apertures. Do not scale up if the source is off-frame.
        if not (3<x<nx-4 and 3<y<ny-4) or model.sum()<=0:return out
        model*=100/model.sum()
        difference=filtered(cut+model,mask)-filtered(cut,mask)
        change=sample(difference,read_xy[ind,0]-xmin,read_xy[ind,1]-ymin)
        out[ind]=np.nan_to_num(change,nan=0.)
        return out


def summarize(raw,groups):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",category=RuntimeWarning)
        delta=raw[:,0]-np.nanmedian(raw[:,1:],axis=1)
        delta[np.sum(np.isfinite(raw[:,1:]),axis=1)<2]=np.nan
        count=np.array([np.isfinite(delta[groups==g]).sum(axis=0) for g in range(3)])
        combined=np.array([np.nanmedian(delta[groups==g],axis=0) for g in range(3)])
    combined[count<3]=np.nan
    return combined,count


def write_rows(path,rows):
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,rows[0].keys())
        w.writeheader()
        w.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-dir",type=Path,required=True)
    parser.add_argument("--cache",type=Path,required=True)
    parser.add_argument("--out",type=Path,default=ROOT)
    parser.add_argument("--manifest",type=Path,default=ROOT/"images.json")
    parser.add_argument("--family",choices=("nominal","centroid-stress"),default="nominal")
    args=parser.parse_args()
    args.cache.mkdir(exist_ok=True,parents=True)
    args.out.mkdir(exist_ok=True,parents=True)
    spec=importlib.util.spec_from_file_location("orbit",args.orbit_dir/'propagate.py')
    orbit=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(orbit)
    inputs=json.loads((args.orbit_dir/'orbit_input.json').read_text())
    if args.family=="centroid-stress":
        params=qmc.Sobol(7,scramble=True,seed=19832006).random_base2(15)
        major,minor=ndtri(params[:,3:5]).T
        pa=np.deg2rad(79.)
        iras=np.stack([18*major*np.sin(pa)+7*minor*np.cos(pa),
                       18*major*np.cos(pa)-7*minor*np.sin(pa)],axis=1)
        radius=40*np.sqrt(params[:,5])
        angle=2*np.pi*params[:,6]
        akari=np.stack([radius*np.cos(angle),radius*np.sin(angle)],axis=1)
        state=orbit.family(inputs,parameters=params[:,:3],angular_offsets=[iras,akari])
    else:
        state=orbit.family(inputs)
    frames=[Frame(r,args.cache) for r in json.loads(args.manifest.read_text()) if r.get("sha256")]
    fm=len(frames)
    coords=np.array([orbit.sky(state,Time(f.t,format="mjd",scale="utc").tt.jd-2451545) for f in frames])
    ngrid=coords.shape[-1]
    rng=np.random.default_rng(SEED)
    subset=rng.choice(ngrid,N_CONTROL)
    ph=rng.uniform(0,2*np.pi,N_CONTROL)
    rad=np.sqrt(rng.uniform(8**2,13**2,N_CONTROL))/60
    displacement=np.array([rad*np.cos(ph)/np.cos(np.deg2rad(-49.7)),rad*np.sin(ph)])
    coords=np.concatenate([coords,coords[:,:,subset]+displacement[None,:,:]],axis=-1)
    inject=rng.choice(ngrid,min(ngrid,N_INJECT),replace=False)
    ntrail=coords.shape[-1]
    groups=np.array([f.group for f in frames])
    # Choose four same-detector reference exposures from each other season,
    # using wavelengths at the median trajectory location of the science exposure.
    refs=[]
    for i,fr in enumerate(frames):
        ra,dec=np.median(coords[i,:,:ngrid],axis=-1)
        info=[]
        for j,rf in enumerate(frames):
            x,y=astro_xy(rf.wcs,ra,dec)
            inside=(-2<x<rf.pixels.shape[1]+2 and -2<y<rf.pixels.shape[0]+2)
            wave=rf.epsfs["CWAVE"][rf.zone(x,y)]
            info.append((not inside,float(wave),j))
        row=[i]
        for g in range(3):
            if g==fr.group:continue
            ranks=sorted((outside,abs(w-info[i][1]),j) for outside,w,j in info if groups[j]==g)
            row.extend(r[2] for r in ranks[:4])
        refs.append(row)
    refs=np.array(refs)
    raw=np.full((fm,refs.shape[1],ntrail),np.nan,dtype='f4')
    response=np.zeros((fm,refs.shape[1],len(inject)),dtype='f4')
    reflected = np.zeros((fm, len(inject)))
    pos, vel, t0, _ = state
    for j,fr in enumerate(frames):
        slots=list(zip(*np.where(refs==j)))
        own=astro_xy(fr.wcs,*coords[j][:,inject])
        # Flux at the spacecraft for a full-phase sphere of radius one Earth
        # radius and unit geometric albedo. Evaluate incident sunlight at
        # the nearest published ePSF wavelength on the native detector.
        t=Time(fr.t,format="mjd",scale="utc").tt.jd-2451545
        sun_to_obs, bary_to_obs=orbit.observers(t)
        start=pos[inject]
        travel=orbit.norm(start-bary_to_obs)/orbit.C
        planet=orbit.kepler(start,vel[inject],t-t0[inject]-travel)
        heliocentric=orbit.norm(planet-bary_to_obs+sun_to_obs)
        geocentric=orbit.norm(planet-bary_to_obs)
        wave=np.array([fr.epsfs["CWAVE"][fr.zone(x,y)] for x,y in own.T],dtype=float)*1e-6
        frequency=c/wave
        solar_flux_jy=1e26*np.pi*(695700e3/astronomical_unit)**2*(2*h*frequency**3/c**2)/np.expm1(h*frequency/(boltzmann*5772.))
        reflected[j]=1e6*solar_flux_jy*(6371e3/astronomical_unit/(heliocentric*geocentric))**2
        read_positions=[]
        for i,k in slots:
            xy=astro_xy(fr.wcs,*coords[i])
            raw[i,k]=sample(fr.response,*xy)
            read_positions.append(xy[:,inject].T)
        reads=np.stack(read_positions)
        for t,xy in enumerate(own.T):
            increments=fr.injected(*xy,reads[:,t,:])
            for a,(i,k) in zip(increments,slots):response[i,k,t]=a
    stats,count=summarize(raw,groups)
    valid=np.all(np.isfinite(stats),axis=0)
    ctrl=stats[:,ngrid:][:,valid[ngrid:]]
    center=np.median(ctrl,axis=1)
    mad=1.4826*np.median(np.abs(ctrl-center[:,None]),axis=1)
    def score(s):
        z=(s-center[:,None])/mad[:,None]
        return np.min(z,axis=0),np.mean(z,axis=0)
    def passes(s):
        low,mean=score(s)
        return (low>2.5)&(mean>4.)
    flags=passes(stats)
    low,mean=score(stats)
    det_rows=[dict(grid_index=int(i),minimum_season_z=float(low[i]),mean_season_z=float(mean[i]),
                   ra_first=float(coords[0,0,i]),dec_first=float(coords[0,1,i])) for i in np.where(flags[:ngrid])[0]]
    if det_rows:
        write_rows(args.out/'triggers.csv',det_rows)
    eligible_idx=np.flatnonzero(valid[:ngrid])
    (args.out/'eligible_grid_indices.json').write_text(json.dumps(eligible_idx.tolist())+'\n')
    best=int(np.nanargmax(np.where(valid[:ngrid],low[:ngrid],np.nan)))
    metrics={"family":args.family,"grid_size":ngrid,"grid_valid":int(valid[:ngrid].sum()),"valid_controls":int(ctrl.shape[1]),
            "controls":N_CONTROL,"control_passes":int(flags[ngrid:].sum()),"grid_passes":int(flags[:ngrid].sum()),
            "exposures":fm,"season_exposures":[int((groups==g).sum()) for g in range(3)],
            "season_mjd_ranges":[[float(min(f.t for f in frames if f.group==g)),float(max(f.t for f in frames if f.group==g))] for g in range(3)],
            "control_centers":center.tolist(),"control_robust_scales":mad.tolist(),
            "best_min_season_track":best,"best_min_season_z":float(low[best]),"best_same_track_mean_z":float(mean[best]),
            "seed":SEED,"injected_tracks":len(inject),"eligible_injected_tracks":int(valid[inject].sum()),
            "per_season_track_min_median_max_usable_exposures":[[float(x) for x in np.percentile(c[:ngrid],[0,50,100])] for c in count]}
    reports={"grid_index":inject, "eligible":valid[inject].astype(int)}
    control_indices=np.flatnonzero(valid[ngrid:])
    (args.out/'eligible_control_indices.json').write_text(json.dumps(control_indices.tolist())+'\n')
    (args.out/'triggered_control_indices.json').write_text(json.dumps(np.flatnonzero(flags[ngrid:]).tolist())+'\n')
    rows=[]
    for mag in FLUX_LEVELS_AB:
        flux=3631e6*10**(-.4*mag)
        injected_raw=raw[:,:,inject]+response*(flux/100)
        ss,_=summarize(injected_raw,groups)
        wins=passes(ss)&valid[inject]
        reports[f"AB{mag:.1f}"]=wins.astype(int)
        rows.append(dict(magnitude_AB=mag,flux_microJy=round(flux,5),tested=len(inject),
                         eligible=int(valid[inject].sum()),recovered=int(wins.sum())))
    write_rows(args.out/'injections.csv',rows)
    model_rows=[]
    with (args.orbit_dir/"brightness.csv").open() as f:
        scenarios={(float(r["mass_Mearth"]),float(r["radius_Rearth"]),float(r["teff_K"]),float(r["geometric_albedo_r"])) for r in csv.DictReader(f)}
    for mass, radius, teff, albedo in sorted(scenarios):
        flux=reflected*radius**2*albedo
        ss,_=summarize(raw[:,:,inject]+response*flux[refs]/100,groups)
        wins=passes(ss)&valid[inject]
        mags=-2.5*np.log10(np.median(flux,axis=0)/3631e6)
        reports[f"M{mass:g}_R{radius:g}_T{teff:g}_A{albedo:g}"]=wins.astype(int)
        model_rows.append(dict(mass_Mearth=mass,radius_Rearth=radius,teff_K=teff,
            geometric_albedo_grey=albedo,trackmedian_AB_min=float(min(mags)),trackmedian_AB_max=float(max(mags)),
            injected=len(inject),eligible=int(valid[inject].sum()),recovered=int(wins.sum())))
    write_rows(args.out/'grey_model_injections.csv',model_rows)
    stream=io.StringIO(newline="")
    writer=csv.writer(stream)
    writer.writerow(reports)
    writer.writerows(zip(*reports.values()))
    (args.out/"injection_trials.csv.gz").write_bytes(gzip.compress(stream.getvalue().encode(),mtime=0))
    (args.out/'summary.json').write_text(json.dumps(metrics,indent=2)+'\n')
    print(json.dumps(metrics))
if __name__=="__main__":main()
