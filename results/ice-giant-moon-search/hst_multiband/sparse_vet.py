"""Freeze sparse-bank survivors and measure them in independent filters and motion frames."""
import argparse
import csv
import gzip
import warnings
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage as ndi
from screen import ROOT, census, reduce
from sparse import rank
from stationary import clean


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def ensure_file(name, cache, mirrors):
    if (cache / name).exists():
        return
    for mirror in mirrors:
        try:
            (cache / name).hardlink_to(mirror / name)
            return
        except FileNotFoundError:
            pass
    census.download({'rootname': name.split('_')[0]}, cache)


def exposure_geometry(frame, path, eph, planet):
    # Rebuild the same WCS/chip/crop as prep, adding physical pixels and a sidereal transform.
    primary = reduce.NAMES[planet][0]
    sky = eph[primary]
    with fits.open(path) as h:
        chip = h['SCI']
        top = left = 0
        if max(chip.shape) > 700:
            candidates = []
            for s in (d for d in h if d.name == 'SCI'):
                w = WCS(s.header, fobj=h)
                x, y = w.all_world2pix([sky], 0, quiet=True)[0]
                ny, nx = s.shape
                if 0 <= x < nx and 0 <= y < ny:
                    candidates.append((min(x, nx-1-x, y, ny-1-y), s, x, y))
            _, chip, x, y = max(candidates, key=lambda i:i[0])
            ny, nx = chip.shape
            left = int(np.clip(round(x)-256, 0, nx-512))
            top = int(np.clip(round(y)-256, 0, ny-512))
        w = WCS(chip.header, fobj=h)
        raw_center = w.all_world2pix([sky], 0, quiet=True)[0]
        # Register the native celestial WCS to the limb/known-moon center used in prep.
        correction = raw_center - frame['c']
        physical0 = np.array([left-chip.header.get('LTV1', 0), top-chip.header.get('LTV2', 0)])
        bodies = np.array(sorted(k for k in eph if k != primary))
        known_xy = w.all_world2pix([eph[k] for k in bodies], 0, quiet=True)-correction
        return dict(wcs=w, correction=correction, physical0=physical0,
                    chip=chip.header['CCDCHIP'], bodies=bodies, known_xy=known_xy)


def load(tag, cache, mirrors, ephemerides):
    planet = tag.split('-')[0].title()
    path = ROOT/'epochs'/tag
    hashes = read(path/'inputs.csv')
    band = tag.split('-')[2].upper()
    sel = [dict(rootname=r['file'].split('_')[0], target=planet.upper(), filter=band,
                visit=r['file'][4:6]) for r in hashes]
    names = {r['rootname'] for r in sel}
    ep = [r for r in ephemerides if r['rootname'] in names]
    with ThreadPoolExecutor(4) as ex:
        list(ex.map(lambda r:ensure_file(r['file'],cache,mirrors), hashes))
    original = reduce.read_csv
    lookup = {'selection.csv':sel,'inputs.csv':hashes,'ephemerides.csv':ep}
    reduce.read_csv = lookup.__getitem__
    try:
        fs, reg, _, _ = reduce.prep(cache, planet, allow_missing_control=True)
    finally:
        reduce.read_csv = original
    by_name = {}
    for row in ep:
        by_name.setdefault(row['rootname'], {})[int(row['body'])] = [float(row['ra']),float(row['dec'])]
    used = [r['rootname'] for r in reg if r['status']=='used']
    for frame, name in zip(fs, used, strict=True):
        frame['rootname'] = name
        frame.update(exposure_geometry(frame,cache/(name+'_flt.fits'),by_name[name],planet))
    return fs


def transformed_bank(params, origin, frames, planet):
    # reduce.orbit parameterizes phase at the *mean emission MJD of its own frames*.
    bank = params.copy()
    dt = (np.mean([f['mid'] for f in frames])-origin)*86400
    bank[:,1] = np.mod(bank[:,1]+reduce.frequency(planet,bank[:,0])*dt,2*np.pi)
    return bank


def sampled(frames, coords):
    vals, valid = [], []
    for f,xy in zip(frames,coords):
        if not np.all(np.isfinite(xy)):
            vals.append(0.); valid.append(False);continue
        rc = np.asarray(xy)[::-1,None]
        vals.append(float(ndi.map_coordinates(f['score'],rc,order=1,mode='constant')[0]))
        valid.append(ndi.map_coordinates(f['support'].astype('f4'),rc,order=1,mode='constant')[0]>.995)
    return np.array(vals), np.array(valid)


def combined(z, ok, visits):
    # Each physical visit receives one normalized score even for three-exposure filters.
    sums, n = [], []
    for v in np.unique(visits):
        use = ok & (visits==v)
        if np.any(use):
            sums.append(np.clip(z[use],-7,7).sum()/np.sqrt(sum(use)))
            n.append(int(sum(use)))
    return sum(sums)/np.sqrt(max(1,len(sums))),len(sums),sum(n),sum(g>2.5 for g in sums)


def morphology(f, xy):
    # Native-electron aperture photometry and the displacement of the smoothed peak in a 7x7 box.
    gy,gx = np.mgrid[-9:10,-9:10]
    stampsys=np.array([gy+xy[1],gx+xy[0]])
    pix=ndi.map_coordinates(f['image'],stampsys,order=1,mode='constant',cval=np.nan)
    good=ndi.map_coordinates(f['good'].astype(float),stampsys,order=0,mode='constant')>.99
    rad=np.hypot(gx,gy);core=(rad<=2.5)&good
    sky=(rad>=5)&(rad<=9)&good&np.isfinite(pix)
    if core.sum()<15 or sky.sum()<50:return (np.nan,)*3
    med=np.median(pix[sky]); noise=1.4826*np.median(abs(pix[sky]-med))
    counts=(pix[core]-med).sum()
    gy7,gx7=np.mgrid[-3:4,-3:4]
    view=ndi.gaussian_filter(np.where(good,pix-med,0),1)[6:13,6:13]
    iy,ix=np.unravel_index(np.nanargmax(view),view.shape)
    return float(counts),float(counts/(max(1,noise)*np.sqrt(core.sum()))),float(np.hypot(gx7[iy,ix],gy7[iy,ix]))


def native_peak(frame, xy):
    # Retain potential cosmic rays in the diagnostic. The screening mask is recorded separately.
    ix,iy=np.rint(xy).astype(int)
    stamp=frame['image'][iy-9:iy+10,ix-9:ix+10].astype(float)
    if stamp.shape!=(19,19) or not np.isfinite(stamp).all():return None
    yy,xx=np.indices(stamp.shape)
    sky=stamp[np.hypot(xx-9,yy-9)>6]
    med=np.median(sky); std=max(1.,1.4826*np.median(abs(sky-med)))
    signal=stamp-med
    py,px=np.unravel_index(np.argmax(signal[5:14,5:14]),(9,9))
    px+=5;py+=5
    near=np.maximum(signal[py-2:py+3,px-2:px+3],0)
    if near.sum()<=0:return None
    wy,wx=np.indices(near.shape)
    center=np.array([ix-9+px-2+np.sum(near*wx)/near.sum(),iy-9+py-2+np.sum(near*wy)/near.sum()])
    return center,signal[py,px]/std,signal[py,px]/near.sum()


def test_anchors(anchors, frames, cleaned, params, common, coords):
    rows=[]
    for anchor in anchors:
        i=anchor['index']
        for j,g in enumerate(frames):
            if g['visit']!=anchor['visit']:continue
            pos_sky=g['wcs'].all_world2pix([anchor['sky']],0,quiet=True)[0]-g['correction']
            prediction={'frozen_orbit':coords[i,j],'fixed_sky':pos_sky}
            if anchor['chip']==g['chip']:prediction['fixed_detector']=anchor['physical']-g['physical0']
            for model,xy in prediction.items():
                z,ok=sampled(cleaned[j:j+1],xy[None,:])
                zraw,_=sampled(frames[j:j+1],xy[None,:])
                counts,ap,_=morphology(g,xy)
                rows.append(dict(source_epoch=common['source_epoch'],template=i+1,a_km=params[i,0],
                    filter=common['filter'],source_root=anchor['name'],source_visit=g['visit'],
                    source_peak_mad=anchor['peak_mad'],source_peak_fraction_5x5=anchor['fraction'],
                    source_peak_offset_pix=anchor['distance'],other_root=g['rootname'],model=model,
                    test_x=xy[0],test_y=xy[1],supports=int(ok[0]),normalized_score=z[0],uncleaned_score=zraw[0],
                    aperture_electrons=counts,aperture_local_contrast=ap))
    return rows


def write(path, rows):
    reduce.write_csv(path,rows)


def vet(tag, a):
    path=ROOT/'epochs'/tag
    planet,proposal,band,mjd=tag.split('-');planet=planet.title()
    with gzip.open(ROOT/(planet.lower()+'_positions.csv.gz'),'rt') as stream:
        ep=list(csv.DictReader(stream))
    survivors=read(path/'sparse_survivors.csv')
    params=np.array([[float(r['a_km']),float(r['phase_rad'])] for r in survivors])
    summary=[]; individual=[]; alternatives=[];anchors=[];anchor_checks=[]
    for filtername in dict.fromkeys([band]+a.filters):
        etag=f'{planet.lower()}-{proposal}-{filtername.lower()}-{mjd}'
        if not (ROOT/'epochs'/etag/'inputs.csv').exists():continue
        fs=load(etag,a.cache,a.mirrors,ep)
        if not fs:continue
        if filtername==band: origin=np.mean([f['mid'] for f in fs])
        bank=transformed_bank(params,origin,fs,planet)
        x,y=reduce.orbit(fs,planet,bank[:,0],bank[:,1])
        coordinates=np.stack([x,y],axis=-1)
        cleaned=clean(fs)
        vs=np.array([f['visit'] for f in fs])
        for i,(par,coords) in enumerate(zip(params,coordinates)):
            zs,ok=sampled(cleaned,coords)
            z0,_=sampled(fs,coords)
            combo,nv,ne,hi=combined(zs,ok,vs)
            # One-degree phase shifts appraise the scale of same-radius noise.
            # Exclude +/-5 degrees around the frozen track. Apply the same visit-coverage floor.
            off=np.arange(6,355,1)
            bg=np.column_stack([np.full(len(off),bank[i,0]),bank[i,1]+np.deg2rad(off)])
            bz,bok=reduce.samples(cleaned,planet,bg)
            null=[c[0] for z,o in zip(bz,bok) if (c:=combined(z,o,vs))[1]>=max(4,nv-1)]
            conf=(1+sum(q>=combo for q in null))/(1+len(null))
            if filtername==band:
                value=rank(zs[None,:],ok[None,:],vs)[0][0]
                if not np.isclose(value,float(survivors[i]['after']),rtol=1e-4):
                    raise ValueError(f'Failed same-filter replay {i}: {value}')
            common=dict(source_epoch=tag,template=i+1,a_km=par[0],phase_at_origin_rad=par[1],origin_emission_mjd=origin,filter=filtername.upper())
            summary.append(dict(**common,test_epoch=etag,registered=len(fs),usable_exposures=ne,usable_visits=nv,visits_gt2p5=hi,
                                normalized_combination=combo,same_radius_tail=conf,background_phases=len(null),
                                uncleaned_combination=combined(z0,ok,vs)[0]))
            for j,(f,xy) in enumerate(zip(fs,coords)):
                separation=np.linalg.norm(f['known_xy']-xy,axis=1)
                nearest=int(np.nanargmin(separation));flux,ap,peak=morphology(f,xy)
                individual.append(dict(**common,rootname=f['rootname'],visit=f['visit'],mid_emission_mjd=f['mid'],
                    x=xy[0],y=xy[1],ccdchip=f['chip'],physical_x=xy[0]+f['physical0'][0],physical_y=xy[1]+f['physical0'][1],
                    supports=int(ok[j]),score_after_static=zs[j],uncleaned_score=z0[j],
                    aperture_electrons=flux,aperture_local_contrast=ap,smoothed_peak_offset_pix=peak,
                    nearest_moon_id=f['bodies'][nearest],nearest_moon_separation_pix=separation[nearest]))
            if filtername!=band or not ok.any():continue
            for j,(f,xy) in enumerate(zip(fs,coords)):
                if not ok[j] or zs[j]<=2.5:continue
                pk=native_peak(f,xy)
                if pk is None:continue
                center,peak_mad,fraction=pk
                anchors.append(dict(index=i,visit=f['visit'],name=f['rootname'],chip=f['chip'],
                    physical=center+f['physical0'],sky=f['wcs'].all_pix2world([center+f['correction']],0)[0],
                    distance=np.linalg.norm(center-xy),peak_mad=peak_mad,fraction=fraction))
            # Anchor at every supported picked exposure; propagating fixed detector / fixed sky
            # and holding its anchor out prevents an excursion from self-confirming.
            for j, (f,xy) in enumerate(zip(fs,coords)):
                if not ok[j]:continue
                physical=xy+f['physical0'];sky=f['wcs'].all_pix2world([xy+f['correction']],0)[0]
                anchored={'orbit':coords,'fixed_detector':[physical-g['physical0'] if g['chip']==f['chip'] else [np.nan]*2 for g in fs],
                          'fixed_sky':[g['wcs'].all_world2pix([sky],0,quiet=True)[0]-g['correction'] for g in fs]}
                for model,locations in anchored.items():
                    zz,support=sampled(cleaned,locations);support &= vs!=vs[j]
                    value,nvis,nexp,nh=combined(zz,support,vs)
                    alternatives.append(dict(**common,anchor=f['rootname'],model=model,other_visits=nvis,
                                             other_exposures=nexp,normalized_other_visits=value,other_visits_gt2p5=nh))
        if filtername!=band:
            anchor_checks.extend(test_anchors(anchors,fs,cleaned,params,common,coordinates))
            if anchor_checks:write(path/'sparse_vet_anchors.csv',anchor_checks)
        print(tag,filtername,'complete',flush=True)
        write(path/'sparse_vet_filters.csv',summary);write(path/'sparse_vet_exposures.csv',individual)
        if alternatives:write(path/'sparse_vet_motion.csv',alternatives)
        del cleaned,fs


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('source_epoch',nargs='+');p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--mirrors',type=Path,nargs='*',default=[]);p.add_argument('--filters',nargs='*',default=['f845m','f763m'])
    args=p.parse_args();args.cache.mkdir(parents=True,exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',category=RuntimeWarning)
        for epoch in args.source_epoch:vet(epoch,args)
