"""Test each retained track against displaced paths and two stationary backgrounds."""
import argparse
import csv
import json
import warnings
from pathlib import Path

import numpy as np
from scipy import ndimage

import screen
from screen import ROOT, geometry, ephemeris, prepare


def tracks(frame, radii, phases, epoch):
    angle = phases + geometry.frequency('Uranus', radii) * (frame['t'] - epoch) * 86400
    u, v = radii * np.cos(angle), radii * np.sin(angle)
    m = frame['proj']
    scale = geometry.ARCSEC / (frame['distance'] + m[2,0]*u + m[2,1]*v)
    return scale * (m[0,0]*u + m[0,1]*v), scale * (m[1,0]*u + m[1,1]*v)


def sample(frame, x, y):
    # Displaced control paths are always evaluated at the same exposure times.
    valid = ndimage.map_coordinates(frame['good'], [y,x], order=1, mode='constant') > .98
    value = ndimage.map_coordinates(frame['score'], [y,x], order=1, mode='constant')
    return np.where(valid.all(-1), value.mean(-1), np.nan)


def local_significance(values):
    control = values[:,9:]
    med = np.nanmedian(control, axis=1)
    mad = 1.4826 * np.nanmedian(abs(control-med[:,None]), axis=1)
    return (values[:,:9]-med[:,None]) / np.maximum(.12, mad[:,None])


def run(args):
    rows = list(csv.DictReader((ROOT/'crossings.csv').open()))
    radii = np.array([float(r['radius_km']) for r in rows])[:,None]
    phases = np.deg2rad([float(r['phase_deg']) for r in rows])[:,None]
    epoch = json.loads((ROOT/'summary.json').read_text())['epoch_emission_jd_tdb']
    frames = []
    for table, path in ((ROOT.parent/'jwst6379_9482/centers_6379.csv',args.long),
                        (ROOT/'short_centers.csv',args.short)):
        for row in csv.DictReader(table.open()):
            file = args.files / (row['exposure_id']+'.fits')
            frame = screen.make_frame(file,row,path)
            from astropy.io import fits
            header = fits.getheader(file,0)
            mjd = .5*(header['EXPSTART']+header['EXPEND'])
            _,_,_,obs = ephemeris.interpolate(path,np.array([[mjd]]))
            frame['sky'] = obs[0,0,0,:2]
            frame['center'] = np.array([float(row['x']),float(row['y'])])
            frames.append(frame)
    # All filters are on the same detector and roll. A small common translation tolerates
    # coarse orbital sampling. At each translation, compare the same independently
    # placed surrounding paths. The grid is in native pixels.
    shifts = [(a,b) for a in (-2.,0.,2.) for b in (-2.,0.,2.)]
    shifts += [(rad*np.cos(k*np.pi/8),rad*np.sin(k*np.pi/8)) for rad in (6,10,15) for k in range(16)]
    dx,dy = np.array(shifts).T[:,:,None]
    all_values=[]
    for j,f in enumerate(frames):
        e,n = tracks(f,radii,phases,epoch)
        x,y = f['map'](e,n)
        raw = sample(f,x[:,None,:]+dx,y[:,None,:]+dy)
        sidereal=[]; relative=[]
        ra,dec = prepare.inverse_sky(*f['sky'],e/geometry.ARCSEC,n/geometry.ARCSEC)
        for i,g in enumerate(frames):
            if abs(float(np.mean(g['t'])-np.mean(f['t'])))<.032:continue
            # Use independently timed frames with a meaningful displacement on the sky.
            fixed = geometry.sky_xy(ra,dec,*g['sky'],1.,[0.,0.])
            sx,sy = g['map'](-fixed[...,0],fixed[...,1])
            px,py = g['map'](e,n)
            sidereal.append(sample(g,sx[:,None,:]+dx,sy[:,None,:]+dy))
            relative.append(sample(g,px[:,None,:]+dx,py[:,None,:]+dy))
        sky_med,planet_med = [np.nanmedian(stack,axis=0) for stack in (sidereal,relative)]
        all_values.append([raw,raw-sky_med,raw-planet_med])
    all_values = np.array(all_values)
    names=['raw','sidereal','planet_fixed']
    criteria=[]
    for k,name in enumerate(names):
        values = all_values[:,k]
        z=[]
        for which,minimum in ((slice(0,10),4),(slice(0,5),2),(slice(5,10),2),(slice(10,12),2)):
            piece=values[which]; support=np.isfinite(piece).sum(0)
            stack=np.where(support>=minimum,np.nanmean(np.clip(piece,-15,15),axis=0),np.nan)
            z.append(local_significance(stack))
        rank=np.minimum.reduce([z[0]/5,z[1]/2,z[2]/2,z[3]/5])
        criteria.append(rank)
        for row,group in zip(rows,rank):
            row[name+'_best_joint_rank']=round(float(np.nanmax(group)),3) if np.isfinite(group).any() else ''
    joint=np.minimum.reduce(criteria)
    for r,vals in zip(rows,joint):
        r['all_controls_joint_rank']=round(float(np.nanmax(vals)),4) if np.isfinite(vals).any() else ''
    screen.dump('vetting.csv',rows)
    # A missing background comparison or masked path is indeterminate, never a detection.
    summary={name:int(np.any(rank>=1,axis=1).sum()) for name,rank in zip(names,criteria)}
    summary.update(tested=len(rows),all_three_same_translation=int(np.any(joint>=1,axis=1).sum()),
                   all_three_evaluable=int(np.any(np.isfinite(joint),axis=1).sum()),
                   highest_joint=float(np.nanmax(joint)))
    (ROOT/'vetting.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(summary,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('files','long','short'):parser.add_argument(key,type=Path)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        run(parser.parse_args())
