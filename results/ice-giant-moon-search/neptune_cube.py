"""Register ice-giant detector tiles and coadd chips in each physical exposure."""
import argparse
import csv
from collections import defaultdict
import json
from pathlib import Path
import tempfile
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from geometry import ARCSEC, sky_xy
from limb_center import measure
from prepare import inverse_sky, prepare


def sequence(path):
    return path.stem.rsplit('_',2)[0]


def run(files,eph,band,out,receipt,control,size=401):
    paths=sorted(p for p in map(Path,files) if fits.getval(p,'FILTER')==band)
    groups=defaultdict(list)
    for p in paths:
        groups[sequence(p)].append(p)
    with np.load(eph) as f:
        radius=f['values'][0,0,2]/.062
    offsets={}; diagnostics=[]
    for p in paths:
        fit=measure(p,radius)
        if fit is None:
            continue
        hdr=fits.getheader(p,'SCI')
        w=WCS(hdr,naxis=2)
        ra,dec=w.all_pix2world(fit['x'],fit['y'],0)
        offsets[sequence(p)]=sky_xy(ra,dec,hdr['MT_RA'],hdr['MT_DEC'],1.,(0.,0.))
        diagnostics.append({**fit,'east_shift_arcsec':-float(offsets[sequence(p)][0]),
                            'north_shift_arcsec':float(offsets[sequence(p)][1])})
    if len(offsets)<2:
        raise ValueError('insufficient enclosed planetary limbs')
    visit_medians={}
    for seq in groups:
        visit=seq.split('_')[0]
        offsets_in_visit=[v for group,v in offsets.items() if group.split('_')[0]==visit]
        if not offsets_in_visit:
            raise ValueError('missing enclosed limb in visit '+visit)
        visit_medians[visit]=np.median(offsets_in_visit,axis=0)
    records=[]; registration=[]
    with tempfile.TemporaryDirectory(dir=Path(out).parent) as tmp:
        centers=Path(tmp)/'centers.csv'
        for group,items in sorted(groups.items()):
            # Apply one moving-target offset from a contemporaneous fitted limb,
            # transferring the same-band median if the limb hits a detector edge.
            dx,dy=offsets.get(group,visit_medians[group.split('_')[0]])
            for p in items:
                h=fits.getheader(p,'SCI')
                ra,dec=inverse_sky(h['MT_RA'],h['MT_DEC'],-dx/ARCSEC,dy/ARCSEC)
                x,y=WCS(h,naxis=2).all_world2pix(ra,dec,0,quiet=True)
                records.append(dict(exposure_id=p.stem,x=float(x),y=float(y),ra_deg=float(ra),dec_deg=float(dec)))
            registration.append(dict(exposure_id=group,source='limb' if group in offsets else 'median of same-band limb fits',
                                     east_shift_arcsec=-float(dx),north_shift_arcsec=float(dy),n_chip_tiles=len(items)))
        with centers.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=records[0]);w.writeheader();w.writerows(records)
        temp=Path(tmp)/'tiles.npz'
        prepare(paths,eph,temp,centers,size,.031,17,{'F140M':.84,'F182M':.96,'F210M':1.06}[band],5.,True,control)
        with np.load(temp) as f:
            cube={k:f[k] for k in f.files}
        names=cube['exposure_id']
        # Keep the first calibrated value if overlapping tile boundary mappings
        # place two valid chips on one output pixel; record any overlaps.
        counts=[]; im=[];var=[];good=[];geom=[]; idx=[]
        for group in sorted(groups):
            selected=np.array([i for i,k in enumerate(names) if k.rsplit('_',2)[0]==group])
            idx.append(selected[0]); mask=cube['valid'][selected]
            weights=mask.sum(0); g=weights>0; counts.append(int((weights>1).sum()))
            z=mask.argmax(0)[None]
            im.append(np.take_along_axis(cube['image'][selected],z,axis=0)[0])
            var.append(np.take_along_axis(cube['variance'][selected],z,axis=0)[0])
            good.append(g);geom.append(cube['geometry_valid'][selected].any(0))
        for k in ['mjd','visit','known_xy','ring_mask','path_ra','path_dec','path_delta','path_mjd','r_au','exptime_s']:
            cube[k]=cube[k][idx]
        cube['known_path']=cube['known_path'][:,idx]
        cube.update(exposure_id=np.array(sorted(groups)),image=im,variance=var,valid=good,geometry_valid=geom)
        np.savez_compressed(out,**cube)
        receipt_data=dict(band=band,unmasked_control_body=control,limbs=diagnostics,
                          registrations=registration,duplicate_valid_pixels_by_exposure=counts,
                          tile_inputs=json.loads(Path(str(temp)+'.inputs.json').read_text()))
        Path(receipt).write_text(json.dumps(receipt_data,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('ephemeris');p.add_argument('band');p.add_argument('output');p.add_argument('receipt')
    p.add_argument('files',nargs='+');p.add_argument('--control',type=int)
    p.add_argument('--size',type=int,default=401)
    a=p.parse_args();run(a.files,a.ephemeris,a.band,a.output,a.receipt,a.control,a.size)
