"""Evaluate paired exposure-stack peaks on native time-resolved detector integrations."""
import argparse
import csv
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory
import warnings

import numpy as np
base=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(base));sys.path.insert(0,str(base/'jwst9482'))
from link_visits import geometry,search,trials
from prepare import inverse_sky
import positive_control as native


def summaries(images):
    arr=np.array(images)
    mid=len(arr)//2
    with warnings.catch_warnings():
        warnings.simplefilter('ignore',RuntimeWarning)
        stacks=[np.nanmedian(z,axis=0) for z in [arr,arr[:mid],arr[mid:]]]
    yy,xx=np.indices(stacks[0].shape)-native.RADIUS
    ann=(xx*xx+yy*yy>=8**2)&(xx*xx+yy*yy<=20**2)
    small=xx*xx+yy*yy<=3**2
    def standardized(im):
        bg=np.nanmedian(im[ann])
        noise=1.4826*np.nanmedian(abs(im[ann]-bg))
        return (im-bg)/noise
    ims=[standardized(im) for im in stacks]
    im=ims[0]
    if not np.isfinite(im[small]).any():
        return [0,0,0,0,0,0,0,0,0],im
    y,x=np.unravel_index(np.nanargmax(np.where(small,im,np.nan)),im.shape)
    cy,cx=[native.RADIUS]*2
    n=int(np.isfinite(arr[:,cy,cx]).sum())
    return [n,im[cy,cx],im[y,x],x-cx,y-cy,ims[1][y,x],ims[2][y,x],
            ims[1][cy,cx],ims[2][cy,cx]],im


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for n in ['early','later','links','meta1','pixels1','meta2','pixels2','out']:
        p.add_argument(n,type=Path)
    a=p.parse_args();paths=[a.early,a.later]
    params=json.loads((a.links/'parameters.json').read_text())
    representatives={}
    with open(a.links/'passing_pairs.csv') as stream:
        for row in csv.DictReader(stream):
            if (row['early_sign'],row['late_sign']) != ('1','1'):continue
            group=(int(row['early_peak']),int(row['late_peak']))
            score=min(float(row['z_first']),float(row['z_second']))
            if group not in representatives or score>representatives[group][0]:
                representatives[group]=(score,row)
    groups=list(representatives)
    elem=np.zeros((len(groups),6))
    for i,key in enumerate(groups):
        row=representatives[key][1]
        elem[i,[0,5]]=[float(row['a_km']),float(row['phase_first_deg'])]
    a.out.mkdir(parents=True,exist_ok=False)
    results={}; receipts=set()
    for j,(path,meta_cache,pixels) in enumerate(zip(paths,[a.meta1,a.meta2],[a.pixels1,a.pixels2],strict=True)):
        cube=search.load_npz(path)
        with TemporaryDirectory() as td:
            tmp=Path(td)
            bank=trials(tmp/'e.csv',path,tmp/'b.npz',elem,params['epoch1_tdb_jd'])
        records=json.loads(Path(str(path)+'.inputs.json').read_text())
        sizes=len(groups)
        all_images=[[] for _ in groups]
        min_known=[[] for _ in groups]; nearest=[[] for _ in groups]
        for f,record in enumerate(records):
            m=native.metadata(meta_cache,record['exposure_id'].removesuffix('_cal'))
            data,x0,y0=native.wide_images(pixels,m)
            # Crop origins remain those of the original header-based byte range.
            # Trail predictions use the independently measured limb center.
            ra,dec=m['wcs'].all_pix2world(record['center_native_x'],record['center_native_y'],0)
            ts=cube['path_mjd'][f]
            pp={}
            mid=m['times']
            for i,key in enumerate(groups):
                x=np.interp(mid,ts,bank['x_path'][i,f])
                y=np.interp(mid,ts,bank['y_path'][i,f])
                scale=cube['scale_arcsec']/geometry.ARCSEC
                ra1,dec1=inverse_sky(ra,dec,-(x-cube['center_xy'][0])*scale,(y-cube['center_xy'][1])*scale)
                xx,yy=m['wcs'].all_world2pix(ra1,dec1,0,quiet=True)
                pp[i]=np.stack([xx,yy],axis=-1)
                body=cube['known_path'][:,f]
                # Same-time separation here covers the 204 s integration rather than the 34 min exposure.
                tr=np.stack([bank['x_path'][i,f],bank['y_path'][i,f]],axis=-1)
                # For each integration include any pair of native path epochs inside the exposure interval.
                mmin=np.full(len(body),np.inf)
                for t in mid:
                    inside=abs(ts-t) <= m['sec']/2/86400+1e-6
                    if inside.any():
                        sep=np.linalg.norm(body[:,inside,:,None]-tr[inside].T[None,None,:,:],axis=2)
                        mmin=np.minimum(mmin,np.min(sep,axis=(1,2)))
                idx=np.argmin(mmin);min_known[i].append(float(mmin[idx]));nearest[i].append(int(cube['known_id'][idx]))
            patches=native.patches_for_positions(data,x0,y0,pp)
            for i in range(sizes):all_images[i].extend(patches[i])
        for i,key in enumerate(groups):
            stats,im=summaries(all_images[i])
            k=np.argmin(min_known[i])
            results[(j,key)]=(stats,nearest[i][k],min_known[i][k])
        receipts.update(native.USED)
    keys=['valid_integrations','center_contrast','peak_contrast','dx','dy','half1_at_peak','half2_at_peak','half1_center','half2_center','nearest_known_id','nearest_integration_trail_px']
    with open(a.out/'integration_vetting.csv','w',newline='') as stream:
        w=csv.writer(stream)
        w.writerow(['early_peak','late_peak','a_km','phase_first_deg','coarse_first_z','coarse_second_z',
                    *['early_'+x for x in keys],*['late_'+x for x in keys]])
        for i,key in enumerate(groups):
            row=representatives[key][1]
            vals=[]
            for j in range(2):
                z,b,d=results[(j,key)];vals += [*z,b,d]
            w.writerow([*key,elem[i,0],elem[i,-1],row['z_first'],row['z_second'],*vals])
    with open(a.out/'range_hashes.csv','w',newline='') as stream:
        w=csv.writer(stream);w.writerow(['range_file','sha256'])
        for f in sorted(receipts):w.writerow([f.name,search.file_hash(f)])


if __name__=='__main__':
    main()
