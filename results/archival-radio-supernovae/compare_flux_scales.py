"""Compare independent compact peaks around each SN in pairs of VLASS images."""
import argparse
import hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from scipy.ndimage import maximum_filter
from scipy.spatial import cKDTree
from screen import download

def peaks(path, row):
    with fits.open(path) as hdus:
        im=hdus[0].data.squeeze().astype(float)*1000
        header=hdus[0].header
    wcs=WCS(header).celestial
    finite=np.isfinite(im)
    median=np.nanmedian(im)
    noise=1.4826*np.nanmedian(abs(im-median))
    positive=np.where(finite, im-median, -1e9)
    mask=(positive==maximum_filter(positive, size=15)) & (positive>8*noise)
    yy,xx=np.where(mask)
    selected=[]
    for y,x in zip(yy,xx):
        if not 12<=x<im.shape[1]-12 or not 12<=y<im.shape[0]-12:continue
        ra,dec=wcs.pixel_to_world_values(x,y)
        east=(ra-row.ra)*np.cos(np.deg2rad(row.dec))*3600
        north=(dec-row.dec)*3600
        if not 15<np.hypot(east,north)<730:continue
        stamp=im[y-9:y+10,x-9:x+10]
        if not np.isfinite(stamp).all():continue
        sky=np.median(stamp)
        dy,dx=np.mgrid[-9:10,-9:10]
        loc=dx*dx+dy*dy<=5**2
        peak=im[y,x]-sky
        area=1.13309*header['BMAJ']*header['BMIN']*3600**2
        compactness=np.sum((stamp-sky)[loc])/(peak*area)
        if .72<compactness<1.45:
            selected.append([east,north,peak,compactness,noise])
    return np.array(selected).reshape(-1,5)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('ledger',type=Path);parser.add_argument('cache',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();rows=pd.read_csv(args.ledger)
    data={}; hashes=[]
    def acquire(row):
        path=args.cache/f'{row.name}-{row.epoch}.fits'
        raw=download(row,.22,path)
        return row,peaks(path,row),hashlib.sha256(raw).hexdigest()
    fields=rows[rows.name.isin(['SN2019bao','SN2019njv'])]
    with ThreadPoolExecutor(max_workers=4) as pool:
        for row,points,sha in pool.map(acquire,fields.itertuples(index=False)):
            data[row.name,row.epoch]=points;hashes.append(dict(name=row.name,epoch=row.epoch,sha256=sha))
    args.output.mkdir(parents=True, exist_ok=True)
    records=[]
    for name in fields.name.unique():
        for first,second in [(1,3),(1,4),(3,4)]:
            a,b=data[name,first],data[name,second]
            if not len(a) or not len(b):continue
            sep,ids=cKDTree(b[:,:2]).query(a[:,:2])
            for p,q,s in zip(a,b[ids],sep):
                if s>2 or min(p[2]/p[4],q[2]/q[4])<10:continue
                records.append(dict(name=name,first=first,second=second,east=p[0],north=p[1],separation=s,first_mjy=p[2],second_mjy=q[2],first_shape=p[3],second_shape=q[3],ratio=q[2]/p[2]))
    pd.DataFrame(records).to_csv(args.output/'compact-field-sources.csv',index=False)
    pd.DataFrame(hashes).to_csv(args.output/'calibration-fields.csv',index=False)

if __name__=='__main__':main()
