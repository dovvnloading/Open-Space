"""Conditional shape comparison on the same DR1 measurements for all families.

Usage: python template_check.py fits.csv output.csv --templates DIR
The three immutable template archives are downloaded and verified at first use.
"""
import argparse
import csv
import hashlib
import io
import tarfile
import zipfile
from pathlib import Path
from urllib.request import urlopen

from frame_inputs import read_frames

import numpy as np

FILES = {
    'known.zip': ('https://zenodo.org/records/19051216/files/spiff_known_bds_binned_csv.zip?download=1',
                  'd8d16ee156074a1ade74cdd9c024e0119347eb11781e83b849680001f1e916ec'),
    'templates.zip': ('https://zenodo.org/records/19051216/files/spiff_templates_raw_csv.zip?download=1',
                      '342bcef7d550acdacb4835491c1bf1ca4f2a9fbcce3d397bc681719e0ad94e9f'),
    'swire_library.tar.gz': ('https://www.iasf-milano.inaf.it/~polletta/templates/data/swire_library.tar.gz',
                            '34cce03a16fde1b3f0da0a39b84ac38ba5c1d526697200d72beed0fed3f9b48c'),
}


def templates(directory):
    directory.mkdir(parents=True,exist_ok=True)
    for name,(url,expected) in FILES.items():
        path = directory/name
        if not path.exists():
            with urlopen(url,timeout=90) as response:
                path.write_bytes(response.read())
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError(f'Template archive SHA-256 mismatch: {name}')
    records=[]
    with zipfile.ZipFile(directory/'templates.zip') as z:
        for subtype in ('T7','T8','T9'):
            name = f'spiff_templates/raw_csv/{subtype}_raw_spherex_standard.csv'
            x,y,_ = np.loadtxt(io.BytesIO(z.read(name)),delimiter=',',skiprows=1).T
            records.append(('late_T_or_Y',subtype,0.,x,y*x*x))
    with zipfile.ZipFile(directory/'known.zip') as z:
        for source,label in [('WISEPA_J173835','Y0_J1738'),('WISE_J035000','Y1_J0350')]:
            name, = (n for n in z.namelist() if n.split('/')[-1].startswith(source))
            x,y,_ = np.loadtxt(io.BytesIO(z.read(name)),delimiter=',',skiprows=1).T
            records.append(('late_T_or_Y',label,0.,x,y*x*x))
    with tarfile.open(directory/'swire_library.tar.gz') as tar:
        for member in tar.getmembers():
            if not member.name.endswith('.sed'):
                continue
            # SWIRE wavelength in Angstrom and F_lambda, convert shape to F_nu.
            x,y = np.loadtxt(tar.extractfile(member)).T[:2]
            x = x / 1e4
            fnu = y*x*x
            for z in np.linspace(0,4,41):
                label=member.name.split('_')[0]
                records.append(('galaxy_or_AGN',label,z,x*(1+z),fnu))
    return records


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('measurements',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--templates',type=Path,required=True)
    args=p.parse_args()
    grid=templates(args.templates)
    fitrows=[r for r in read_frames([args.measurements]) if r['status']=='fit' and 3.1<=float(r['wavelength_um'])<=4.97]
    out=[]
    for target in sorted({r['target'] for r in fitrows}):
        data=[r for r in fitrows if r['target']==target]
        w,d,f,sig=np.array([[float(r[k]) for k in ('wavelength_um','bandwidth_um','flux_ujy','sigma_ujy')] for r in data]).T
        # Uniform top-hat approximation to each detector response.
        sw=w[:,None]+d[:,None]*np.linspace(-.5,.5,7)
        models=[]
        for family,name,z,x,y in grid:
            valid=np.isfinite(x)&np.isfinite(y)
            xx,yy=x[valid],y[valid]
            curve=np.interp(sw,xx,yy,left=np.nan,right=np.nan)
            over=np.searchsorted(xx,sw).clip(1,len(xx)-1)
            # Reject material holes in empirical brown-dwarf spectra.
            if family=='late_T_or_Y':
                curve[np.abs(xx[over]-xx[over-1]) > .20]=np.nan
            models.append(np.trapezoid(curve,axis=1)/6)
        models=np.array(models)
        keep=np.all(np.isfinite(models),axis=0)
        if not np.any(keep):continue
        data_snr=f[keep]/sig[keep];templates_snr=models[:,keep]/sig[None,keep]
        scale=np.maximum(0,(templates_snr@data_snr)/(templates_snr**2).sum(axis=1))
        chi=((data_snr-scale[:,None]*templates_snr)**2).sum(axis=1)
        for family in ('late_T_or_Y','galaxy_or_AGN'):
            indices=np.array([j for j,r in enumerate(grid) if r[0]==family])
            j=indices[np.argmin(chi[indices])]
            _,name,z,_,_=grid[j]
            out.append(dict(target=target,n_frames=int(keep.sum()),n_before_common=int(len(w)),
                        family=family,template=name,z=round(z,2),chi2=round(chi[j],3),
                        null_chi2=round(float(data_snr@data_snr),3),
                        n_templates=len(indices)))
    with args.output.open('w') as f:
        wr=csv.DictWriter(f,fieldnames=list(out[0]));wr.writeheader();wr.writerows(out)


if __name__=='__main__':
    main()
