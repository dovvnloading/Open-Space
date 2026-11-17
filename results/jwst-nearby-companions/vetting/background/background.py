"""JADES GOODS-S surface counts and Sonora Bobcat comparison for the L 788-34 source.

Run python background.py CACHE OUTPUT. Statistical control is defined independently
of flux: at least twenty sources with positive F200W and F444W exposure in each
30x30 arcsec cell. Counts and the geometrical denominator use the identical cells.
"""
import argparse
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import tarfile
from urllib.request import urlopen
import numpy as np
from scipy.stats import chi2

FILES = {
    "gskronc.dat.gz":"https://vizier.cfa.harvard.edu/ftp/cats/V/159/gskronc.dat.gz",
    "gsflags.dat.gz":"https://vizier.cfa.harvard.edu/ftp/cats/V/159/gsflags.dat.gz",
    "bobcat.tar.gz":"https://zenodo.org/api/records/5063476/files/evolution_and_photometery.tar.gz/content",
}


def load_data(cache, out):
    cache.mkdir(parents=True, exist_ok=True)
    receipt = []
    for name,url in FILES.items():
        p=cache/name
        if not p.exists():
            with urlopen(url, timeout=120) as r:
                p.write_bytes(r.read())
        receipt.append(dict(filename=name,url=url,bytes=p.stat().st_size,
                            sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    (out/'inputs.json').write_text(json.dumps(receipt,indent=2)+'\n')


def field(cache, out):
    # CDS V/159 byte ranges are one-indexed inclusive in the ReadMe.
    # KRON_CONV includes a F444W-homogenized point-spread function in all bands.
    data=[]
    with gzip.open(cache/"gskronc.dat.gz","rt") as k,gzip.open(cache/"gsflags.dat.gz","rt") as g:
        for kr,fl in zip(k,g,strict=True):
            assert kr[0:11] == fl[0:11]
            if min(float(fl[244:256]),float(fl[335:347])) <= 0:
                continue
            ra,dec = float(kr[12:23]),float(kr[24:35])
            sw,lw = float(kr[192:204]) / 1000, float(kr[465:477]) / 1000
            data.append((ra,dec,sw,lw,int(fl[790:796])))
    ra,dec,sw,lw,flag=np.array(data).T
    x=(ra-53.12)*np.cos(np.radians(-27.8))*3600
    y=(dec+27.8)*3600
    # Require dense coverage in four quarters: each 30" square has at least
    # five dual-exposure catalogued sources in each independent 15" quadrant.
    cells=np.floor(np.stack([x,y],axis=1)/30).astype(int)
    quarter=np.floor(np.stack([x,y],axis=1)/15).astype(int)
    uq, counts=np.unique(quarter,axis=0,return_counts=True)
    qdict={tuple(q):count for q,count in zip(uq,counts)}
    accepted={tuple(cell) for cell in np.unique(cells,axis=0)
              if all(qdict.get(tuple(cell*2+[a,b]),0)>=5 for a in [0,1] for b in [0,1])}
    use=np.array([tuple(c) in accepted for c in cells])
    area=len(accepted)*.25
    color=np.full(len(sw),np.nan)
    good=(sw>0)&(lw>0)
    color[good]=-2.5*np.log10(sw[good]/lw[good])
    rows=[]
    for level in [20,25,35,40]:
        for label,mask in [('all',np.ones(len(sw),bool)),('broad_blue_AB', (color>=-1.5)&(color<=-.5)), ('catalogue_galaxy', flag==0), ('blue_catalogue_galaxy', (color>=-1.5)&(color<=-.5)&(flag==0))]:
            n=int(np.sum(use & (lw>=level) & mask))
            density=n/area
            rows.append(dict(min_F444W_uJy=level,selection=label,count=n,control_area_arcmin2=area,
                             number_per_arcmin2=density,conditional_expect_r3p21=density*np.pi*3.21**2/3600,
                             poisson_lower95_per_arcmin2=chi2.ppf(.025,2*n)/2/area if n else 0.,
                             poisson_upper95_per_arcmin2=chi2.ppf(.975,2*n+2)/2/area))
    with (out/'control_counts.csv').open('w') as f:
        wr=csv.DictWriter(f,list(rows[0]));wr.writeheader();wr.writerows(rows)


def models(cache,out):
    d=1000/138.228
    rows=[]
    with tarfile.open(cache/'bobcat.tar.gz') as tar:
        for metallicity in ['-0.5','+0.0','+0.5']:
            lines=tar.extractfile('photometry_tables/flux_table_JWST'+metallicity).read().decode().splitlines()
            hdr=next(i for i,s in enumerate(lines) if s.strip().startswith('Teff'))
            keys=lines[hdr].split()
            bands=keys[keys.index('F070W'):]
            tab=np.loadtxt(io.StringIO('\n'.join(s for s in lines[hdr+1:]
                                                  if s.strip() and s.strip()[0].isdigit() and '*' not in s)))
            cols=[6+bands.index(k) for k in ['F200W','F444W']]
            flux=1000*(10/d)**2*10**tab[:,cols]  # catalog parallax; mJy at 10pc -> observed µJy
            for model,(sw,lw) in zip(tab,flux):
                if 17<=lw<=45:
                    rows.append(dict(metallicity=metallicity,teff=int(model[0]),logg=model[1],F200W_uJy=sw,
                                     F444W_uJy=lw,color_AB=-2.5*np.log10(sw/lw)))
    with (out/'models_at_host_distance.csv').open('w') as f:
        wr=csv.DictWriter(f,list(rows[0]));wr.writeheader();wr.writerows(rows)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('cache',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();a.out.mkdir(exist_ok=True,parents=True)
    load_data(a.cache,a.out);field(a.cache,a.out);models(a.cache,a.out)
