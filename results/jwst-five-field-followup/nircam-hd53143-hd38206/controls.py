"""Conditional F444W source counts and Bobcat photometry at the two stellar distances.
Reuse the input reader and the original GOODS-S field calculation as a control on
the mask. Extends its flux sampling to each threshold crossing and this annulus.
"""
import argparse,csv,gzip,io,json,sys,tarfile
from pathlib import Path
import numpy as np
from scipy.stats import chi2

def run(parent,cache,out):
    sys.path.insert(0,str(parent.resolve()/'vetting/background'))
    import background
    background.load_data(cache,out)
    # Retain the original calculation as a reproducibility check for the new threshold list.
    background.field(cache,out)
    with gzip.open(cache/'gskronc.dat.gz','rt') as k,gzip.open(cache/'gsflags.dat.gz','rt') as f:
        dat=[]
        for a,b in zip(k,f,strict=True):
            assert a[:11]==b[:11]
            if min(float(b[244:256]),float(b[335:347]))>0:dat.append([float(a[12:23]),float(a[24:35]),float(a[465:477])/1000])
    ra,dec,flux=np.array(dat).T
    quarter=np.floor(np.column_stack([(ra-53.12)*3600*np.cos(np.deg2rad(-27.8)),(dec+27.8)*3600])/15).astype(int)
    coords,n=np.unique(quarter,axis=0,return_counts=True);qcounts=dict(zip(map(tuple,coords),n))
    cells=quarter//2
    good={tuple(p) for p in np.unique(cells,axis=0) if all(qcounts.get(tuple(2*p+np.array([a,b])),0)>=5 for a in [0,1] for b in [0,1])}
    use=np.array([tuple(p) in good for p in cells]);area=len(good)/4
    for row in csv.DictReader((out/'control_counts.csv').open()):
        if row['selection']=='all':
            assert int((use&(flux>=float(row['min_F444W_uJy']))).sum())==int(row['count']) and area==float(row['control_area_arcmin2'])
    sets={}
    with tarfile.open(cache/'bobcat.tar.gz') as t:
        for met in ['-0.5','+0.0','+0.5']:
            lines=t.extractfile('photometry_tables/flux_table_JWST'+met).read().decode().splitlines();j=next(i for i,s in enumerate(lines) if s.strip().startswith('Teff'));keys=lines[j].split();bands=keys[keys.index('F070W'):]
            tab=np.loadtxt(io.StringIO('\n'.join(s for s in lines[j+1:] if s.strip() and s.strip()[0].isdigit() and '*' not in s)))
            sets[met]=(tab,bands)
    rows=[]
    for p in csv.DictReader((out/'paired_color_morphology.csv').open()):
        visit=p['visit'];lw=float(p['long_flux_uJy']);b=p['filter'];isa=visit.startswith('hd53143');dist=1000/(54.52523915853711 if isa else 14.146997)
        magcolor=[]
        for met,(tab,bands) in sets.items():
            vals=1000*(10/dist)**2*10**tab[:,[6+bands.index(b),6+bands.index('F444W')]]
            mask=(vals[:,1]>=.75*lw)&(vals[:,1]<=1.25*lw)
            for v in vals[mask]:magcolor.append(v[0]/v[1])
        n=int((use&(flux>=lw)).sum());angle=np.pi*(5.5**2-2.5**2)/3600
        rows.append(dict(visit=visit,x=p['x'],y=p['y'],flux_uJy=lw,blank_GS_count=n,control_area_arcmin2=area,expected_full_annulus_interlopers=n/area*angle,
                         poisson_low=.5*chi2.ppf(.025,2*n)*angle/area if n else 0.,poisson_high=.5*chi2.ppf(.975,2*n+2)*angle/area,
                         bobcat_age_unrestricted_nodes=len(magcolor),short_filter=b,lw_match_tolerance_fraction=.25,
                         model_min_AB_color=-2.5*np.log10(max(magcolor)) if magcolor else '',model_max_AB_color=-2.5*np.log10(min(magcolor)) if magcolor else ''))
    with (out/'controls_and_models.csv').open('w') as f:w=csv.DictWriter(f,list(rows[0]));w.writeheader();w.writerows(rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('cache',type=Path);p.add_argument('out',type=Path);a=p.parse_args();run(a.parent,a.cache,a.out)
