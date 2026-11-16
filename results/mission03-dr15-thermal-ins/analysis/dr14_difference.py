"""Search for soft DR15 sources without a DR14 slim position within 15 arcsec."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from scipy.spatial import cKDTree
from select_dr15 import COLUMNS, PROBS, scalar


def sphere(ra,dec):
    r,d=np.radians(ra),np.radians(dec)
    return np.column_stack((np.cos(d)*np.cos(r),np.cos(d)*np.sin(r),np.sin(d)))


def run(dr15,dr14,out):
    out.mkdir(exist_ok=True,parents=True)
    with fits.open(dr15,memmap=True) as new, fits.open(dr14,memmap=True) as old:
        d=new[1].data; p=old[1].data
        if len(d)!=818656 or len(p)!=692109:
            raise ValueError('Expected source-level DR15 and DR14 slim inputs')
        cuts={'source_rows':len(d)}
        keep=np.ones(len(d),dtype=bool)
        def cut(label,condition):
            keep[:]&=condition
            cuts[label]=int(keep.sum())
        cut('flag_0_or_1',np.isin(d['SUM_FLAG'],[0,1]))
        cut('point_like',d['EXTENT']==0)
        cut('nn_ge_20',d['DIST_NN']>=20)
        cut('stack_ml_ge_30',d['STACK_DET_ML']>=30)
        cut('ep_cts_ge_70',d['EP_CTS']>=70)
        cut('stat_position_gt_0_le_3',(d['RADEC_ERR']>0)&(d['RADEC_ERR']<=3))
        f=np.array([d[f'EP_{i}_FLUX'] for i in range(1,6)],dtype=float)
        e=np.array([d[f'EP_{i}_FLUX_ERR'] for i in range(1,6)],dtype=float)
        fsoft=f[:3].sum(axis=0)
        es=np.sqrt(np.sum(e[:3]**2,axis=0))
        cut('finite_soft_ge_3e-14',(np.isfinite(f[:3])&(f[:3]>=0)).all(axis=0)&(fsoft>=3e-14))
        hr2,e2=d['EP_HR2'],d['EP_HR2_ERR']
        cut('hr2_upper_lt_m0.5_or_boundary',((abs(hr2)<1)&(e2>0)&(hr2+e2<-.5)) | (hr2==-1))
        cameras=np.zeros(len(d),dtype=np.int8)
        for camera in ('PN','M1','M2'):
            seen=np.zeros(len(d),dtype=bool)
            for band in (1,2):
                cf,ce=d[f'{camera}_{band}_FLUX'],d[f'{camera}_{band}_FLUX_ERR']
                seen |= np.isfinite(cf)&np.isfinite(ce)&(ce>0)&(cf>=3*ce)
            cameras+=seen
        cut('two_cameras_soft_ge_3sigma',cameras>=2)
        cut('hard_band_veto',np.all(np.isfinite(f[3:])&np.isfinite(e[3:])&(f[3:]>=0)&(e[3:]>0)&(f[3:]<3*e[3:]),axis=0))
        cut('no_intra_or_large_significant_inter_variability',~np.array(d['VAR_FLAG'],dtype=bool)&(~(d['VAR_PROB']<1e-5) | ((d['FRATIO']>0)&(d['FRATIO']<2))))
        idx=np.flatnonzero(keep)
        old_xyz=sphere(p['SC_RA'].astype(float),p['SC_DEC'].astype(float))
        dist,nn=cKDTree(old_xyz).query(sphere(d['RA'][idx],d['DEC'][idx]))
        arc=7200*np.degrees(np.arcsin(np.minimum(1,dist/2)))
        rows=[(i,int(j),float(theta)) for i,j,theta in zip(idx,nn,arc) if theta>15]
        cuts['no_dr14_slim_source_within_15_arcsec']=len(rows)
        has=np.logical_or.reduce([d[k]>=0.8 for k in PROBS])
        rows.sort(key=lambda row:(bool(has[row[0]]),-float(fsoft[row[0]]),int(d['SRCID'][row[0]])))
        cuts['no_release_match_ge_0.8']=sum(not bool(has[i]) for i,_,_ in rows)
        with (out/'difference.csv').open('w',newline='') as stream:
            extra=['RANK','HAS_RELEASE_MATCH','F_SOFT','F_SOFT_ERR','N_CAM_B12_SN3','DR14_NEAREST_SRCID','DR14_NEAREST_SEP_ARCSEC']
            w=csv.DictWriter(stream,fieldnames=extra+COLUMNS,lineterminator='\n')
            w.writeheader()
            for k,(i,j,theta) in enumerate(rows,1):
                row={}
                for c in COLUMNS:
                    v=d[c][i]; null=new[1].columns[c].null
                    row[c]='' if null is not None and v==null else scalar(v)
                w.writerow(dict(row,RANK=k,HAS_RELEASE_MATCH=int(has[i]),F_SOFT=fsoft[i],
                    F_SOFT_ERR=es[i],N_CAM_B12_SN3=int(cameras[i]),DR14_NEAREST_SRCID=int(p['SRCID'][j]),DR14_NEAREST_SEP_ARCSEC=theta))
        (out/'difference_cutflow.json').write_text(json.dumps(cuts,indent=2)+'\n')
        print(cuts)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('dr15','dr14','output'): parser.add_argument(name,type=Path)
    a=parser.parse_args()
    run(a.dr15,a.dr14,a.output)
