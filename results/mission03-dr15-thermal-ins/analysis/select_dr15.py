"""Select measured ultrasoft 5XMM-DR15 sources and retain the HR2=-1 boundary.

Use the native, uncompressed 5XMM_DR15cat_source_v1.0.fits file for low memory
use. Outputs contain catalogue data, observed band sums, and ranking fields.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.io import fits


N_SOURCE = 818656
PROBS = ('GAIA_MATCH_PROB', 'WISE_MATCH_PROB', 'OM_MATCH_PROB')
COLUMNS = [
    'SRCID', 'IAUNAME', 'RA', 'DEC', 'RADEC_ERR', 'ELL_MAJOR', 'ELL_MINOR',
    'ELL_PA','LII','BII','N_OBS','N_CONTRIB','MJD_FIRST','MJD_LAST',
    'DIST_NN','N_BLEND','EXTENT','SUM_FLAG','STACK_DET_ML','EP_CTS',
    *[f'EP_{band}_FLUX{s}' for band in range(1,6) for s in ('','_ERR')],
    *[f'EP_HR{band}{s}' for band in range(1,5) for s in ('','_ERR')],
    'PN_CTS','M1_CTS','M2_CTS','STACK_GAMMA','STACK_NH','VAR_FLAG',
    'VAR_PROB','FRATIO','FRATIO_ERR','APPROX_SOURCE_VAR',
    'GAIADR3_GMAG','GAIADR3_PARALLAX_OVER_ERROR','WISE_W1MAG','WISE_W2MAG',
    'OM_U_AB_MAG','OM_B_AB_MAG','OM_V_AB_MAG',*PROBS,
    'SPEC_GAMMA_PL','SPEC_NH_PL','SPEC_PVALUE_PL','SPEC_FLAG_PL','REDSHIFT_ZSP'
]


def scalar(value):
    if isinstance(value,(np.floating,float)) and not np.isfinite(value):
        return ''
    if isinstance(value,bytes):
        return value.decode('ascii').strip()
    return value.item() if isinstance(value,np.generic) else value


def run(path, out):
    out.mkdir(parents=True,exist_ok=True)
    with fits.open(path, memmap=True) as hdus:
        d=hdus[1].data
        if len(d)!=N_SOURCE or '421' != str(hdus[1].header['TFIELDS']) or '5XMM ' not in d['IAUNAME'][0]:
            raise ValueError('Expected native 818656-row, 421-column 5XMM-DR15 source table')
        # Promote before multiplying or squaring cgs fluxes.
        bands=np.array([d[f'EP_{i}_FLUX'] for i in range(1,6)], dtype='float64')
        errors=np.array([d[f'EP_{i}_FLUX_ERR'] for i in range(1,6)], dtype='float64')
        fsoft=bands[:3].sum(axis=0)
        esoft=np.sqrt((errors[:3]**2).sum(axis=0))
        counts={'source_rows':len(d)}
        mask=np.ones(len(d),dtype=bool)

        def cut(stage, criterion):
            nonlocal mask
            mask &= criterion
            counts[stage]=int(np.count_nonzero(mask))

        cut('quality_flag_0_or_1',np.isin(d['SUM_FLAG'],[0,1]))
        cut('point_source',d['EXTENT']==0)
        cut('neighbor_at_least_20_arcsec',d['DIST_NN']>=20)
        cut('stack_likelihood_ge_30',d['STACK_DET_ML']>=30)
        cut('counts_ge_100',d['EP_CTS']>=100)
        cut('formal_position_error_gt_0_le_3_arcsec',(d['RADEC_ERR']>0)&(d['RADEC_ERR']<=3))
        cut('valid_soft_flux_ge_1e-14',np.isfinite(bands[:3]).all(axis=0)&(bands[:3]>=0).all(axis=0)&(fsoft>=1e-14))
        hr1,hr2=d['EP_HR1'],d['EP_HR2']
        e1,e2=d['EP_HR1_ERR'],d['EP_HR2_ERR']
        pre_hardness=mask.copy()
        cut('measured_hr1_plus_error_lt_0.4',(abs(hr1)<1)&(e1>0)&(hr1+e1<0.4))
        common=mask.copy()
        measured=(abs(hr2)<1)&(e2>0)&(hr2+e2<-.65)
        boundary=hr2==-1
        hard_ok=np.all(np.isfinite(bands[3:])&np.isfinite(errors[3:])&(bands[3:]>=0)&(errors[3:]>0)&(bands[3:]<3*errors[3:]),axis=0)
        var_ok=~np.asarray(d['VAR_FLAG'],dtype=bool)&~(d['VAR_PROB']<1e-5)
        # Source-row SPECTRA may be false even when observation products exist.
        # Missing match probabilities are retained for independent cone searches.
        has_match=np.logical_or.reduce([d[p]>=0.8 for p in PROBS])
        cameras=np.zeros(len(d),dtype=np.int8)
        for camera in ('PN','M1','M2'):
            cf,ce=d[f'{camera}_1_FLUX'],d[f'{camera}_1_FLUX_ERR']
            cameras+=np.isfinite(cf)&np.isfinite(ce)&(ce>0)&(cf>=3*ce)
        hard_tail=np.sum(bands[2:]+1.645*errors[2:],axis=0)
        rescue=(((abs(hr1)<1)&(e1>0)&(hr1+e1<-.5))|(hr1==-1))
        rescue &= (cameras>=2)&np.all(np.isfinite(errors[2:])&(errors[2:]>0),axis=0)
        rescue &= (hard_tail<.5*fsoft)&~(common&(measured|boundary))
        # A high-column source can have HR1>0.4. Require detections in two
        # cameras in at least one of bands 1 or 2 of each camera.
        soft_cameras=np.zeros(len(d),dtype=np.int8)
        for camera in ('PN','M1','M2'):
            seen=np.zeros(len(d),dtype=bool)
            for band in (1,2):
                f,e=d[f'{camera}_{band}_FLUX'],d[f'{camera}_{band}_FLUX_ERR']
                seen |= np.isfinite(f)&np.isfinite(e)&(e>0)&(f>=3*e)
            soft_cameras+=seen
        hr3,e3=d['EP_HR3'],d['EP_HR3_ERR']
        hr3_measured=(abs(hr3)<1)&np.isfinite(e3)&(e3>0)
        hr2_bright=((abs(hr2)<1)&(e2>0)&(hr2<-.35)&(hr2+e2<-.25)) | (hr2==-1)
        bright=(fsoft>=1.5e-13)&hr2_bright&((~hr3_measured)|(hr3+e3<0))&(soft_cameras>=2)
        bright &= ~(common&(measured|boundary))&~rescue
        for label,branch in [('measured',measured),('boundary',boundary),('rescue',rescue),('bright',bright)]:
            mask=pre_hardness.copy() if label in ('rescue','bright') else common.copy()
            cut(f'{label}_color_gate',branch)
            cut(f'{label}_hard_band_screen',hard_ok)
            cut(f'{label}_no_reported_variability_veto',var_ok)
            indices=np.flatnonzero(mask)
            indices=sorted(indices,key=lambda i:(bool(has_match[i]),-float(fsoft[i]),int(d['SRCID'][i])))
            counts[f'{label}_no_release_match_ge_0.8']=int(np.sum(mask & ~has_match))
            with (out/f'{label}.csv').open('w',newline='') as stream:
                writer=csv.DictWriter(stream,fieldnames=['RANK','HAS_RELEASE_MATCH','F_SOFT','F_SOFT_ERR','N_CAM_B1_SN3','N_CAM_B12_SN3','F345_PLUS_1P645ERR']+COLUMNS,lineterminator='\n')
                writer.writeheader()
                for j,i in enumerate(indices,1):
                    fields={}
                    for k in COLUMNS:
                        value=d[k][i]
                        null=hdus[1].columns[k].null
                        fields[k]='' if null is not None and value==null else scalar(value)
                    writer.writerow(dict(RANK=j,HAS_RELEASE_MATCH=int(has_match[i]),F_SOFT=fsoft[i],F_SOFT_ERR=esoft[i],N_CAM_B1_SN3=int(cameras[i]),N_CAM_B12_SN3=int(soft_cameras[i]),F345_PLUS_1P645ERR=scalar(hard_tail[i]),**fields))
        (out/'cutflow.json').write_text(json.dumps(counts,indent=2)+'\n')
        print(json.dumps(counts,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('catalogue',type=Path)
    parser.add_argument('output_directory',type=Path)
    args=parser.parse_args()
    run(args.catalogue,args.output_directory)
