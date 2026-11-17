"""DR1-only fit and zero-refit QR3 transfer for the UGC 8856 3.3-micron PAH alternative.

Target and check chosen after viewing coadd/native outcomes. This is retrospective
shape validation of a common astrophysical hypothesis, not an independent line ID.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear

ID='J135607.93+300452.8'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('study',type=Path)
    args=p.parse_args()
    field=args.study/'wide'
    data=pd.read_csv(args.study/'visits/profile_J135607.csv')
    # Match the native-profile experiment's window, retaining every good
    # field exposure in it, independently of measured source brightness.
    data=data.loc[(data.status=='ok') & data.wave_um.between(3.07,3.73)].copy()
    redshifts=pd.read_csv(field/'catalogs/model_redshifts.csv').set_index('source_id')
    z=redshifts.loc[ID,'redshift']
    wave=data.wave_um.to_numpy()
    center=3.29*(1+z)
    # Analytic gaussian LSF with per-exposure CBAND. Peak response includes
    # narrower higher-amplitude correspondence at a fixed line integral.
    width=np.hypot(.045*(1+z),data.bandwidth_um)
    profile=(.045*(1+z)/width) * np.exp(-4*np.log(2)*((wave-center)/width)**2)
    continuum=np.column_stack([np.ones(len(wave)),wave-center])
    sigma=np.hypot(data.sigma,.03*abs(data.flux)).to_numpy()
    train=data.release.to_numpy()=='dr1'
    y=data.flux.to_numpy()
    out={'source_id':ID,'z':z,'assumed_rest_um':3.29,'assumed_rest_fwhm_um':.045,
         'expected_um':center,'training_n':int(train.sum()),'transfer_n':int((~train).sum())}
    for name,mat in [('linear',continuum),('linear_plus_pah',np.column_stack([continuum,profile]))]:
        lo=np.full(mat.shape[1],-np.inf)
        if name=='linear_plus_pah':lo[-1]=0
        fit=lsq_linear(mat[train]/sigma[train,None],y[train]/sigma[train],
                       bounds=(lo,np.full(mat.shape[1],np.inf)))
        pred=mat@fit.x
        data[name+'_prediction']=pred
        out[name]={'dr1_fit':fit.x.tolist(),'dr1_diagnostic_chi2':float(sum(((pred[train]-y[train])/sigma[train])**2)),
                   'qr3_no_refit_chi2':float(sum(((pred[~train]-y[~train])/sigma[~train])**2))}
    outdir=args.study/'visits/wide'
    data[['release','mjd','wave_um','bandwidth_um','flux','sigma','linear_prediction','linear_plus_pah_prediction']].to_csv(
        outdir/'pah_transfer.csv',index=False,float_format='%.8g')
    (outdir/'pah_transfer.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':
    main()
