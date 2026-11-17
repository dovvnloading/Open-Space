"""Archive measured tables and derive area, completeness and calibration diagnostics."""
import argparse
import gzip
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd

def summarize(run,out,pixscale):
    out.mkdir(parents=True,exist_ok=True)
    for name in ['detections','calibrators','coverage','injections']:
        raw=(run/f'{name}.csv').read_bytes()
        (out/f'{name}.csv.gz').write_bytes(gzip.compress(raw,mtime=0))
    shutil.copyfile(run/'model.npz',out/'model.npz')
    d=pd.read_csv(run/'detections.csv');c=pd.read_csv(run/'coverage.csv')
    i=pd.read_csv(run/'injections.csv').groupby(['mag_ab','sign'])[['trials','recovered']].sum().reset_index()
    i['fraction']=i.recovered/i.trials
    i.to_csv(out/'recovery.csv',index=False)
    d.query('screen_pass == 1').to_csv(out/'screened.csv',index=False)
    pooled=i.groupby('mag_ab')[['trials','recovered']].sum()
    frac=pooled.recovered/pooled.trials
    src=pd.read_csv(run/'calibrators.csv').query('used_psf == 1')
    shift=src[['dx','dy']].median()
    depth={}
    for frac0 in [.5,.9]:
        below=np.flatnonzero(frac.values<frac0)
        if not len(below) or below[0]==0: continue
        j=below[0]; depth[str(frac0)]=float(np.interp(frac0,frac.iloc[[j,j-1]],frac.index[[j,j-1]]))
    data=dict(scale_arcsec=pixscale,area_arcmin2=c.n_pix.sum()*pixscale**2/3600,pixels=int(c.n_pix.sum()),
              peaks=len(d),shape_pass=int(d.screen_pass.sum()),calibrators=len(src),dx_pix=float(shift.dx),dy_pix=float(shift.dy),
              radial_registration_scatter_median_pix=float(np.median(np.hypot(src.dx-shift.dx,src.dy-shift.dy))),
              median_calibrator_ratio=float(src.ratio.median()),empirical_psf_AB_depth=depth)
    (out/'stats.json').write_text(json.dumps(data,indent=2)+'\n');print(out,data,flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('run',type=Path);p.add_argument('out',type=Path);p.add_argument('arcsec_per_pix',type=float)
    a=p.parse_args();summarize(a.run,a.out,a.arcsec_per_pix)
