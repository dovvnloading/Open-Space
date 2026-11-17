"""Audit RV constancy and stellar luminosity constraints in two transverse candidates."""
from pathlib import Path
import csv, json,sys
import pandas as pd
import numpy as np
from scipy.stats import chi2
from astropy.io import fits
D=Path(__file__).resolve().parent
# Koposov et al. 2026: bright-program multi-night RV uncertainty floor is ~1 km/s.
stars=pd.read_csv(D/'epochs.csv',dtype={'targetid':str})
coadds=pd.read_csv(D.parent/'queries/coadds.csv',dtype={'targetid':str,'source_id':str})
rows=[]
for tid,group in stars.groupby('targetid'):
    parents=coadds[coadds.targetid==tid]
    q=group[(group.rvs_warn==0)&(group.fiberstatus==0)&(group.success==True)&(group.vrad_err>0)]
    # Within-night statistical error on weighted mean, plus one shared calibration floor per night.
    nights=[]
    for night,h in q.groupby('night'):
        w=1/h.vrad_err**2;nights.append((int(night),float(sum(w*h.vrad)/sum(w)),float(np.sqrt(1/sum(w)+1.0))))
    a=np.array(nights); w=1/a[:,2]**2;mean=sum(a[:,1]*w)/sum(w);c2=sum(w*(a[:,1]-mean)**2);dof=len(a)-1
    rows.append({'targetid':tid,'source_id':parents.source_id.iloc[0],'exposures':len(q),'distinct_nights':len(nights),'baseline_days':float(q.mjd.max()-q.mjd.min()),'nights_mean_error':nights,'constant_rv':float(mean),'chi2':float(c2),'dof':dof,'p_constancy':float(chi2.sf(c2,dof))})
(D/'multi_epoch_tests.json').write_text(json.dumps(rows,indent=2)+'\n')
s=fits.open(sys.argv[1],memmap=True)[1].data
fields='SOURCE_ID TARGETID FLAG_GOOD TEFF_COR LOGG_COR G0 BP0 RP0 MG_16 MG_50 MG_84 QFLAG_INPUT QFLAG_COMMENTS'.split()
with (D/'luminosity.csv').open('w') as out:
    wr=csv.writer(out);wr.writerow(fields+['BP0_RP0','spdist_kpc'])
    for r in s[np.isin(s['TARGETID'],[int(x['targetid']) for x in rows])]:
        wr.writerow([r[k] for k in fields]+[r['BP0']-r['RP0'],10**((r['G0']-r['MG_50']-10)/5)])
