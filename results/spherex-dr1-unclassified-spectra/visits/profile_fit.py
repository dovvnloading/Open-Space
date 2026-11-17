"""Estimate the D4 band center independently in DR1 and QR3, using the native scan."""
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parent
REDSHIFT = .0310622
EXPECTED = 3.30*(1+REDSHIFT)
rows = list(csv.DictReader((ROOT/'profile_J135607.csv').open()))
result = dict(source_id='J135607.93+300452.8',external_redshift=REDSHIFT,expected_3p30_obs_um=EXPECTED)
for release in ('dr1','qr3'):
    allrows=[r for r in rows if r['release']==release]
    accepted=[r for r in allrows if r['status']=='ok' and 3.07<=float(r['wave_um'])<=3.73]
    group=dict(retrieved=len(allrows),accepted_in_window=len(accepted))
    if len(accepted)>=6:
        w,b,y,e=[np.array([float(r[key]) for r in accepted]) for key in ('wave_um','bandwidth_um','flux','sigma')]
        noise=np.hypot(e,.03*np.abs(y))
        def model(params):
            intercept,slope,amplitude,center=params
            # Approximate both the detector response and intrinsic aromatic complex
            # as Gaussians. Fix intrinsic FWHM = 0.04 um in the rest frame.
            sigma_inst=b/2.355
            sigma_intr=.04*(1+REDSHIFT)/2.355
            sigma=np.hypot(sigma_inst,sigma_intr)
            return intercept+slope*(w-3.4)+amplitude*np.exp(-.5*((w-center)/sigma)**2)/sigma
        solution=least_squares(lambda p:(model(p)-y)/noise,
             [np.median(y),0.,.2,EXPECTED],bounds=([-np.inf,-np.inf,0.,3.25],[np.inf,np.inf,np.inf,3.55]),
             loss='soft_l1')
        center=float(solution.x[3])
        group.update(date_count=len({r['date_obs'][:10] for r in accepted}),center_um=center,
                     rest_center_um=center/(1+REDSHIFT),difference_from_expected_um=center-EXPECTED,
                     median_response_fwhm_um=float(np.median(b)),offset_in_native_widths=float((center-EXPECTED)/np.median(b)),
                     fit_parameters=solution.x.tolist(),min_wave=float(w.min()),max_wave=float(w.max()))
    result[release]=group
(ROOT/'profile_fit.json').write_text(json.dumps(result,indent=2)+'\n')
