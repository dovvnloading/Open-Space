"""Measure off-season recurrence for catalog-short alerts and named follow-ups."""
from pathlib import Path
import csv,json
from collections import Counter
import numpy as np
from fetch import OUT,write_csv
CACHE=Path('.survey-cache')
NAMES={'MOA-2024-BLG-111','MOA-2025-BLG-0093','MOA-2025-BLG-0198','MOA-2025-BLG-0211','MOA-2025-BLG-0301','PRIME-2025-BLG-0415','PRIME-2024-BLG-171'}

def jd_year(year):
    return 2440587.5+(np.datetime64(f'{year}-01-01')-np.datetime64('1970-01-01')).astype(int)


def main():
    phot=list(csv.DictReader((OUT/'photometry.csv').open()))
    selection=[]
    yearly=[]
    for r in phot:
        try: te=float(r['page_te_d']); u=float(r['page_u0'])
        except (TypeError,ValueError):continue
        if not (0<te<.5 or r['event'] in NAMES):continue
        t,f,e=np.load(CACHE/'curves'/(r['event']+'.npz'))['data'].T
        # Numerically stable evaluation of A(u)-1; half-excess at u_half solves
        # A=1+0.5*(A(u0)-1), using u^2=2/[sqrt(A^2-1)(A+sqrt(A^2-1))].
        if u>0:
            excess=4/(u*np.sqrt(u*u+4)*(u*u+2+u*np.sqrt(u*u+4)))
            m=excess/2; root=np.sqrt(m*(m+2)); uh2=2/(root*(1+m+root))
            width=2*te*np.sqrt(max(0,uh2-u*u))
        else: width=None
        prior=[]
        for year in range(2005,2027):
            ix=(t>=jd_year(year))&(t<jd_year(year+1))
            tt,ff,ee=t[ix],f[ix],e[ix]
            if len(tt)<3:continue
            med=np.median(ff);mad=1.4826*np.median(abs(ff-med))
            z=(ff-med)/np.hypot(ee,mad)
            tri=np.minimum.reduce([z[:-2],z[1:-1],z[2:]])
            tri[tt[2:]-tt[:-2]>1]=-np.inf
            if not np.isfinite(tri).any():continue
            k=np.argmax(tri)
            row=dict(event=r['event'],year=year,n=len(tt),median_flux=round(med,3),mad_sigma=round(mad,3),
                points_gt5=int((z>5).sum()),triple_weakest_z=round(tri[k],3),
                triple_middle_jd=round(tt[k+1],6),triple_span_d=round(tt[k+2]-tt[k],6))
            yearly.append(row)
            if year!=int(r['event'].split('-')[1]) and tri[k]>5:prior.append(str(year))
        selection.append(dict(event=r['event'],assessment=r['assessment'],catalog_te_d=te,catalog_u0=u,
                              catalog_pspl_fwhm_d=round(width,5) if width is not None else '',
                              points_within_half_day=r['points_within_half_day'],
                              offyear_three_gt5_within_day='|'.join(prior)))
    write_csv(OUT/'short_alert_checks.csv',selection)
    write_csv(OUT/'seasonal_recurrence.csv',yearly)
    # Avoid a field-count bias when stating a survey cadence: publish per-alert counts
    # and explicitly summarize by field/season. This is conditional on sampled alerts.
    for key in ('input_bytes','raw_rows','usable'):
        assert all(int(r[key])>0 for r in phot)
    grouped={}
    for r in phot:
        group=tuple(r['event'].split('-')[:2])+(r['field'],)
        grouped.setdefault(group,[]).append(r)
    out=[]
    for (survey,year,field),rows in sorted(grouped.items()):
        gap=[float(r['positive_gap_lt_halfday_median_min']) for r in rows if r['positive_gap_lt_halfday_median_min']]
        out.append(dict(survey=survey,season=year,field=field,alerts=len(rows),usable_epoch_sum=sum(int(r['usable']) for r in rows),
                        median_of_within_night_gap_medians_min=round(float(np.median(gap)),3) if gap else ''))
    write_csv(OUT/'fields.csv',out)

if __name__=='__main__':main()
