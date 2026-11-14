"""Refit primary-secondary separation with OGLE-III and Gaia DR3 catalog periods."""
import json
from pathlib import Path
import numpy as np
from timing import SPLITS, fit, signed_phase, profile

BASE=Path(__file__).resolve().parent
def main():
    initial_data=json.loads((BASE/'measurements.json').read_text())
    gaia=json.loads((BASE/'gaia-validation/validation.json').read_text())
    periods={'OGLE-IV':initial_data['catalog_period_d'],'OGLE-III':1.1696692,
        'Gaia DR3':gaia['gaia_published_model']['period_days']}
    zero=initial_data['catalog_epoch_hjd']
    results={}
    for label,period in periods.items():
        epochs=[]
        for name,source,low,high in SPLITS:
            t,m,e=np.loadtxt(BASE/'timing-inputs'/source).T
            if source=='ogle3.dat': t+=2450000
            use=np.isfinite(t+m+e)&(e>0)&(e<=.1)&(t>=low)&(t<high)
            t,m,e=t[use],m[use],e[use]
            ph=signed_phase((t-zero)/period)
            # Allow a common phase shift accumulated from changing the fold period.
            phase_shift=float(signed_phase((t.mean()-zero)*(1/period-1/periods['OGLE-IV'])))
            ph=signed_phase(ph-phase_shift)
            par=fit(ph,m,e,[np.quantile(m,.25),.4,0,.023,.4,.51,.023],'gaussian')
            centers=[par[2],par[5]]
            core=[abs(signed_phase(ph-center))<.05 for center in centers]
            mean_t=[float(t[s].mean()) for s in core]
            epochs.append({'name':name,'separation':par[5]-par[2],
                'core_rows':[int(s.sum()) for s in core], 'core_mean_hjd':mean_t,
                'secondary_minus_primary_mean_core_date_days':mean_t[1]-mean_t[0]})
        results[label]={'period_days':period,'epochs':epochs,
            'OGLE-IV_last_minus_first':epochs[5]['separation']-epochs[2]['separation'],
            'OGLE-III_last_minus_first':epochs[1]['separation']-epochs[0]['separation']}
    og=results['OGLE-IV']['epochs']
    dtchange=og[5]['secondary_minus_primary_mean_core_date_days']-og[2]['secondary_minus_primary_mean_core_date_days']
    maxfreq=max(abs(1/p-1/periods['OGLE-IV']) for p in periods.values())
    (BASE/'period_check.json').write_text(json.dumps({'catalog_periods':results,
         'approximation_note':'For the catalog ephemeris only, calculate mean date of observations within 0.05 cycles of each fitted center. The change in their mean-date difference multiplied by catalog frequency discrepancy approximates a common period error artifact. This is a diagnostic; an actual refit at each alternative period appears above.',
         'change_in_mean_core_sampling_date_difference_early_to_late_d':dtchange,
         'max_abs_catalog_frequency_difference_per_day':maxfreq,
         'approximate_magnitude_common_period_error_shift_cycles':abs(dtchange)*maxfreq,
         'OGLE_III_period_table':'https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/ecl.dat'},indent=2)+'\n')
if __name__=='__main__': main()
