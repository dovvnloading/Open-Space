"""Summarize the COSMOS-022 exposure history and the host-conditional difference luminosity."""
import json
import numpy as np
import pandas as pd
from astropy.cosmology import Planck18
from astropy.time import Time
from campaign import HERE

KEY='cosmos_medium-022'


def difference_magnitude(flux,z):
    return 31.4-2.5*np.log10(flux)-Planck18.distmod(z).value+2.5*np.log10(1+z)


def main():
    a=pd.read_csv(HERE/'individual_apertures.csv').query('candidate==@KEY').copy()
    a['utc_day']=[Time(m,format='mjd').isot[:10] for m in a.mjd]
    a['detector']=a.exposure.str.split('_').str[-1]
    stats=a.groupby(['utc_day','band','detector']).agg(
        n=('fixed_njy','count'),median_fixed_njy=('fixed_njy','median'),
        min_fixed_njy=('fixed_njy','min'),max_fixed_njy=('fixed_njy','max'),
        median_centered_njy=('centered_njy','median'),min_coverage=('coverage','min'))
    stats.to_csv(HERE/'cosmos022_daily.csv',float_format='%.6g')
    a=pd.read_csv(HERE/'pixel_apertures.csv')
    r=a[(a.candidate==KEY)&(a.band=='f162m')&(a.radius_arcsec==.16)].iloc[0]
    tab=pd.read_csv(HERE/'host_neighbors.csv')
    h=tab[(tab.candidate==KEY)&(tab.catalog=='cosmos2020')].sort_values('sep_arcsec').iloc[0]
    d,e=r.delta_njy,r.empirical_err_njy
    result=dict(identifier=KEY,band='F162M',radius_arcsec=r.radius_arcsec,
        delta_aperture_njy=d,empirical_error_njy=e,host_id=int(h.id),host_sep_arcsec=h.sep_arcsec,
        lephare_redshift=h.lp_zpdf,lephare_68=[h.lp_zpdf_l68,h.lp_zpdf_u68],
        eazy_redshift=h.ez_z500,eazy_68=[h.ez_z160,h.ez_z840],
        spec_z=0.5067,spec_z_quality=3,spec_z_id=414071,spec_z_sep_arcsec=.0682,
        spec_z_source='COSMOS Spec-Z Compilation DR1.1, survey 103; host galaxy',
        observed_band_nominal_um=1.62,rest_band_nominal_um=1.62/(1+.5067),
        spectroscopic_host_absolute_difference_ab=difference_magnitude(d,.5067),
        flux_only_envelope_ab=[difference_magnitude(d+e,.5067),difference_magnitude(d-e,.5067)],
        photometric_host_absolute_difference_ab=difference_magnitude(d,h.lp_zpdf),
        photometric_host_envelope_ab=[difference_magnitude(d+e,h.lp_zpdf_u68),difference_magnitude(d-e,h.lp_zpdf_l68)],
        assumptions='Planck18; host association; fixed aperture; redshifted equivalent filters; no extinction, total-flux correction, or peak-time estimate.')
    (HERE/'cosmos022_host_magnitude.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
