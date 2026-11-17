"""Fetch Gaia DR3 and every SDSS-V visit attached to the DR19 screened targets."""
from concurrent.futures import ThreadPoolExecutor
import csv
import io
import hashlib
from urllib.request import urlopen
from urllib.parse import urlencode
from fetch_survey import HERE,WORK,query
GAIA='https://gea.esac.esa.int/tap-server/tap/sync'
GAIA_COLS=('source_id,ra,dec,parallax,pmra,pmdec,ra_error,dec_error,parallax_error,pmra_error,pmdec_error,'
           'ra_dec_corr,ra_parallax_corr,ra_pmra_corr,ra_pmdec_corr,dec_parallax_corr,'
           'dec_pmra_corr,dec_pmdec_corr,parallax_pmra_corr,parallax_pmdec_corr,pmra_pmdec_corr,'
           'ruwe,visibility_periods_used,duplicated_source,astrometric_params_solved,'
           'phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,nu_eff_used_in_astrometry,pseudocolour,'
           'ecl_lat,ipd_frac_multi_peak,ipd_gof_harmonic_amplitude,astrometric_excess_noise_sig,'
           'radial_velocity,radial_velocity_error,rv_nb_transits,rv_chisq_pvalue,rv_amplitude_robust,'
           'teff_gspphot,logg_gspphot,ag_gspphot,phot_bp_rp_excess_factor,non_single_star')
SDSS_COLS = {'mwm_boss_allvisit':('sdss_id,gaia_dr3_source_id,spectrum_pk,mjd,fieldid,telescope,'
                                     'snr,zwarning_flags,xcsao_v_rad,xcsao_e_v_rad,xcsao_rxc,xcsao_teff'),
             'mwm_apogee_allvisit':('sdss_id,gaia_dr3_source_id,spectrum_pk,mjd,plate,fiber,telescope,'
                                      'snr,spectrum_flags,v_rad,e_v_rel,xcorr_v_rad,doppler_flags,'
                                      'doppler_rchi2,n_components,doppler_teff,on_target,valid')}
def fetch(job):
    table,ids=job
    if table=='gaia':
        adql=f'SELECT {GAIA_COLS} FROM gaiadr3.gaia_source WHERE source_id IN ({ids})'
        args={'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','QUERY':adql}
        with urlopen(GAIA+'?'+urlencode(args),timeout=180) as resp:
            raw=resp.read()
        if not raw.startswith(b'source_id,'):raise ValueError(raw[:300])
    else:
        sql=f"SELECT {SDSS_COLS[table]} FROM {table} WHERE sdss_id IN ({ids}) AND release='sdss5' ORDER BY sdss_id,mjd,spectrum_pk"
        raw=query(sql).split(b'\n',1)[1]
    return raw
if __name__=='__main__':
    rows=list(csv.DictReader((HERE/'prefilter_tail.csv').open()))
    for table in ['gaia',*SDSS_COLS]:
        key='gaia_dr3_source_id' if table=='gaia' else 'sdss_id'
        ids=sorted({r[key] for r in rows},key=int)
        jobs=[(table,','.join(ids[i:i+75])) for i in range(0,len(ids),75)]
        with ThreadPoolExecutor(5) as pool: blocks=list(pool.map(fetch,jobs))
        lines=[blocks[0].splitlines()[0]]
        for block in blocks: lines.extend(block.splitlines()[1:])
        p=HERE/f'{table}_tail.csv'
        body=b'\n'.join(lines)+b'\n'
        p.write_bytes(body)
        print(table,len(lines)-1,hashlib.sha256(body).hexdigest(),flush=True)
