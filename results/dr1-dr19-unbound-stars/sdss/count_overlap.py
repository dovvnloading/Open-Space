"""Count overlap between the two eligible coadd tables using the exact base masks."""
from fetch_survey import HERE,query
q="""SELECT COUNT_BIG(*) AS n_both FROM mwm_boss_allstar b JOIN mwm_apogee_allstar a ON b.sdss_id=a.sdss_id
WHERE b.release='sdss5' AND b.gaia_dr3_source_id>0 AND b.snr>=15 AND b.e_v_rad>0 AND b.e_v_rad<=25 AND ABS(b.v_rad)<5000
AND a.release='sdss5' AND a.gaia_dr3_source_id>0 AND a.snr>=15 AND a.n_good_rvs>=1 AND ABS(a.v_rad)<5000"""
if __name__=='__main__':
    p=HERE/'audit'/'overlap.csv'
    p.parent.mkdir(exist_ok=True)
    p.write_bytes(query(q).split(b'\n',1)[1])
