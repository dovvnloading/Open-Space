"""Validate OGLE-SMC-ECL-06090 at the independently observed Gaia G epochs.

Run with Python, NumPy, SciPy and pyerfa. Paths are relative to this script.
The OGLE ephemeris is fixed to its OGLE-IV catalog entry, 06090.
"""
import csv
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode

import erfa
import numpy as np
from scipy.optimize import least_squares

HERE = Path(__file__).resolve().parent
PERIOD = 1.1696719
EPOCH = 2460000.0369
RA, DEC = np.radians([21.06275, -73.32594444444444])
LINE_OF_SIGHT = np.array([np.cos(RA)*np.cos(DEC), np.sin(RA)*np.cos(DEC), np.sin(DEC)])
SOURCE_ID = '4686449235706841984'
BASE = 'https://gea.esac.esa.int'
URLS = {
    'gaia-epochs.csv': BASE + '/data-server/data?' + urlencode({
        'RETRIEVAL_TYPE':'EPOCH_PHOTOMETRY', 'ID':'Gaia DR3 '+SOURCE_ID, 'VALID_DATA':'true', 'FORMAT':'CSV'}),
    'gaia-catalog-model.csv': BASE + '/tap-server/tap/sync?' + urlencode({
        'REQUEST':'doQuery', 'LANG':'ADQL', 'FORMAT':'csv',
        'QUERY':'SELECT * FROM gaiadr3.vari_eclipsing_binary WHERE source_id = '+SOURCE_ID}),
    'ogle4-I.dat': 'https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/phot/I/OGLE-SMC-ECL-06090.dat',
}


def to_hjd_utc(bjd_tcb):
    # ERFA epv00 gives position of Earth measured heliocentrically, barycentrically.
    # Their difference is the Sun's barycentric position. Subtract its projection.
    # Treat TDB as TT in the final clock conversion (sub-2-ms periodic error).
    a, b = erfa.tcbtdb(2451545., np.asarray(bjd_tcb)-2451545.)
    helio, bary = erfa.epv00(a, b)
    sun_projection = ((bary['p']-helio['p']) @ LINE_OF_SIGHT) * 499.004783836 / 86400
    tai1, tai2 = erfa.tttai(a,b)
    utc1, utc2 = erfa.taiutc(tai1,tai2)
    return utc1 + utc2 - sun_projection


def phase(hjd):
    return ((hjd-EPOCH)/PERIOD) % 1


def distance(ph, center):
    return (ph-center+.5) % 1-.5


def model(ph, params):
    base, d1, c1, w1, d2, c2, w2 = params
    return base + d1*np.exp(-.5*(distance(ph,c1)/w1)**2) + d2*np.exp(-.5*(distance(ph,c2)/w2)**2)


def fit_ogle(t, m, e):
    x = [np.median(m), .4,.01,.025,.4,.51,.025]
    lo = [17, .01,-.10,.008,.01,.40,.008]
    hi = [19,2,.10,.13,2,.60,.13]
    return least_squares(lambda z: (model(phase(t), z)-m)/np.hypot(e,.012), x,
                         bounds=(lo,hi), loss='soft_l1',f_scale=1.8).x


def run():
    data={}
    for name,url in URLS.items():
        blob=(HERE / name).read_bytes()
        data[name]={'url':url,'bytes':len(blob),'sha256':hashlib.sha256(blob).hexdigest()}
    rows=list(csv.DictReader((HERE/'gaia-epochs.csv').open()))
    catalog=next(csv.DictReader((HERE/'gaia-catalog-model.csv').open()))
    assert all(r['source_id']==SOURCE_ID for r in rows) and catalog['source_id']==SOURCE_ID
    clean=[r for r in rows if r['variability_flag_g_reject']=='false' and r['rejected_by_photometry']=='false'
           and r['g_transit_mag'] and r['g_transit_time'] and float(r['g_transit_flux_over_error'])>0]
    bjd=np.array([float(r['g_transit_time'])+2455197.5 for r in clean])
    hjd=to_hjd_utc(bjd)
    ph=phase(hjd)
    g=np.array([float(r['g_transit_mag']) for r in clean])
    t,m,e=np.loadtxt(HERE/'ogle4-I.dat').T
    ok=(t>=hjd.min())&(t<=hjd.max())&(e>0)&(e<=.10)
    fit=fit_ogle(t[ok],m[ok],e[ok])
    d1=np.abs(distance(ph,fit[2])); d2=np.abs(distance(ph,fit[5]))
    out=(d1>=.10)&(d2>=.10)
    base=np.median(g[out])
    window_summary={}
    for name,dist in [('primary',d1),('secondary',d2)]:
        inside=dist<.035
        window_summary[name]={'n':int(inside.sum()), 'median_magnitude':np.median(g[inside]),
                              'median_excess_vs_baseline_mag':np.median(g[inside])-base,
                              'observations':[{'hjd_utc':hjd[i],'ogle_phase':ph[i],
                                               'G':g[i],'transit_id':clean[i]['transit_id']} for i in np.flatnonzero(inside)]}
    f=float(catalog['frequency'])
    reference=float(catalog['reference_time']) + 2455197.5
    p1,p2=[float(catalog[f'geom_model_gaussian{i}_phase']) for i in (1,2)]
    e1,e2=[float(catalog[f'geom_model_gaussian{i}_phase_error']) for i in (1,2)]
    mid=to_hjd_utc(reference + np.array([p1,p2])/f)
    return {'inputs': data,'source_id':SOURCE_ID,'ogle_period_days':PERIOD,'ogle_primary_epoch_hjd':EPOCH,
            'clock':'Gaia times are BJD(TCB)-2455197.5. ERFA TCB to TDB, treat TDB as TT (within ~2 ms), convert TT to UTC, and subtract the projected barycentric position of the Sun. OGLE published HJD treated as UTC.',
            'gaia_rows_total':len(rows),'gaia_G_rejection_flag_count':sum(r['variability_flag_g_reject']=='true' for r in rows),
            'gaia_rows_used':len(clean),'gaia_hjd_range': [hjd.min(),hjd.max()],
            'gaia_bjd_tcb_minus_hjd_utc_range_seconds': [(86400*(bjd-hjd)).min(),(86400*(bjd-hjd)).max()],
            'contemporaneous_ogle_I_n':int(ok.sum()),'ogle_I_window_selection':'Error in (0, 0.10] mag and time between first and last usable Gaia transit. Fixed catalog ephemeris. Free constant plus two periodic Gaussians, robust soft-L1 fit with 0.012 mag noise floor and f_scale=1.8.',
            'ogle_I_fit_baseline_d1_c1_width1_d2_c2_width2':fit.tolist(),
            'ogle_I_separation':fit[5]-fit[2],
            'gaia_G_contrasts':{'response':'Within 0.035 cycle of respective OGLE-fitted center. Baseline at least 0.10 cycle from both centers. Window choice is exploratory. Gaia model uses the same Gaia photons as these counts.',
                               'baseline_n':int(out.sum()), 'baseline_magnitude':base,**window_summary},
            'gaia_published_model':{'model':catalog['model_type'],'period_days':1/f,
                'HJD_UTC_two_model_eclipse_centers_near_reference':mid.tolist(),
                'OGLE_phase_of_Gaia_centers':phase(mid).tolist(),'separation_in_Gaia_cycles':(p2-p1)%1,
                'separation_error_assuming_zero_covariance':float(np.hypot(e1,e2)),
                'G_depths_mag':[float(catalog['geom_model_gaussian1_depth']),float(catalog['geom_model_gaussian2_depth'])],
                'reduced_chi_squared':float(catalog['reduced_chi2'])},
            'limits':'Gaia measures an intermediate epoch. Its sparse core sampling and static catalog fit do not independently trace decades of drift. The quoted spacing error combines marginal errors without a covariance term.'}

if __name__=='__main__':
    (HERE/'validation.json').write_text(json.dumps(run(),indent=2)+'\n')
