"""Derive conditional Einstein mass and light for the Q1 systems.

python infer.py MODEL_TABLE_DIRECTORY PHZ_CSV APERTURES_CSV [--out PATH]
See README.md for geometry, passbands, selection, and systematic errors.
"""
import argparse
import ast
import csv
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.integrate import quad

HERE = Path(__file__).resolve().parent
BANDS = ('vis', 'y', 'j', 'h')
BAND_KEY = dict(vis='VIS', y='NIR_Y', j='NIR_J', h='NIR_H')
PC_M = 3.085677581491367e16
C_KMS = 299792.458
ARCSEC_RAD = np.pi / 648000
# M = c^2 theta^2 D_l/(4 G beta), using CODATA G and IAU nominal solar GM.
M_SCALE = (C_KMS * 1e3)**2 * (1e6 * PC_M) * ARCSEC_RAD**2 / (4 * 1.3271244e20)


def read(path):
    with path.open() as f:
        return list(csv.DictReader(f))


def write(path, rows):
    if not rows:
        return
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow({k: ((f'{v:.9f}' if k in ('ra','dec','ra_deg','dec_deg','ra_center','dec_center') else f'{v:.7g}') if isinstance(v, (float, np.floating)) and np.isfinite(v)
                            else '' if isinstance(v, (float, np.floating)) else v) for k, v in row.items()})


def sep(ra, dec, ra2, dec2):
    dra, ddec = np.deg2rad(ra - ra2), np.deg2rad(dec - dec2)
    h = np.sin(ddec/2)**2 + np.cos(np.deg2rad(dec)) * np.cos(np.deg2rad(dec2)) * np.sin(dra/2)**2
    return np.rad2deg(2 * np.arcsin(np.sqrt(np.clip(h, 0, 1)))) * 3600


@lru_cache(maxsize=8192)
def chi(z):
    return C_KMS / 70 * quad(lambda x: (0.3 * (1+x)**3 + 0.7)**-0.5, 0, z)[0]


def mass_floor(theta, z):
    return M_SCALE * chi(z)/(1+z) * theta**2


def luminosity(mag, z, sun_mags, band):
    distance_pc = 1e6 * chi(z) * (1+z)
    solar = np.interp(z, sun_mags['z'], sun_mags[band], left=np.nan, right=np.nan)
    abs_ab = mag - 5*np.log10(distance_pc/10) + 2.5*np.log10(1+z)
    return 10**(0.4 * (solar - abs_ab))


def phz_bounds(row):
    if 'phz_95_int' in row:
        return ast.literal_eval(row['phz_95_int'])
    return float(row['phz_90_int1']), float(row['phz_90_int2'])


def run(data, phz_path, aperture_path, out):
    out.mkdir(parents=True, exist_ok=True)
    rojas = read(HERE / 'rojas_spectroscopy.csv')
    desi = read(HERE / 'desi_matches.csv') if (HERE / 'desi_matches.csv').exists() else []
    cat = {r['id_str']: r for r in read(data / 'q1_discovery_engine_lens_catalog.csv')}
    phz = {r['object_id']: r for r in read(phz_path)}
    flux = {r['id_str']: r for r in read(data / 'modeling_mge_magnitude.csv')}
    models = read(data / 'modeling_lens_mass.csv')
    ap = {(r['id_str'], r['band']): r for r in read(aperture_path) if r['aperture'] == 'einstein'}
    sun = np.genfromtxt(HERE / 'sun_shifted_band_ab.csv', delimiter=',', names=True)
    systems, bands, pairs = [], [], []
    for m in models:
        key = m['id_str']
        obj = key.split('_')[-1].replace('NEG','-')
        coords = cat.get(key)
        if not coords:
            continue  # orphan in the released model catalogue
        ra, dec = float(coords['right_ascension']), float(coords['declination'])
        matched = [r for r in rojas if sep(ra, dec, float(r['ra_deg']), float(r['dec_deg'])) < 1.]
        assert len(matched) <= 1
        p = phz.get(obj)
        zp = float(p['phz_median']) if p and p['phz_median'] else np.nan
        priors = [r for r in desi if r['id_str'] == key and r['adoptable'] == '1']
        zl, zsource, lens_type, source_type = np.nan, np.nan, '', ''
        if matched:
            r = matched[0]
            zl, lens_type = float(r['zl_spec']), 'Rojas_spec'
            if r['zs_spec']:
                zsource, source_type = float(r['zs_spec']), r['zs_quality']
        elif len(priors) == 1 and np.isfinite(zp) and abs(float(priors[0]['z'])-zp)/(1+zp) < 0.1:
            zl, lens_type = float(priors[0]['z']), 'DESI_primary_PHZ_consistent'
        elif p and p['phz_classification'] == '2' and zp > 0:
            zl, lens_type = zp, 'merged_PHZ'
        theta = float(m['einstein_radius_effective_median_pdf'] or 'nan')
        mf = mass_floor(theta, zl) if zl > 0 else np.nan
        beta = 1-chi(zl)/chi(zsource) if zsource > zl > 0 else np.nan
        lo, hi = phz_bounds(p) if p and p['phz_median'] else (np.nan, np.nan)
        systems.append(dict(id_str=key,ra_deg=ra,dec_deg=dec,grade=coords['grade'],
                            zl=zl,zl_basis=lens_type,zl_merged_PHZ=zp,PHZ_lo=lo,PHZ_hi=hi,
                            zs=zsource,zs_basis=source_type,beta=beta,
                            theta_eff_arcsec=theta,M_beta1_Msun=mf,M_spectra_Msun=mf/beta))
        if not (zl > 0 and theta > 0 and np.isfinite(theta)):
            continue
        for band in BANDS:
            obs = ap.get((key, BAND_KEY[band]))
            prefix = 'vis' if band == 'vis' else 'nir_' + band
            mag_model = float(flux[key][prefix + '_lens_magnitude_ab_median_pdf'] or 'nan')
            ltot = luminosity(mag_model, zl, sun, band)
            snr, lap, lum_limit, apmag = np.nan, np.nan, np.nan, np.nan
            if obs and obs['source_subtracted_snr'] and obs['noise_only_5sigma_mag']:
                counts = float(obs['source_subtracted_counts'])
                snr = float(obs['source_subtracted_snr'])
                sigma = 10**(.4*(float(obs['zp_ab'])-float(obs['noise_only_5sigma_mag'])))/5
                flux_upper = max(0.,counts) + 5*sigma
                lum_limit = luminosity(float(obs['zp_ab']) - 2.5*np.log10(flux_upper), zl, sun, band)
                if snr >= 5 and counts > 0:
                    apmag = float(obs['source_subtracted_mag'])
                    lap = luminosity(apmag, zl, sun, band)
            bands.append(dict(id_str=key,band=band,zl=zl,zl_basis=lens_type,zs=zsource,theta_eff_arcsec=theta,
                              mag_MGE_total=mag_model,L_total_Lsun=ltot,mag_E=apmag,SNR_E=snr,L_E_Lsun=lap,
                              L_E_upper5_Lsun=lum_limit,M_beta1_Msun=mf,ML_total_beta1=mf/ltot,
                              ML_E_beta1=mf/lap,ML_E_conservative5_beta1=mf/lum_limit,ML_E_spec=mf/lap/beta))
            if np.isfinite(beta) and np.isfinite(lap):
                pairs.append(dict(id_str=key,band=band,zl=zl,zs=zsource,zs_basis=source_type,
                                  M_E_Msun=mf/beta,L_E_Lsun=lap,ML_E=mf/beta/lap,
                                  ML_at_radius_low_flux_high=mass_floor(float(m['einstein_radius_effective_lower_1_sigma']),zl)
                                  /beta/lap/(1+1/snr),ML_at_radius_high_flux_low=mass_floor(float(
                                      m['einstein_radius_effective_upper_1_sigma']),zl)/beta/lap/(1-1/snr)))
    write(out/'redshift_geometry.csv', systems)
    write(out/'mass_light.csv', bands)
    write(out/'spectroscopic_controls.csv', pairs)
    bright = []
    for row in read(HERE/'ecker_bright_q1.csv'):
        if not row['theta_E'] or (not row['zl_spec'] and (not row['zl_phot'] or '(' in row['zl_phot'])):
            continue
        zl = float(row['zl_spec'].split('^')[0] or row['zl_phot'])
        theta, mag = float(row['theta_E']), float(row['mag_IE'])
        mass, lum = mass_floor(theta,zl), luminosity(mag,zl,sun,'vis')
        bright.append(dict(name=row['name'],ra_deg=row['ra_deg'],dec_deg=row['dec_deg'],
                           zl=zl,zl_basis='spec' if row['zl_spec'] else 'PHZ',theta_arcsec=theta,
                           mag_VIS_catalog=mag,M_beta1_Msun=mass,L_catalog_Lsun=lum,
                           ML_catalog_beta1=mass/lum))
    write(out/'ecker_bright_geometry.csv', bright)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('models', type=Path)
    parser.add_argument('phz', type=Path)
    parser.add_argument('apertures', type=Path)
    parser.add_argument('--out', type=Path, default=HERE)
    args = parser.parse_args()
    run(args.models, args.phz, args.apertures, args.out)
