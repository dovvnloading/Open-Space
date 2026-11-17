# SDSS-V DR20 high-speed white-dwarf classifications checked against spectra

All five objects in the DR20 tail have BOSS spectra consistent with ordinary F or K stars. The two marginal `DZ:` labels correspond to luminous, reddened objects with K4 III spectral templates and low surface gravities from BOSS Net. Neither has demonstrated LP 40-365-like chemistry. The five objects come from a white-dwarf classification pipeline, with no requirement that they were originally targeted as white dwarfs. Two were flux standards. [SDSS Snow White pipeline and limitations](https://www.sdss.org/dr20/mwm/astra/pipelines-in-astra/snowwhite/) · [SDSS advice on stellar parameters](https://www.sdss.org/dr20/mwm/data/parameters/)

| Gaia DR3 source | Snow White | Corrected tangential speed (km/s) | Nominal M_G | BOSS optical template | BOSS Net T_eff (K) / log g | Balmer-core line shift (km/s) |
|---|---|---:|---:|---|---:|---:|
| 1542290255978428800 | DA: | 681 | 4.69 | F2 III | 6093 / 4.52 | −215 ± 29 |
| 6504108468006561408 | DA/DAZ | 602 | 4.17 | F3/F5 V | 6169 / 4.56 | +75 ± 6, +66 ± 6 |
| 6706016839105528064 | DA: | 589 | 0.87 | F2 III | 6437 / 3.41 | −83 ± 6 |
| 6060880483145981568 | DZ: | 588 | −0.03 | K4 III | 4839 / 2.24 | +318 ± 9 |
| 4106812887663593088 | DZ: | 571 | 0.72 | K4 III | 4917 / 2.05 | +4 ± 15 |

Nominal absolute magnitudes use the corrected inverse parallax and the observed G magnitude before dereddening. The two K giants have corrected parallax S/N of 3.48 and 6.02. Their model-dependent Gaia GSP-Phot G-band extinctions are 1.95 and 2.44 mag, which would make their intrinsic magnitudes still brighter. Their Gaia GSP-Phot temperatures are about 6100 K. The differences from the BOSS Net temperatures illustrate the dependence of the BP/RP temperature estimates on reddening and the atmospheric model. Both have Galactic latitude within 6° of the plane. BOSS observed both in its metal-poor halo programme (`mwm_halo_mp_xp_boss_single`). BOSS Net estimates [Fe/H] of −1.07 and −1.64 respectively. The reported metallicities use a normal-star model.

The nominal DZ objects show absorption in the Hβ window of **0.98 ± 0.29** and **1.10 ± 0.29 Å** respectively and Ca II K window absorption of 6.70 ± 1.53 and 7.77 ± 0.52 Å. The 606088… spectrum also has Hα at 1.29 ± 0.24 Å. Their Mg I 4703 measurements are 0.15 ± 0.24 and 0.05 ± 0.23 Å. Interstellar Ca/Na absorption and blends complicate the elemental interpretation at these low latitudes. Unusable line windows remain blank. The luminosities, giant templates, low fitted gravities and ordinary-star classification jointly account for the spectra. A catalogue flag or a single metal line would not determine an elemental enrichment.

SIMBAD exact-ID and 10 arcsec cone searches identify 1542290255978428800 with published F5 V star SDSS J122601.09+462800.7 (catalogued RV −219 km/s) and 6706016839105528064 as an RR Lyrae variable. No positional association was returned for the other three. Absence from SIMBAD is insufficient evidence of novelty.

## Reproduce

Run `python check_tail.py` in an environment with NumPy, SciPy and Astropy. The script uses `../sdss/tail.csv`, the local `../spectra/fit_control.py`, Gaia DR3, SIMBAD, the public Valis API and BOSS reduction v6_2_1 daily coadds. The metadata verifies the exact Gaia DR3 identifier for all five sources. Six distinct nightly spectra were downloaded, and every input URL and SHA-256 is recorded in `receipts.json`. Overlapping epoch coadds are omitted. The BOSS redshift warning and BOSS Net result flags are zero for all six.

The independent fitter operates on valid, unmasked pixels. It fits the cores of available Hγ, Hβ and Hα with one Doppler shift, local linear continua, and independent widths and amplitudes, seeded by the survey redshift. This is a consistency check on a specified set of lines, not a blind radial-velocity discovery. Formal errors omit calibration, atmospheric profile mismatches and stellar variability. Pseudo-equivalent widths integrate the measured per-pixel flux deficit around each Doppler-shifted laboratory wavelength with local continuum uncertainty. The continuum and the preselected line windows can contain other transitions. No abundance ratios or chemical upper limits follow directly from them. [SDSS-V DR20](https://sdss.org/dr20/data_access/) · [Valis](https://sdss.org/dr20/data_access/zora_and_valis/) · [SIMBAD](https://simbad.cds.unistra.fr/simbad/)
