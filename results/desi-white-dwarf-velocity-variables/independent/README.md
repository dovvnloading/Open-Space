# Independent checks of two DESI radial velocity candidates

The independent LAMOST spectrum extends the velocity record of WDJ023809.44+273329.79 back to **2014 November 2**. Neither star has a secure photometric orbital period in these data.

## A second spectrograph for J0238+2733

LAMOST observed J023809.47+273327.9 in plate GAC040N27M1, spectrograph 3, fiber 93, observation 257903093. Its position is within 1 arcsecond of the Gaia proper-motion track. DR12 classifies it as a white dwarf. The coadd header gives a heliocentric vacuum redshift of −0.00039822 ± 0.00015747, equivalent to **−119.4 ± 47.2 km/s**. The header records three integrations and EXPTIME = 1800 s. Their start minutes differ by 33 and 34 min. Temporal averaging matters for short orbits.

The LAMOST velocity differs from the DESI blue-arm −350.1 ± 20.9 km/s state by 230.7 km/s (4.47 combined formal standard errors) and lies within one combined standard error of −79.6 ± 22.2 km/s. Attributing the formal significance to orbital motion requires accounting for template offsets. DR5 fitted the **same exposure** with subtype A7IV and returned −60.8 ± 39.3 km/s. DR9 and DR12 returned the same WD redshift. In a separate exploratory fit to the DR12 wavelength pixels, changing the fitted Balmer lines and continuum windows shifted the joint result between approximately −115 and −7 km/s. The re-releases provide one observational epoch. The uncertainty in line modeling prevents a precise orbital constraint. [DR12 spectrum](https://www.lamost.org/dr12/v2.0/spectrum/fits/257903093) · [DR5 spectrum](https://dr5.lamost.org/spectrum/fits/257903093)

## Gaia and ZTF

The Gaia DR3 non-single-star flag is zero for both objects. Both report VarFlag = NOT_AVAILABLE, no epoch photometry, an empty Gaia radial velocity and zero percent multiple image peaks. RUWE is 1.049 for J0421−1325 and 0.926 for J0238+2733. These fields describe Gaia's catalogue coverage. The DR3 reference positions and proper motions in the attached VizieR receipts explain the offset from the J2000 designations. [Gaia DR3 catalogue](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=I/355/gaiadr3)

| Target | Clean g/r epochs | Best exploratory frequency (d⁻¹) | Period (min) | g semi-amplitude (mmag) | r semi-amplitude (mmag) | Per-target nominal search FAP |
|---|---:|---:|---:|---:|---:|---:|
| J0421−1325 | 819 / 913 | 54.323502 | 26.508 | 7.07 ± 1.62 | 7.22 ± 1.65 | 0.015 |
| J0238+2733 | 506 / 854 | 47.353499 | 30.409 | 5.36 ± 1.27 | 6.65 ± 1.43 | 0.012 |

The joined g/r ZTF DR24 search covers 0.2–144 cycles/day, after subtracting each object ID's median magnitude, removing duplicate exposures and filtering quality as implemented in verify.py. The uncertainty on amplitude is conditional on the selected frequency. The Baluev analytic false-alarm probabilities assume independent errors and apply to that target and that frequency search. We also inspected single-band, time-split and alias alternatives. The separate g and r searches select different global maxima. The frequencies in the table are unconfirmed. The 64 and 76 independently selected Pan-STARRS1 epochs give fixed-frequency semi-amplitudes of 6.2 ± 3.6 and 2.4 ± 4.9 mmag respectively. These uncertainties allow weak signals. ZTF source catalogues also omit undetected epochs. These searches do not exclude narrow deep eclipses. [ZTF DR24](https://irsa.ipac.caltech.edu/Missions/ztf.html) · [PS1 DR2](https://catalogs.mast.stsci.edu/panstarrs/)

## A bright trail over J0238+2733

At MJD 59147.267512 a ZTF i-band point appears 0.940 mag brighter than the object's median. In the native 50 × 50 pixel science cutout a bright diagonal trail crosses the measured white dwarf position. A principal-axis fit to 94 pixels more than 100 image counts above the median and outside a six-pixel radius around the source gives a trail spanning 62.6 pixels with a cross-track scatter of 0.53 pixels. Its axis passes within **0.22 pixels** of the source centroid. This exposure cannot measure an intrinsic flare. [ZTF exposure cutout](https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/2020/1025/267512/ztf_20201025267512_001650_zi_c02_o_q2_sciimg.fits?center=39.5394,27.557849&size=50arcsec&gzip=false)

Run `python verify.py` in this directory with NumPy, pandas, SciPy and Astropy installed. Add `--scan` to recompute the two full periodograms. The manifest records archived URLs and exact input byte hashes. Current delivery endpoints may serve changed releases. Omitted input files can be fetched using the script. A fresh VizieR query includes a new query timestamp and consequently has a different byte hash.

LAMOST is the Guo Shou Jing Telescope, a National Major Scientific Project built by the Chinese Academy of Sciences, funded by the National Development and Reform Commission and operated and managed by the National Astronomical Observatories, Chinese Academy of Sciences. ZTF and Pan-STARRS1 project and funding acknowledgments are provided by their respective release archives. This work has made use of data from the European Space Agency Gaia mission processed by the Gaia Data Processing and Analysis Consortium.
