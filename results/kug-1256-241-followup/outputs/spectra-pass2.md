# Spectral comparison and sensitivity

| Product | MJD | Rest-frame feature | Fixed-width flux ± formal error | 3σ sensitivity |
|---|---:|---|---|---:|
| spec-5991-56076-0484.fits | 56076.0 | Halpha | 2173.771 ± 51.86664 | 155.5999 |
| spec-5991-56076-0484.fits | 56076.0 | Hbeta | 233.7096 ± 48.73641 | 146.2092 |
| spec-5989-56312-0965.fits | 56312.0 | Halpha | 1962.439 ± 77.05221 | 231.1566 |
| spec-5989-56312-0965.fits | 56312.0 | Hbeta | 525.7214 ± 67.90993 | 203.7298 |
| coadd-0-496-20210506.fits-39628351622615360 | 59341.18021683 | Halpha | 495.1905 ± 59.93486 | 179.8046 |
| coadd-0-496-20210506.fits-39628351622615360 | 59341.18021683 | Hbeta | -226.0455 ± 62.00422 | 186.0127 |
| coadd-sv3-bright-28115.fits-39628351622615360 | 59341.18021683 | Halpha | 506.7686 ± 51.37946 | 154.1384 |
| coadd-sv3-bright-28115.fits-39628351622615360 | 59341.18021683 | Hbeta | -237.9817 ± 55.15819 | 165.4746 |
| coadd-0-504-20210529.fits-39628351622615360 | 59364.17969121 | Halpha | 521.8418 ± 56.35795 | 169.0738 |
| coadd-0-504-20210529.fits-39628351622615360 | 59364.17969121 | Hbeta | -249.414 ± 64.34264 | 193.0279 |

The common template has FWHM 5,000 km/s, centered on each vacuum Balmer transition. It uses each native observed flux density unit × observed Å. BOSS and DESI use 10⁻¹⁷ erg cm⁻² s⁻¹ Å⁻¹. LAMOST flux scales may be relative. Both BOSS spectra precede the 2014 slew. Coverage and dates for any later spectrum appear in the pixel index.

The continuum is linear, with a shared narrow velocity and width and a fixed doublet ratio for [O III] or [N II]. The error comes from weighted regression at the best narrow parameters, inflated by residual scatter. The 3σ sensitivity is the formal error times three. The upper95 column adds 1.645 errors to the positive part of the fitted amplitude. Both are conditional on this profile, the continuum model and the recorded aperture. The Hβ template can absorb Fe II, stellar absorption and mismodeled narrow lines. Use the free-profile fit and spectra as well when interpreting changes. The continuum-normalized plot permits visual line-shape comparison across different flux units.

Pixel arrays, masks expressed as zero inverse variance, flux units, input provenance and observation dates are deposited. Run `python spec_extra.py --replay` from the study folder with the pinned dependencies. This repeats the fixed-width measurements without querying archives. The original full FITS bytes are identified by SHA-256 in the receipts.
