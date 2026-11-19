# Referee checklist

The machine [independent replay](referee.py) must verify freeze hash, field and image identities, non-overlapping detector layers and dates, good-pixel areas, trigger counts, and result status. A detector-layer point-coverage test does not prove every pixel's layer provenance.

VIS/NISP mosaics are not independent calibrated exposure-level frames. Gaussian PSF matching and linear WCS interpolation can generate artefacts. Source masks and flags must be verified for the DR1 products. Matched-star flux scaling, nearby blends, streaks, cosmic rays, persistence and negative residuals require independent review.

The injection test is conditional on a held-fixed filter and photometric scale, and does not calibrate historical false-alarm rates, a survey volume, high-z/SLSN/TDE light curves, or control time. Even a zero-candidate search here is not a validated volumetric or temporal rate limit without those calibrations.

Gaia, MPC/Horizons, TNS, host redshift and literature checks must be completed and independently repeated before any candidate is classed as a supernova/TDE or labelled [REMARKABLE]. An unqueried or failed endpoint is never a clean association test.
