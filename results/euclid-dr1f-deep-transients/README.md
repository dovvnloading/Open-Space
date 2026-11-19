# Euclid DR1-Foundation deep-field transient differences

**Claim first: Awaiting the fresh DR1-Foundation pixel scan. No validated transient or calibrated rate limit is asserted.** This is a small, frozen, three-field experiment, not the full 53 deg² deep survey.

The [frozen protocol](freeze.json) predeclares one **2 × 2 arcmin** region each in North, South and Fornax (total nominal **0.003333 deg²**) and a deterministic choice of the first three lexical MER science observations per VIS, NISP Y, J, H band. Field/visit lists and product URLs with acquisition SHA256 receipts are frozen in [output/frozen_manifest.json](output/frozen_manifest.json) before looking at pixels. The North patch overlaps the prior [development search](../euclid-foundation-eligibility/README.md), whose two raw excursions were a streak and a blend, so North is not a blinded discovery sample.

On a networked GitHub Actions runner, [the pipeline](run.py) uses the IRSA deep SIA collection for science/RMS/flag/detector-layering object identities, then the public S3 bucket with FITS **byte-range sections** for science pixels. It never assumes that unsupported IBE cutout query parameters succeeded. It requires disjoint detector-layer exposures at the center, time separation over one day, joint good-pixel mask and stellar PSF/relative-flux calibration before differencing. A Gaussian PSF match to the worse epoch precedes both signed peak screens. An empirical positive/negative control-noise scale, PSF-shape/streak cuts and synthetic source injection at 7/10/15 control sigma measure a *conditional difference-map* efficiency, not exposure-level or physical-light-curve completeness.

Raw triggers remain unvalidated until detector persistence, cosmic-ray and ghost checks, moving-body ephemerides (JPL Horizons/MPC), Gaia variables, host and photometric redshift (Euclid/DESI), colour and light curve, TNS and literature, and independent revisits have been examined. Missing service results are not negative evidence. No temporal/volumetric transient rate limit follows without a calibrated light-curve control time and completeness. Never sum overlapping band-pair areas.

[Referee requirements](REFEREE.md) · [Machine independent replay](referee.py) · [Pipeline tests](tests.py).

This work has made use of the Data Release (DR1) data from the Euclid mission of the European Space Agency (ESA). See [IRSA DR1 overview](https://irsa.ipac.caltech.edu/data/Euclid/docs/overview_dr1.html) and [MER image products](https://irsa.ipac.caltech.edu/euclid-archive-documentation/data-products/).
