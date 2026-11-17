# Optical and catalog follow-up of the published 100–600 K subset

The temperature selection includes all 16 stars in Contardo & Hogg's table with published single-blackbody temperatures in the study range. Their fitted infrared luminosities and temperatures were inferred from unresolved photometry. [Contardo & Hogg 2024](https://doi.org/10.3847/1538-3881/ad6b90)

Run `python followup.py` with NumPy. The script compares the epoch-propagated Gaia J2000 positions with SIMBAD objects within 15 arcsec, preserving classifications and angular separations. An object in the cone may be a neighbor. A missing catalog identification provides no constraint on undiscovered dust, galaxies or companions.

ZTF DR24 r-band measurements use the same time-dependent astrometric and photometric cuts as [the Hephaistos test](../optical/ztf.py). Each retained object ID has at least five accepted measurements. Nightly median magnitudes are compared to that object's median night, flagging drops of at least 0.2 mag. A flagged night requires image and comparison-star vetting before interpretation. Sparse epochs cannot exclude short, shallow or constant obscuration. The final table separates archive errors, empty results and qualifying curves. Response hashes, query URLs, retained measurements by object and the compressed archive responses accompany the summary.

The SIMBAD neighborhood search also includes Hephaistos A–J. Each sightline uses its own indexed 15 arcsec cone and records any failed mirror requests. The displayed coordinate separation refers to Gaia propagated to epoch J2000 and the published SIMBAD coordinates. A zero-row success and a query error have separate receipt states.
