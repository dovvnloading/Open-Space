## Public spectroscopy and astrometry searches on 2026-11-14

No additional spectrum or epoch radial velocity for TIC 96100077 was recovered in the searches below. The exact HTTP requests, returned counts and response digests are in [archive-searches.json](archive-searches.json). The digests identify the bytes received on the audit date. Archive responses can change.

| Source | Search | Result |
| --- | --- | --- |
| ExoFOP-TESS | Exact TIC 96100077 | Target identified, with zero spectroscopy entries, time series or uploaded files. |
| Gemini Observatory Archive and SALT Science Database | Pointing centers within 120 arcsec | Each returned zero files. |
| NOIRLab Astro Data Archive | File centers in a box extending approximately 2 arcmin on each side | Zero files. A separate CHIRON science-object query without the position restriction returned five records at its five-row limit. |
| Data Central spectral service and the GAVO FEROS and Flash/Heros services | SSA position 264.526816°, −49.349215° with the SIZE parameter set to 0.03333333° | Each returned a successful VOTable with zero spectra. |
| ESO raw science observations | Search GIRAFFE and UVES pointing centers in a rectangle containing all locations within 12.5 arcmin of the target, allowing for the full 25-arcmin FLAMES field. Repeat across all instruments. | Zero GIRAFFE or UVES exposures. The all-instrument search returned three GRIPS19 exposures labeled SKY MAP / SKY. A focused FEROS query also returned zero rows. [ESO FLAMES field specification](https://www.eso.org/sci/facilities/paranal/instruments/flames/overview.html) |
| DESI DR1, 6dFGS, SB9 and CRVAD-2 | DESI mean fiber coordinates in a box approximately 2 arcmin wide, and 30-arcsec cones on the three VizieR catalogs | Zero matching rows in each search. |
| Gaia DR4 prerelease | Check every source in the public 2026-06-26 epoch-astrometry example by sky position. | The table contains 1,008 entries from 12 distinct sources. The nearest source to the TIC position is 39.854° away. [ESA prerelease](https://www.cosmos.esa.int/web/gaia/dr4-prerelease) |

For queries based on pointing centers, a null result covers the stated search geometry. Multi-object spectra with a distant pointing center or spectra outside these public holdings can still exist. The raw AAT archive and CADC spatial queries timed out. The CDS associated-spectra endpoint returned HTTP 503. These three services remain unresolved.
