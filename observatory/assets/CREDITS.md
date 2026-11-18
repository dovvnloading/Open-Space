# Observatory image plates

The [plate manifest](plates.json) records the image file, attribution, context and SHA-256 digest for each cover. The three JPEGs preserve the bytes returned by their archives. The dossiers separately link to the observations used for their measurements.

## Mars

NASA/JPL-Caltech assembled 102 Viking Orbiter images into a perspective mosaic of the Valles Marineris hemisphere. The site displays NASA's 1,024-pixel rendition. MAVEN measured the ions and electrons in the August 2025 dossier from orbit. [NASA image and credit](https://science.nasa.gov/resource/mars-planet-globe/)

## The Sun

The 1,024-pixel browse frame credits STEREO/SECCHI/NRL. The image carries the STEREO Ahead COR2 instrument name and the timestamp 2021-07-13 19:08:45 UTC. The occulting disk blocks the central Sun. DONKI lists coronal mass ejections starting at 16:12 and 18:53 UTC and describes overlap in COR2 imagery. This frame supplies coronal context. The direction used to compare possible source regions comes from Solar Orbiter/STIX hard X-rays. [COR2 daily archive](https://stereo-ssc.nascom.nasa.gov/browse/2021/07/13/ahead/cor2/1024/) · [Deposited DONKI entries](../../results/2021-july-13-farside-eruption/inputs/donki_extract.json)

## J0826+7002

The Pan-STARRS1 Surveys and STScI/MAST supply the i, r and g stacks mapped to red, green and blue. The request takes 160 native pixels per side at approximately 0.25 arcseconds per pixel and renders a 512-pixel JPEG. The requested center (126.5341175°, +70.0456447222222°) is the mean of the two catalog positions in the deposited target table. Their great-circle separation is 5.8295 arcseconds. A matching r-band FITS cutout places A at (76.83, 68.39) and B at (82.92, 90.91) on the zero-indexed native grid. Its center is 0.101 arcseconds from the requested midpoint. The survey image shows the field. The timing comparison uses ZTF and Gaia light curves. [Target positions](../../results/j0826-lens-delay/targets.csv) · [Matching r-band FITS cutout with celestial coordinates](https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?ra=126.5341175&dec=70.0456447222222&size=160&format=fits&red=/rings.v3.skycell/2550/058/rings.v3.skycell.2550.058.stk.r.unconv.fits) · [PS1 archive and survey acknowledgements](https://outerspace.stsci.edu/spaces/PANSTARRS/pages/298812201/Pan-STARRS1%2Bdata%2Barchive%2Bhome%2Bpage) · [Pan-STARRS1 Surveys, Chambers et al.](https://arxiv.org/abs/1612.05560)

## WISE J0458+6434

Open-Space created this brown-dwarf illustration with OpenAI ImageGen. The cloud bands, colors and globe are artistic choices. The investigation follows an unresolved radio source near the binary. Its evidence figure shows the archived radio images. [Original artwork record](https://github.com/dovvnloading/Open-Space/blob/4a2f6d0936cd2e54335590bac41d1b62ae90aa53/observatory/assets/CREDITS.md) · [Radio investigation](../../results/wise-j0458-radio/README.md)
