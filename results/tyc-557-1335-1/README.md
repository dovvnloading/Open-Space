# TYC 557-1335-1: occultation double conflicts with a shared 59-day orbit

On 11 September 2022, (59592) 1999 JW58 briefly blocked two sources toward **TYC 557-1335-1**, also known as **UCAC4 479-129046**. Gaia spectroscopy independently gives this line of sight a 59.392 ± 0.066 day orbit. At the measured distance and eccentricity, putting the two occulted sources on that orbit requires a minimum total mass of **3,080 solar masses**. The most permissive test, allowing any bound eccentricity and shifting all input measurements against the inference, still requires **more than 109 solar masses**. **The third star remains unconfirmed.** The weaker second occultation dip was observed at **one station only** and appears in just **two full exposures** in the processed display video. The Gaia single-lined spectroscopic orbit relies on only **11 usable radial velocities** and its period has no independent confirmation. If both the secondary dip and the spectroscopic orbit are genuine, a third stellar component is one explanation; an unrelated optical companion, asteroid satellite, spurious secondary dip or an aliased radial-velocity period remain possible. A physical outer orbit and the identity of the spectroscopic pair have not been established.

## Identification and selection

Sposetti's observing report, figure 4 of the 2023 IOTA meeting report, and the European occultation archive identify **UCAC4 479-129046**. Its UCAC fit magnitude is 11.367 and its coordinates match Gaia DR3 **2697389613944690688**. The current IOTA XML also uses this Gaia ID. The observers discovered the double in 2022. The new inference is that its two resolved sources cannot account for the independently measured radial-velocity orbit as a single ordinary stellar pair. [IOTA meeting report, figure 4](https://iota-es.de/JOA/JOA2023_4.pdf#page=32) · [European event archive](https://euraster.ericfrappa.com/results/2022/index.html) · [UCAC4 catalogue row](inputs/ucac4.csv)

The 7 November 2026 Occult archive member contains 17,326 event records. All 357 records carrying a DoubleStar element were retained, including repeat records and entries without a positive two-dimensional separation. They give 355 distinct Gaia IDs. An exact-ID join to Gaia DR3's two-body-orbit table returns 18 solutions for 17 sources. TYC 557-1335-1 is the highest-significance single-lined spectroscopic solution in this matched sample. Its four geometries all have positive separations and reported errors. We selected it after inspection. No population incidence or search-wide probability is inferred. [All 357 rows](inputs/screen.csv) · [Archive provenance](inputs/sources.json)

## Independent reduction of the recorded dips

Sposetti recorded from Gnosca, Switzerland with a 42-cm telescope at 20 ms integration. The original analysis gives durations of 0.242 ± 0.020 s and 0.051 ± 0.020 s. His observer report gives approximate drops of 1.2 and 0.6 mag. The first dip was followed by a return to baseline before the shorter dip. Different depths favor occultations of stars with different fluxes over a single point source successively hidden by two opaque bodies, although a grazing or partially sampled second event remains possible. [Original report, light curve and video](https://www.occultations.ch/occultations_2022.html)

We independently extracted aperture flux from the first 120 display frames. The GPS stamp advances monotonically throughout the selected segment. Later frames contain slow-motion repeats. A three-pixel aperture in the resized 480 × 300 video gives the following mean fluxes relative to off-event baseline. The references are the three brightest isolated stars above the target in the stacked display field.

| Source | Twelve full primary frames | Two full secondary frames |
|---|---:|---:|
| Target | 0.191 | 0.413 |
| Reference 1 | 1.049 | 1.113 |
| Reference 2 | 1.057 | 0.914 |
| Reference 3 | 1.066 | 0.969 |

No off-event consecutive pair among 95 overlapping target windows is as faint as the secondary at any of radii 2, 3, or 4 pixels. The measured secondary flux spans 0.413–0.483 over those apertures. The independent reduction uses the same observation, a lossy processed display video and only two full exposures for the secondary. These counts support the presence of a second localized dip. They are unsuitable for calibrated component photometry or a discovery false-alarm probability. [Frame photometry](inputs/video_photometry.csv) · [Video reduction](reduce_video.py)

## Which stars can share the 59-day orbit?

The Gaia SB1 has 11 usable radial-velocity observations, semi-amplitude **22.386 ± 0.247 km/s**, eccentricity 0.0836 and parallax **1.12835 ± 0.03168 mas** from the Gaia source table. Its spectroscopic mass function is 0.0683 solar masses. Bashi et al.'s published equation 7 and table 1 give score 0.761, above their clean-sample cut of 0.587. The score is a classifier and a short series of 11 velocities admits aliases. The period lies within 3.6 days of the Gaia scan precession period. Holl et al. give a geometrical RVS wavelength shift of 0.146 km/s per milliarcsecond along scan. At the widest measured separation, the maximum full intercomponent shift is **0.83 km/s**, 27 times smaller than the observed semi-amplitude. This tests simple scan-angle displacement at the occultation epoch. Gaia observed earlier, and unrelated contaminating spectra, blending of physically offset lines and RV period aliases remain possible. A 60-arcsecond Gaia cone contains five other catalogued objects. The closest is 14.40 arcseconds away and 8.05 G magnitudes fainter. One Gaia transit was deblended. A neighbour unresolved on milliarcsecond scales remains possible. No independent spectrum was obtained for this work. [Gaia inputs](inputs/gaia.csv) · [Neighbourhood](inputs/field_60arcsec.csv) · [Bashi et al.](https://doi.org/10.1093/mnras/stac2928) · [Gaia scan-angle systematics](https://doi.org/10.1051/0004-6361/202245353)

For projected separation s = ρ/ϖ in au and period P in years, an orbit with eccentricity e must obey s ≤ a(1+e). Kepler's law then gives M ≥ [s/(1+e)]³/P². Each row below represents a fit to the same positive chord.

| Position angle | Separation | Projected distance | Minimum total mass using measured e |
|---:|---:|---:|---:|
| 239.2° | 5.3 ± 0.6 mas | 4.70 au | 3,080 M☉ |
| 212.3° | 5.7 ± 0.5 mas | 5.05 au | 3,831 M☉ |
| 229.0° | 5.3 ± 0.6 mas | 4.70 au | 3,080 M☉ |
| 256.0° | 5.7 ± 0.5 mas | 5.05 au | 3,831 M☉ |

For the most adverse calculation we reduced each separation by three quoted errors, increased the period and parallax by three quoted errors, and used r < 2a for **any** bound eccentricity. The smallest resulting floor is 109.97 solar masses. The formal uncertainty stress test shifts the three inputs simultaneously. Inverting the same relation, even a 5-solar-mass system at the smallest nominal separation requires P > 1.61 years at arbitrary eccentricity. [Geometry](inputs/geometry.csv) · [Calculated results](analysis.json)

Searches within 30 arcseconds in the 2023 update of Tokovinin's Multiple Star Catalog and the VizieR Washington Double Star table returned no rows. SIMBAD classifies the source as SB*, and its aliases include TIC 319470579 and 2MASS J21534385+0542140. Searches by UCAC4 identifier, Gaia ID and principal aliases recovered the previously reported visual duplicity but no earlier connection to this SB1 orbit. These searches cannot establish absolute priority. This single positive station also leaves an asteroid satellite, an unrelated optical companion and a spurious secondary dip as alternatives. Differential measurements over more than one epoch and phase-resolved spectroscopy would distinguish the physically bound hierarchy from them. [Catalogue query URLs](inputs/catalog_checks.csv) · [SIMBAD record](https://simbad.cfa.harvard.edu/simbad/sim-id?Ident=TYC+557-1335-1)

Run `python analyze.py` for every dynamical bound, aperture comparison and classifier value. To repeat the photometry, download the AVI at the URL in `inputs/sources.json` and run `python reduce_video.py PATH_TO_AVI` with NumPy and FFmpeg installed. The reducer verifies its SHA-256 before decoding.
