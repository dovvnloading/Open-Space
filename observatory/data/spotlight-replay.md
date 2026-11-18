# Replaying the featured observations

This audit used the research code and input receipts at [revision 4a2f6d0](https://github.com/dovvnloading/Open-Space/tree/4a2f6d0936cd2e54335590bac41d1b62ae90aa53). Each catalog build also checks the twelve displayed readouts against the tables at its selected revision. A build-time table check and a fresh retrieval from an instrument archive answer different questions. The results of the fresh retrievals and the remaining uncertainties are recorded below.

## Mars

**Question.** Did MAVEN measure a sustained particle enhancement at Mars across multiple spacecraft orbits?

All four downloaded Parquet files matched the byte counts and SHA-256 digests in the survey receipt. Replaying the 1,189,114 SEP timestamps reproduced the complete summary and all 99 rows of the daily table. The August event selected 125 of 166 clock hours and spanned 46 orbits. Peak hourly medians divided by their respective survey references give 157,000 / 30.55 = 5,139.12 for ions and 616,500 / 25.2125 = 24,452.16 for electrons. The display rounds these ratios to the nearest ten. [Inputs, eligibility and full event set](../../results/maven-2025-sep-events/summary.json) · [Daily measurements](../../results/maven-2025-sep-events/daily.csv)

The original NASA PDS table and XML label were retrieved separately. The 45,698,463-byte table matched its archive MD5 and its deposited SHA-256. All 450 peak-hour timestamps matched between PDS and the mirror, and the median and 95th-percentile fractional differences were zero for each of the 16 flux and uncertainty columns. Both products encode the same instrument measurements. The two competing solar-eruption forecasts reproduced from the deposited forecast comparison. Identifying which eruption supplied the particles requires shock and spectral evidence. [Primary PDS receipt and limits](../../results/maven-2025-sep-events/primary_pds_audit.json) · [Forecast comparison](../../results/maven-2025-sep-events/forecast_alignment.json)

Reproduction uses the [seasonal survey](../../examples/maven_sep_final_season.py), [PDS decoder and comparison](../../examples/maven_pds_kp_primary_audit.py), and [forecast comparison](../../examples/maven_cme_forecast_alignment.py). Run their commands from the [investigation report](../../results/maven-2025-sep-events/README.md) with outputs directed to a new directory.

## Sun

**Question.** Which numbered solar region best matches the direction of the published STIX flare centroid?

The two original STIX FITS products, the full STIX imaging catalog and the NOAA region summary matched all four deposited SHA-256 digests. Rebuilding the FITS extract reproduced all 2,716 count-rate records and both sets of header values. The science product reached 286.35 counts/s in the 15–25 keV band at 18:27:38 spacecraft UTC. The Earth-equivalent time was 18:28:53 UTC. Quick-look and science products share the same photon observations. [Archive URLs and digests](../../results/2021-july-13-farside-eruption/inputs/sources.json) · [FITS extraction](../../results/2021-july-13-farside-eruption/extract_fits.py)

The published centroid and observer coordinates agreed with the frozen flare record. The geometry replay reproduced the deposited results for all seven region centers. The measured position angle is 104.99°. The two nearest directions have residuals of 2.92° for 12849 and 14.00° for 12847. Under the report’s illustrative coordinate perturbations, those ranges are 0.00–6.01° and 10.37–17.96°. Source height, coordinate drift and unnumbered regions leave the footpoint open. Independent reconstruction of the X-ray centroid remains pending. [Geometry replay](../../results/2021-july-13-farside-eruption/reproduce.py) · [All candidate regions](../../results/2021-july-13-farside-eruption/measurements.json)

## Quasar

**Question.** Do two resolved images of J0826+7002 show similar brightness changes at the same positive time delay across surveys?

The independent scalar scan in the observatory recovers 210 days from 418 B and 489 A ZTF g-band nights, and 205 days from 37 B and 39 A Gaia G-band nights. The image coordinates give 5.82949 arcseconds of separation. Running the original NumPy replay separately reproduced all 83 archived output lines, including zero-delay, reversed-order, alias, omitted-season and red-noise checks. The six follow-up tests passed. [Nightly measurements and original replay](../../results/j0826-lens-delay/README.md) · [Published output](../../results/j0826-lens-delay/output.txt)

Fresh Gaia A and B service replies had SHA-256 digests `ba02ee42a4ecae6b8cae12b8d322dea7f47e7421f49b49c0561c09acd66d59d3` and `3ef205a67365f3d69891f2bed2f71fce7c473a2334eda4455a880f869a24778a`. Their bytes differed from the archived response receipts. Reapplying the published cuts and aggregation nevertheless recovered all 219 deposited Gaia nightly rows across G, BP and RP exactly. The IRSA ZTF cone query timed out on two attempts. The successful ZTF replay used the deposited nightly product, SHA-256 `481c8362147bb00e76e4fda626d48e04d5611742f0dbadef18748f5d65f04975`. [Source receipts](../../results/j0826-lens-delay/sources.csv) · [Archive reduction](../../results/j0826-lens-delay/rebuild.py)

The locked-model future holdouts confirm the sensitivity to epoch. At the three published split dates the score difference between zero and +210 days is approximately +21.14, +5.50 and −0.26, with positive values favoring the delayed model. The latest block slightly favors zero over exactly +210. Changes to the intrinsic and microlensing models also shift the preferred lag. A narrower delay interval requires more independently resolved variability. [Holdout definitions, competing delays and follow-up evidence](../../results/j0826-lens-delay/followup/README.md)

## Dwarf

**Question.** Does motion of the radio peak identify WISE J0458+6434 as its emitter?

All four deposited VLASS FITS cutouts matched their SHA-256 receipts. Refitting the restoring beam recovered the four recorded forced flux densities and the three available centroids. The 2017 and 2025 radio centroids are separated by 3.28 arcseconds. The predicted infrared positions move 2.93 arcseconds over the same visits. The 2020 fitted flux density was −0.150 mJy/beam against a 0.148 mJy/beam local control scatter. All 8,260 translated tracks reproduced the archived rank counts. These are the deposited cutout pixels from the survey’s quick-look images. [Cutout origins and hashes](../../results/wise-j0458-radio/sources.json) · [Refit](../../results/wise-j0458-radio/replay.py) · [Control checks](../../results/wise-j0458-radio/check.py)

The native-pixel follow-up reproduced every row of the competing-template table. The combined endpoint contrast in favor of a moving beam over a beam fixed at the low-frequency catalog position ranks first among the target and 64 translated controls. Two variable stationary beams still fit the intermediate image, and the first survey epoch has imaging caveats. Reimaging the raw interferometer data and obtaining a new epoch with the target farther from both stationary positions remain open tests. [Native-pixel comparison](../../results/wise-j0458-heldout/native-pixel-tests.csv) · [Models and limitations](../../results/wise-j0458-heldout/README.md)
