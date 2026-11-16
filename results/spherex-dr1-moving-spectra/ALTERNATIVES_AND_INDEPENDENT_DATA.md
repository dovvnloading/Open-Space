# Fixed physical alternatives and independent-data search

## Alternative templates
Same official PSG/CEM configuration provenance and generic geometry as `PSG_PHYSICAL_SHAPE_REPORT.md`; fixed100 K, pure GSFC[CO] or GSFC[H2O]. Molecular lifetimes1300000s and77000s are taken directly from the corresponding entries in the original official67P template. Wavelength coverage extended to4.10–5.0µm. Successful returned headers explicitly identify the requested molecule. Only final molecular column is used.

`fit_psg_alternatives.py` convolves each with the same exact shifted response and fits nonnegative amplitude plus signed constant/linear continuum per epoch. Input is the frozen `psg_fit_input_snapshot.json`. Full output is `psg_alternatives_fit.json`.

At aperture2.5 with formal errors and constant continuum:

|Fixed model|December training chi-square|January2 chi-square|
|---|---:|---:|
|CO100K|1014.95|1707.64|
|H2O100K|1924.09|615.40|
|CO2 100K(previous test)|6.47|4.879|

CO has zero output in4.10–4.45µm in this particular model; its major modeled band is near4.7µm. Integrating the instrument's small offband tails can give nonzero template entries, so any large fitted CO normalization or apparent improvement from such tails is **not a meaningful CO production constraint**. The model's4.10–4.45 flux fraction relative to4.10–5.0 is0 for CO and0.03590 for H2O. This excludes these fixed pure alternatives as explanations of the observed feature, not their presence elsewhere in the comet. It does not rule out all temperatures, species, mixtures, contaminants, or instrumental artifacts.

Linear-continuum alternative fits remain poor (Jan2CO216.02,H2O249.06), but parameter rank may be degenerate and the output records rank rather than blindly counting columns. Formal chi-square values are diagnostics, not calibrated discovery significances. The assumed10% floor tests remain sensitivity calculations only.

## Independent observational evidence: bounded search, not a priority certificate
Search performed November16 2026 using exact designation, SWAN25B, CO2, infrared, SPHEREx, JWST, IRTF/iSHELL, CRIRES and spectroscopy/December/January combinations. No verified independent infrared spectrum or previous CO2 measurement of this target was located. Search-engine absence is not an archive completeness proof. MAST/ESO/IRTF moving-target holdings have **not** been comprehensively searched; no claim of zero holdings is warranted.

New useful primary-source lead:
- https://neoexchange.lco.global/target/156360/ identifies C/2025R2 as LOOK Project target and lists observation windows overlapping Dec21–22 and Dec24–25.
- https://neoexchange.lco.global/block/51401/ identifies LTP2025A-004,1m0 telescope,8x110s cadence, and last observation2025-12-30 19:37. This verifies actual observing activity in the relevant season, **not** that every scheduled window succeeded, nor that this is infrared/spectroscopy. No calibrated frames were retrieved here.
- Public archive query attempted: https://archive-api.lco.global/frames/?OBJECT=C%2F2025%20R2&start=2025-12-15&end=2026-01-05&limit=5 . Initial request timed out; a curl retry returned `count:16`, `count_estimated:true`, but an empty `results` list and a next-page link. This is inconclusive public metadata, not16 verified accessible exposures and not evidence of no data. The raw response is retained in `lco_archive_query.json`. No outreach or authentication attempted.

Previously verified optical observations remain context only:
- Jena Oct18 optical spectrum: https://www2.astro.uni-jena.de/Users/markus/C2025R2.html (C2,[OI],NH2; not aCO2 measurement).
- TRAPPIST-North observations on **2025-11-27 UT**, reported in [ATel 17538](https://orbi.uliege.be/bitstream/2268/338585/1/ATel%2317538.pdf) on December 10, provide independent **optical activity context**, not infrared CO2 confirmation. The preliminary OH, CN, C2 and C3 production estimates and phase-corrected dust proxies use profiles at 10,000 km and a Haser model with parent and daughter speeds of 1 km/s. These model assumptions and spatial sampling differ from our SPHEREx apertures, and the observations precede our December 22/25 spectra by 25/28 days. Do not combine them directly into a CO2/H2O abundance ratio or treat them as contemporaneous validation of the infrared feature.

Next productive independent check is to retrieve actual LOOK images near the SPHEREx dates for coma morphology/dust evolution, while separately auditing infrared archive holdings with target aliases and ephemeris-aware footprints. Optical data alone cannot confirm a4.26µm molecular identification.

## Reproduce
```sh
for m in co h2o; do
 curl -A 'Mozilla/5.0' --max-time 150 -sS -d type=rad \
 --data-urlencode file@psg/${m}_alternative_100K.cfg \
 https://psg.gsfc.nasa.gov/api.php -o psg/${m}_alternative_100K_response.txt
done
python fit_psg_alternatives.py --input psg_fit_input_snapshot.json \
 --output psg_alternatives_fit.json --calibration-root calibration --model-root psg
```
Requests should be sequential. The saved spectra and hashes permit reproduction without regenerating live PSG outputs.
