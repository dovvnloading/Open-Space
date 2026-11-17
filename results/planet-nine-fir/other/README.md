# Additional published candidates

Phan et al. report 13 pairs passing initial cuts and one candidate after image inspection. The published final pair is treated in the [main report](../README.md). [Phan et al.](https://doi.org/10.1017/pasa.2025.10024)

We transcribed the 20 peaks printed in Tables 3 and 4 of the ACT search. The authors found no significant excess over the sign-reversed search. Their ten Planet Nine-like peaks and ten general peaks are follow-up hypotheses. The reported sky coordinates are heliocentric on MJD 57688. The script extrapolates the authors' tangent-plane motions and distances and subtracts Earth's position at each date. Output precision is computational. The published table rounds the reference coordinates, distance and velocity, and supplies no per-peak covariance. [Naess et al.](https://arxiv.org/abs/2104.10264)

Chen et al. report separate July 2006 and July 2007 AKARI targets. They are different candidates, each observed on one day. We propagated example circular orbits with 72 headings and 11 distances from 400 to 900 au per target, including annual parallax. Rowan-Robinson gives an approximate parallactic center (ecliptic J2000 longitude 13.62°, latitude 71.64°), distance 225 ± 15 au and very uncertain position angle. We used the journal coordinates, 7 distances, all 72 headings and a July 1, 1983 reference date. The resulting annulus omits center uncertainty (quoted ±0.5°) and eccentric solutions. None of these grids has posterior weights. Rowan-Robinson found no convincing sequence in Pan-STARRS and documented cadence and cirrus limitations. [Chen et al.](https://arxiv.org/abs/2506.12854) · [Rowan-Robinson, journal version](https://academic.oup.com/mnras/article/510/3/3716/6428411)

Run candidates.py to regenerate the monthly bounding coordinates for each target in 2023–2026. Install NumPy and pyerfa and keep the adjacent orbits directory. The component-wise extrema for an annulus include a large interior that no sampled orbit reaches.

## Source receipts

| File read, November 17, 2026 | SHA-256 |
|---|---|
| https://arxiv.org/pdf/2104.10264 | ee8ec8d2ad587cf0b7f3b135732a6c88a47b3d1577c4eb7149af8c16e07212b0 |
| https://arxiv.org/pdf/2506.12854 | a613fde9ccc24904e98c15a49961fd30588df24b0c9875a94ed5858a18b7fb1b |
| https://arxiv.org/pdf/2111.03831 | 538e559c365592293db156eecd135541fe2cb1a65e52e7b3e33d91a7235e0a90 |

## Literature boundaries

Sedgwick & Serjeant found cirrus contamination in their final IRAS/AKARI sample. Liseau et al. initially reported a one-epoch submillimeter source near Alpha Centauri. Follow-up ALMA observations redetected it at essentially the same sky position. The change from December 2016 to March 2017 is only 0.088″. At any heliocentric distance in 500–1,000 au throughout that interval, the projected Earth displacement would require the source to move farther than a bound body can travel. Vlemmings et al. discussed a 12–25 au bound-object interpretation for Gna and an unbound interpretation at larger distances. The authors withdrew the claim pending further detections. [ALMA image fits and orbital-motion tests](../extra_candidates/README.md) · [ALMA follow-up](https://arxiv.org/abs/1904.03043) · [Gna version history](https://arxiv.org/abs/1512.02650) Search terms included IRAS/AKARI candidate, Planet Nine follow-up, ACT Planet 9 candidates and ALMA Gna / Alpha Centauri. This is a finite literature inventory. [Liseau et al.](https://arxiv.org/abs/1608.02384) · [Vlemmings et al.](https://arxiv.org/abs/1512.02650)

Rice & Laughlin's 17 TESS detections were reported as candidate distant TNOs. Their authors state the realized search did not reach Planet Nine. Socas-Navarro & Trujillo's January 2026 optical parallax search reports no credible candidate for targeted follow-up. Sedgwick & Serjeant attributed their final sample to possible cirrus. We have not assigned new object-level exclusions to those records. [TESS paper](https://doi.org/10.3847/PSJ/abc42c) · [Authors’ TESS summary](https://www.astro.yale.edu/malenarice/researchhighlight4.html) · [2026 optical search](https://doi.org/10.33232/001c.155004) · [Sedgwick & Serjeant](https://doi.org/10.1093/mnras/stac2044)

## Euclid long-VIS footprint inventory

We downloaded the 26,724 DR1 long-VIS footprint records from ESA. We propagated each secondary candidate to the exposure midpoint after a 1.6° daily shortlist. Projected footprint polygons include four exposures on the nominal ACT Table 4 peak 6 trajectory. The other 19 ACT tracks and the screened Chen and Rowan-Robinson circular models have no polygon hits in this product selection. This geocentric inventory omits a spacecraft parallax correction near boundaries, individual CCD masks, shorter images and NISP. A separate [native-pixel search of the four ACT4_06 images](../act_vis/README.md) measures the sensitivity for that track and its offset patch. Run `python euclid_coverage.py /path/to/cache.csv` in this directory. The query receipt gives the original response hash; the intersection record retains matched polygons, epochs and trajectory IDs. [Footprints](euclid_footprints.csv) · [Query receipt](euclid_footprint_query.json) · [Intersections](euclid_intersections.jsonl.gz)
