# IceCube reported this GRB coincidence on 2024-04-24

IceCube's discovery circular explicitly names GRB 240419B inside the 90% region of IC240424A and reports an offset of about 1.7°. The BAT burst preceded the Gold track by 4.5601119 days. The published containment statement establishes a known positional coincidence. The report supplies no evidence of a causal association or post-trials significance. [IceCube GCN 36283](https://gcn.nasa.gov/circulars/36283)

| Measurement | ICRS RA, Dec (deg) | 90% radius | Offset to track best fit | Offset to final XRT |
|---|---|---|---:|---:|
| Frozen Swift GRB-table BAT position | 325.853, 4.299 | 3.3 arcmin | 1.742° | 4.662 arcmin |
| Mission BAT catalogue reanalysis | 325.8913, 4.232397 | 1.456 arcmin | 1.668° | 0.351 arcmin |
| Final astrometrically enhanced XRT afterglow | 325.89670, 4.23468 | 2.2 arcsec | 1.666° | 0 |
| Nordic Optical Telescope afterglow | 325.89704, 4.23481 | Unspecified | 1.665° | 1.30 arcsec |

The XRT centroid falls outside the early BAT radius and inside the later mission BAT radius. The optical identification agrees with the final XRT position. The paper search uses the catalogue position chosen before matching. These follow-up positions serve only for source verification. [Early BAT analysis](https://gcn.nasa.gov/circulars/36292) · [BAT catalogue](https://swift.gsfc.nasa.gov/results/batgrbcat/GRB240419B/web/GRB240419B.html) · [Final XRT localization](https://www.swift.ac.uk/xrt_positions/01223072/) · [NOT photometry](https://gcn.nasa.gov/circulars/36191)

| Epoch relative to BAT trigger | Band | Measurement |
|---|---|---|
| 2024-04-19 12:22:52.33 UTC | BAT | Multi-peaked burst lasting about 150 s, circular T90 132.16 ± 26.29 s in 15–350 keV. |
| +141 s onward | XRT | Early bright emission, then a fading X-ray afterglow. The initial refined analysis found a late decay index of 1.02 with uncertainties +0.12, −0.11. |
| +159 s onward | UVOT | Optical/UV upper limits, including white >20.4 mag reported for the initial series. |
| +0.698 d, +0.701 d | NOT optical | r = 22.01 ± 0.10 and z = 20.81 ± 0.19 AB mag. The counterpart is absent in deeper pre-burst images. These are different filters, so the two points alone do not measure a decay slope. |
| Last XRT bin ends +1.39384 d | XRT | The deposited XRT curve contains 165 bins. Its final bin is centered at +1.38431 d. |
| +4.56011 d | IceCube | No simultaneous XRT bin. The neutrino fell well after the prompt gamma-ray emission. |

A prompt gamma-ray/neutrino coincidence is excluded by the measured delay. The available points cannot test activity at the neutrino epoch or establish a delayed physical origin. The independent IceCube follow-up reported no significant multiplet excess in its tested windows. Fermi-LAT examined other possible counterparts in the neutrino region, including a gamma-ray transient, so this field admits alternatives. [Swift discovery](https://gcn.nasa.gov/circulars/36180) · [XRT follow-up](https://gcn.nasa.gov/circulars/36192) · [UVOT](https://gcn.nasa.gov/circulars/36197) · [MITSuME](https://gcn.nasa.gov/circulars/36227) · [IceCube multiplet search](https://gcn.nasa.gov/circulars/36309) · [Fermi-LAT](https://gcn.nasa.gov/circulars/36349)

The joint-designation search on 2026-11-17 included GCN, indexed ATels, and arXiv. It recovered IceCube's own prior report. No positive claim depends on an absence of search hits. The CSV lists related circulars for the two event names, including reports about alternative counterparts. The parent receipts hash the complete GCN archive. Local receipts identify XRT and BAT products. Run `python verify.py` for separations, the delay and the light-curve time coverage. [GCN search](https://gcn.nasa.gov/circulars?query=240419B)
