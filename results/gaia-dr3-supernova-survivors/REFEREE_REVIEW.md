# Independent referee record — 17 November 2026

**Disposition: accept as a bounded catalogue null with explicit low completeness, not as a Galactic survivor constraint. No validated new source, new chemical abundance pattern or new supernova-remnant association.** The lead Gaia DR3 6094201393481158656 remains unresolved. No [REMARKABLE] label or addition to the remarkable findings log is justified.

## Executed replay

An independent GitHub Actions runner checked out the study on `research/gaia-dr3-supernova-survivors` at commit `f7b35c32c38f0265a5df8f4924748b7cf94389bb`. [The two-job referee run](https://github.com/dovvnloading/Open-Space/actions/runs/44050032886) passed; the subsequently amended prose and root catalogue do not change the numerical data. Python 3.12 used NumPy 2.3.5, SciPy 1.17.0 and gaiadr3-zeropoint 0.1.0, plus Astropy 8.0.2, Gala 1.12.0 and galpy 1.12.1 for orbits.

The `frozen-inputs` job ran `sha256sum -c SHA256SUMS`, the newly independent `referee_check.py`, `control/recovery.py`, `search/analysis.py` and `search/inject.py`, and finally both checksum and independent assertions again. The `orbit-and-remnant` job reran `origin/run.py` on the included Gaia tables, then fetched the SHA-256-pinned October 2024 Green supernova-remnant catalogue and reran `origin/remnants.py`. **Both jobs regenerated their deposited outputs byte for byte, with no hash differences.** No simulation was interpreted as a real detection.

The independent check verifies every deposited entry in the SHA-256 file manifest, hashes of the locally deposited Gaia TAP responses and query receipts, the valid VOTable ID enumeration, the original two orbit batches, and each source's astrometric threshold using a separate covariance/velocity-root calculation. The Green archive SHA-256 was checked anew by the remnant script (source hash `0dcffa17be68c3f2467d5e80eae40c6dddba68195b5a4765556c6b4e9f966396`). Query receipts keep exact ADQL, endpoints, times, counts and checksums. Release catalogues and Gaia bulk batches that were **not checked into Git** are identified by input hashes and a rerun recipe, but this review did **not** independently re-download and re-crossmatch all >113,000 external source rows: its union counts are recomputed from the fifteen deposited, mutually exclusive survey-overlap records. This is a **receipt-backed internal replication**, not an independent second survey reduction.

## Claim checks

| Quantity checked | Independently computed result | Status |
| --- | ---: | :--- |
| Spectroscopic union, after deduplicating four release membership sets | 113,166 distinct Gaia IDs | Pass |
| Quality and corrected parallax S/N ≥5 | 67,433; 0 above 1,000 km/s | Pass |
| Quality and corrected parallax S/N ≥3 | 75,684; 0 above 1,000 km/s | Pass |
| Largest measured transverse speed at either cutoff (Gaia DR3 3688712561723372672) | 704.555 km/s | Pass |
| All-sky Gaia parent before speed cut; independent TAP enumeration after it | 1,023,709; 12 sources | Pass |
| Above 1,000 km/s at a one-sided 95% lower bound | D6-1 and D6-3 only | Pass |
| Fastest lower bound among the other nine, Gaia DR3 6094201393481158656 | 810.750 km/s (nominal corrected point 1,117.448 km/s) | Pass |
| Same target's closest 3-arcsec Gaia neighbour | 1.3492 arcsec | Pass |
| Initial published-control recovery | D6 3/11; LP40 0/8; He 0/1 | Pass |
| D6 objects missed at raw Gaia parallax S/N >3 | 8/11 (all G ≈19–20.5) | Pass |
| 1,500-km/s injected survivors at 17 ≤ G <18 | 56.2302% at 1 kpc, 36.2970% at 2 kpc, 1.8880% at 4 kpc | Pass |
| D6-2 to SNR G70.0−21.5 minimum projected separation | 0.90469°, 93,750 years ago; shell radius 2.75° | Pass |
| Galactocentric median speed of D6-1 / D6-2 / D6-3, geometric distances | 1,991.22 / 1,125.55 / 2,303.29 km/s | Pass |

The 67,433-source zero gives an exact-binomial one-sided 95% upper fraction of **44.4 per million measured exceedances**, *only under an additional i.i.d./exchangeable-observation model*. The raw sample and the spectral catalogues are nonrepresentative: many D6 runaways have weak or negative parallaxes and almost all of the confirmed D6 population lies outside the four specialized white-dwarf spectral tables. The all-sky selection also rejects radial-dominated motions and sources with proper motion ≤35 mas/yr, and its absolute-magnitude proxy does not extinction-correct. This is therefore **a weak all-sky null** and is not a bound on thermonuclear-supernova rates.

The high-speed unresolved object is not an exception: its 95% lower limit of 811 km/s does not pass the proposed threshold. It has a fainter neighbour 1.35 arcsec away and Gaia's ordinary-cool-star spectral interpretation. The decisive follow-up is high-resolution imaging to quantify the blend, spatially resolved two-epoch optical spectroscopy to measure systemic RV and chemistry, and better five-/six-parameter astrometry (Gaia DR4 or independent). A spectral resemblance or apparent high raw velocity by itself does not prove a compact runaway.

## Additional cross-checks and limits

The checked-in chemistry and hot-subdwarf summaries agree with the top-level README's 1,452 metal/carbon-selected objects, 2,791 hot-subdwarf entries with catalogue RVs, and 247 single-lined systems with complete local comparison velocities. Four DESI Mg lines for J1109+0001 have a published fit of +96.5 ± 7.5 km/s in the retained spectral tables. These additional fits and the full upstream archive census were **read and compared to the deposited summaries, not independently rerun** during this referee pass. As already stated in the research README, the shell intersection is a **projected** association reported previously by Shen et al., not a new causal birthplace. Spectroscopic RVs, Gaia distances and subsequent supernova-remnant orbits are not measured for the nine new astrometric leads.

This review requests **no change** to `docs/REMARKABLE_FINDINGS.md`, and the PR title remains an ordinary bounded search rather than a discovery claim.
