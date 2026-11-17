# Independent referee review of PR #142

**Verdict: Accept as a bounded null, not a discovery.** As of 2026-11-17 the ten-field search identified no validated uncatalogued planet or brown-dwarf companion at the reported tested separations and flux ratios. This does not rule out closer/fainter companions, other bands, unsearched epochs or every target in the wider archive. No remarkable-finding label or log entry is warranted.

## Independent numerical audit

The [pure-standard-library referee tests](../../tests/test_jwst_nearby_referee.py) reconstruct key claims from the deposited records without taking the narrative summary as input.

| Audit | Regenerated finding |
|---|---|
| Frozen MAST inventory | 9,056 primary exposure filesets in 2,653 visit/optical-element groups from 240 programs. The nearby-host table has 792 matching groups, 331 distinct Gaia hosts in 70 programs. The object/title literature screen has 214 nonempty matches; absence of a title match is not absence of a paper. |
| L 788−34 | Fourteen Pan-STARRS detections on nine nights. Median stationary residual 0.252157 arcsec; host-following residual 11.133507 arcsec. The remaining source is extended in JWST and VISTA; its morphology and position support a background galaxy interpretation rather than a new planet. |
| G29−38 | Pinned Spitzer IRAC2 2004 and 2009 centroids offset by 0.069706 and 0.055268 arcsec from the JWST fixed-sky locus; host-following predictions differ by 10.003 and 7.834 arcsec. There are 24 Gaia registration ties per epoch. |
| GJ 205 | Two masks were checked. Gaia field-star relative prediction (4.043889, 3.705548) arcsec vs measured (4.019, 3.661) arcsec; 0.051029 arcsec residual. Host and field-star parallaxes 175.31 and 1.14 mas, respectively. |
| Ross 619 | Two F444W peaks score 5.766 and 6.269 with primary references but only 1.307 and 3.073 with the separate reference set. Neither meets the alternate-library threshold. |
| Two white-dwarf controls | WD 1202−232: four stationary and seven unmatched peaks; WD 2105−82: five stationary. No peak tracked its host in the qualifying repeat test; unmatched is not a background or planet classification. |
| Contrast and recovery | Five NIRCam F444W fields: six radii each and 2,376 trials each (11,880 total). Both GJ 205 MIRI masks: six radii and 3,168 trials each (6,336 total). Also checked 1,602 G29−38 trials, 1,220 Ross 154 trials, 1,920 white-dwarf control trials and their tabulated 90% thresholds. All recomputed 95%-recovery levels matched the underlying NIRCam and GJ 205 trial ledgers. |

## Independent clean-machine pixel replay

The [GitHub Actions JWST FITS referee replay](https://github.com/dovvnloading/Open-Space/actions/runs/44043395741) completed **four passing jobs** after two checker-only fixes (numeric-grid rank comparison instead of exact float-string keys, and directing one alternative-reference script to its scratch output folder). The runner re-fetched FITS via pinned source URLs, compared SHA-256 with the committed manifests, and re-executed:

- Spitzer image registration and forced-centroid motion measurements for G29−38 (2004 and 2009 IRAC2).
- The MIRI F2100W/F2550W G29−38 reduction, source count and injection curves.
- The NIRCam F444W L 788−34 and Ross 619 RDI searches, injection curves and alternate-reference source checks.
- The GJ 205 MIRI F1065C/F1550C PSF-subtracted searches and injection curves.

The first attempt also demonstrated **no image-reduction mismatch**; its two failing checks were defects in the new referee wrapper, corrected and rerun to all-pass. The [replay code](referee_pixel_replay.py) and the [CI workflow](../../.github/workflows/jwst-nearby-fits-referee.yml) are committed for repeatability.

## Receipt and scope qualifications

MAST search receipts preserve per-page offsets, row counts, release-date predicates and response digests. Input manifests preserve source URLs, sizes and SHA-256 digests. Where downloaded in the clean-runner replay, FITS SHA-256s were rechecked against the frozen manifests. **Not all inputs for every field were independently downloaded again:** V374 Peg, G 99−49, the MCC124 control, Ross 154 and the two white-dwarf mosaics were independently checked at the *intermediate numerical table/trial* level only. Their original FITS hashes remain claims in the deposited manifests, not freshly authenticated by this review. Original HTTP response SHA-256 receipts were reviewed structurally, not live re-fetched byte for byte; a later MAST query may change as the archive updates. Off-axis calibrations with blank rights flags are *unknown metadata*, not public-rights claims.

No bootstrap confidence intervals, survey-wide companion occurrence-rate limit or mass limit follows from these contrast curves alone. The methods and acceptance regions vary by band, field geometry, mask, usable azimuths and recovery criterion. Multi-epoch claims require same-source identification, not merely detecting a faint peak. The source ledger retains unresolved/unmatched controls as pending checks.

## Prior-paper exception and next work

The latest target-specific literature is not exhausted by MAST Papertrack or a frozen SIMBAD title query: in particular published searches already cover HR 2562, HD 92945 and HD 107146. The original snapshot returned papers only through June 2026. These fields are **not** first-priority unpublished targets. The five highest-value *pending* unreduced targets are prioritized in the [claim-first README](README.md#highest-value-unreduced-targets--pending-not-null-results): LAWD 66, HD 53143, HD 38206, HD 21997 and HD 202628 (the last requires an explicit JWST-paper applicability check).

Root catalog contains one plain-language *Null result* entry for this work, with the count updated. No edit was made to docs/REMARKABLE_FINDINGS.md; no [REMARKABLE] title.
