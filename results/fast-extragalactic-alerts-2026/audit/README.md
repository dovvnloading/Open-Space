# Thirty-two alert dispositions and referee audit

The primary cohort is strictly first-broker-detection objects in the **2026-11-01 through 2026-11-16 UTC** Rubin and ZTF public alert streams, with the two-positive-same-band rapid-pair cuts frozen before the search. This is an audit of **alert leads**, not a complete census of astrophysical objects. No discovery-level fast extragalactic event is validated.

| Evidence disposition | Rubin | ZTF | Total |
|---|---:|---:|---:|
| Pre-existing point source | 15 | 0 | 15 |
| Bogus / flagged subtraction pair | 6 | 2 | 8 |
| Likely stellar flare / variable | 0 | 5 | 5 |
| Asteroid (two IDs for one body) | 0 | 2 | 2 |
| Unresolved | 2 | 0 | 2 |
| **Total** | **23** | **9** | **32** |

[Every object, coordinates, specific evidence, Gaia/PS1 separations, TNS and GCN comparisons](dispositions.csv) · [Query receipts and raw-download SHA-256 values](query_receipts.json) · [Bounded replay script](validate.py).

**Interpretation limits.** A Rubin label of “pre-existing point source” means an independently catalogued positional counterpart within 1 arcsec (or, for source 171198361324486700, the recorded Legacy Surveys counterpart 0.11 arcsec away). A low-significance Gaia parallax does **not** establish stellar or extragalactic origin. “Bogus” rejects a rapid *pair* because at least one selected Rubin difference detection has an isdipole flag, or (ZTF) a documented repeated subtraction/dipole residual; it does not claim the underlying sky source is fictitious. The five ZTF stellar categories are probable, not spectroscopically measured for every object. Full ZTF circumstances appear in [ZTF.md](../ztf/ZTF.md).

**Two unresolved Rubin objects.** 171167577750896661 (RA 64.8702045°, Dec −46.9643867°) has no Gaia match in a 10-arcsec cone and lies outside PS1 coverage. The official 2026-11-17 00:00 UTC TNS public snapshot identifies **AT 2026aidn**, 0.212 arcsec away, but neither a spectral class nor an extragalactic distance is established. The independent broker TNS alias agrees. 171198363518632040 (RA 348.5270562°, Dec −11.8414068°) has no Gaia match; its nearest PS1 detection is 1.017 arcsec away, slightly beyond the conservative one-arcsec counterpart threshold. The same official TNS snapshot has **no object within 10 arcsec** of this alert. Neither unresolved designation is a claimed orphan afterglow or Cow-like event.

**Trigger coincidence.** The archived November 1–16 GCN inventory and chosen localization/time windows produce 234 sky-time comparisons for the 32 candidates. No supported GCN association emerged. Two ZTF positions on asteroid (18150) enter a broad GBM probability region but lie 13.46 degrees from the narrower SVOM 90% localization. This is not an exhaustive search of every gamma-ray, X-ray or neutrino notice or sub-threshold event. [Trigger limits and receipts](../coincidences/README.md).

**Independent referee arithmetic and receipts.** Recomputing the selection cuts from the stored pairs confirms **47 Rubin pairs / 23 objects** and **nine ZTF pairs / nine objects**, with no passing-pair threshold failures. The stored funnel values and ZTF program-1 checks agree. The replay verifies all **34 SHA-256 input hashes** in the turnover audit, and the Rubin raw HTTP-reply digests matched their gzipped CSVs. The 555-member ZTF evidence supplement contains **319 receipt-labelled members** and a receipts/ztf/ directory with the original TAP query receipts and reply digests. The separate ZTF program-verification ZIP is also hashed in [referee_checks.json](referee_checks.json). The full official TNS ZIP (213,398 data rows) and every Gaia, PS1 and secondary broker request have status, retrieval time, byte count and SHA-256 recorded in [query_receipts.json](query_receipts.json). The official TNS archive is mutable and is not copied into this repository.

**Recovery and pending rate denominator.** At AB 18 and a 2-mag/day rise/fade, the two-point injection recovery was 4.27% (Rubin) and 0.436% (ZTF) *conditional on the recorded alert epochs*. For a measured rise **and** fade within a day, only 244 Rubin and nine ZTF cohort objects have the needed temporal visits; no observed triple passed. These fractions are **not survey exposure efficiencies**. No independent footprint, area-time/volume-time, trigger-generation or host-vetting completeness is measured, so **no numerical effective exposure or transient rate limit can responsibly be computed**. Those calibrations, longer-baseline photometry and classification of the two unresolved Rubin objects remain pending. [Recovery design](../recovery/design.md).

**Referee verdict:** reproducible, scoped null finding within the frozen broker-new window; no discovery-level designation or claim.
