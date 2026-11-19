# Auger open-data UHECR sources: reproducibility-first search

**Claim (initial, unvalidated): no new cosmic-ray source or energy-ordered multiplet has been established.** This directory tracks an attempted independent search; no unrun comparison may be described as a null result.

## Dataset and scope

- Auger open data: [current-version DOI 10.5281/zenodo.4487612](https://doi.org/10.5281/zenodo.4487612); latest release listed by the Auger portal is **release 3, 2024-03-20, DOI 10.5281/zenodo.10488964**.
- The release is approximately **10% of quality-selected Auger cosmic-ray events through August 2018** (not 10% of all subsequent Phase-I or Phase-II measurements); additionally includes **the separately published UHECR event catalog**. The event subsets must not be mixed or treated as independently sampled without checks.
- The 2012 Auger multiplet reference used >20 EeV, a >45 EeV seed, a 20-degree extent, transverse half-width <1.5 degrees, and inverse-energy directional Pearson r>0.9, and found no significant multiplet (12-plet, 6% isotropic p).

## Status (update on every step)

- **Completed:** DNS reachability survey in [host-connectivity receipt](receipts/host-connectivity.md) — all seven relevant hostnames failed DNS in the isolated computation environment; read the connectivity limitation before any attempted inference.
- **Pending:** independently acquire and SHA-256 hash the event files into computation environment; validate columns, cuts and declination exposure; reproduce dipole and Centaurus A / starburst excess controls; run energy-ordered multiplets; compare LHAASO, magnetars, young pulsars, local-volume galaxies and 4LAC-DR3 against exposure-preserving sky scrambles and energy permutations, with family-wise trial calibration; independent referee check.
- **Not a calibrated null:** as of this initial commit, no real-event hypothesis tests have run.

## Reproducibility receipts

The host test was executed locally, **not in GitHub Actions**. All claims depend on future public-data hashes, independent selection masks, exposure accounting, published reference comparisons, Monte Carlo settings and output receipts. Do not apply `[REMARKABLE]` unless a discovery survives validation under [repository policy](../../docs/REMARKABLE_FINDINGS.md).
