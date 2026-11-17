# Radial-velocity search for US 708 analogues

**US 708 is the only ≥700 km/s absolute heliocentric RV in the 2,791 measured entries of the Culpan et al. hot-subdwarf catalogue. Its revised Keck/ESI velocity is +835 ± 10 km/s (El-Badry et al. 2026).** A DESI DR1 MWS crossmatch adds two nominal high-RV targets. Both have flagged, conflicting fits and neither qualifies as a confirmed survivor. In the 2026 compilation within 500 pc, the largest measured Galactic-rest-frame speed is 302 km/s among 247 objects.

| Independently evaluated product | Coverage | Observed result |
| --- | --- | --- |
| Culpan et al. known hot subdwarfs | 2,791 of 6,616 rows have at least one SDSS or LAMOST RV | US 708 at +793 ± 16 km/s. Next largest absolute values are 510 and 506 km/s. |
| DESI DR1 MWS | 781 spectra of 668 exact-ID matches, 415 objects have at least one RVS_WARN=0 fit | Two nominal hits with warning flags. Maximum clean absolute RV is 397.2 km/s. |
| SDSS-V DR19 and DR20 Astra star products | 21 DR19 matches, zero with usable RVs. DR20 has 345 matched stars and 304 with finite RV and positive uncertainty and number of good RVs. | No ≥700 km/s hit. The maximum DR20 absolute RV is 396.9 km/s. The releases overlap. |
| Dawson et al. 500 pc sample, CDS tables 2 and 4 | 253 tabulated RVs, 247 complete Galactic-rest-frame velocity triples | Maximum absolute RV is 281 km/s (BD+39 3226). Maximum 3D speed is 301.86 km/s (PG 1519+640). |

Culpan contains 1,958 SDSS RVs and 1,316 LAMOST RVs with 483 overlaps. Every row with an RV is distinct by its catalogue record number. Among the 2,791 rows with RVs, 2,790 have distinct Gaia identifiers and one has no identifier. Across the full catalogue, 6,593 rows carry a Gaia source ID and 6,552 of those IDs are distinct. The newer survey joins used the distinct IDs. All three nominal high-RV objects each have a unique Culpan row. The original class assignments can include contaminants. The SDSS comparison uses the source-level RV fields in the Astra SnowWhite product irrespective of whether the white-dwarf classifier produced a label.

## Flagged DESI sources and crosschecks

- Gaia 3670512036510958208 has adopted velocities +1500 ± 64 and −792 ± 260 km/s from two DESI targets with warning flags 14 and 5. The pure RVSpecFit column of the +1500 fit itself allows ±430 km/s. An independent, position-matched legacy SDSS spectrum (plate 531, MJD 52028, fibre 564) gives cz ≈ −71 ± 57 km/s. SIMBAD lists this as a DA: white dwarf. Its measured Gaia proper motion is 2.91 mas/yr.
- Gaia 1391315863395191168 has adopted values −833 ± 797 and −100 ± 58 km/s in two coadds, flags 5 and 8. Legacy SDSS plate 1052, MJD 52466, fibre 619 gives cz ≈ +128 ± 113 km/s. SIMBAD lists it as a DA white dwarf. Its measured proper motion is 2.69 mas/yr.
- All these pipeline RVs can be affected by template mismatch. The discordant estimates provide insufficient evidence for a persistent high centre-of-mass velocity, peculiar composition or a new remnant. No new object advances to a Galactic orbit.

## US 708 spectrum and control limitations

We retrieved US 708's SDSS legacy spectrum, plate 870, MJD 52325, fibre 575. We fitted a local linear continuum and independent nonnegative Gaussian absorption depth and width for each line, sharing a relativistic Doppler shift. Eleven starting velocities sample the interval −1,500 to +1,500 km/s. The SDSS coadd supplies inverse variance and its log vacuum wavelengths in the heliocentric frame. The air-to-vacuum conversion and relativistic factor come from the LP 40-365 control fitter. The eight selected He I/He II features yield +772 ± 18 km/s formal. The four He II Pickering lines near Balmer wavelengths give +828 ± 23 km/s. Individual line solutions and the SHA-256 of the spectrum are in us708_fit.json.

**The current Keck/ESI benchmark is +835 ± 10 km/s.** El-Badry et al. (2026, Appendix “An Updated Radial Velocity for US 708”) re-reduced ten archival ESI exposures and measured +818 ± 15 km/s from a newer ESI observation. They suspected an air-to-vacuum conversion error behind the roughly 80 km/s historical discrepancy and revised the inferred ejection speed to 802 km/s. Geier et al. (2015) had reported +898 ± 30 km/s on SDSS and +917 ± 7 km/s on ESI, and showed that hydrogen-rich templates can confuse the Pickering and Balmer series. Our eight-line Gaussian result falls 63 km/s below the updated ESI value. The Pickering-only fit agrees with it to within formal errors. These fits recover the threshold control and demonstrate line-model systematics that the formal uncertainties omit. Gaia's parallax S/N is 0.33, and its measured proper motion is only 5.73 mas/yr, explaining the miss in the transverse-speed search. [El-Badry et al. (2026)](https://arxiv.org/abs/2606.11293) · [Geier et al. (2015)](https://arxiv.org/abs/1503.01650)

## Local conditional limit

For the 247 measured stars in Dawson et al. table 4, zero have Galactic-rest speed at or above 1,000 km/s. Inverting a binomial model gives a one-sided 95% upper fraction of 1 − 0.05^(1/247) = **1.21%** for an exchangeable population with the same object definition and complete recognition above that threshold. The catalogue gives a finite census of the reported velocities. The statistical interpretation additionally assumes representative sampling. The analysis selected apparently single-lined systems and excluded composite spectra, and six of the 253 table-2 objects lack a published table-4 vector. Extreme RVs may affect the original spectroscopic classification. Individual binary motions complicate estimates of some systemic RVs. This bound applies to the selected local hot-subdwarf sample and cannot be converted here to a Galactic ejection rate or a constraint on the D6 WD population. The authors used Gaia and spectroscopic measurements to derive U, V and W under their stated solar motion. We evaluated their tabulated vectors. [Dawson et al. (2026), section 5.2](https://arxiv.org/html/2601.02250v1)

## Reproduction

Run the following with Python, NumPy, SciPy and Astropy. The working-directory cache stores the full Culpan table and the legacy spectrum. SURVIVOR_CACHE can point to another cache. VizieR stamps its output header with the access date, so an otherwise identical future response can have a different byte hash. The first command checks all published record numbers. The crossmatch queries, input hashes and threshold outputs are included here.

```sh
python screen.py --download
python desi.py --download
python nearby.py
python fit_us708.py
python legacy.py
python sdss.py /path/to/verified/astra/cache
```

The last command uses the verified files supplied by the adjacent SDSS study's download step. It verifies the SHA-256 of the compressed products before reading their extracted counterparts. The archived Gaia and SIMBAD input queries and receipt files give their exact-ID checks. Survey inputs come from [Culpan et al. (2022)](https://doi.org/10.1051/0004-6361/202243337), [DESI DR1 at Data Lab](https://datalab.noirlab.edu/tap), [SDSS DR20](https://dr20.sdss.org/) and [Dawson et al. (2026)](https://doi.org/10.1051/0004-6361/202558123).
