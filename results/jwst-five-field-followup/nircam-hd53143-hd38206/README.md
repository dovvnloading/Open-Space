# F444W imaging of HD 53143 and HD 38206

**HD 53143 recovers a previously rejected background object. HD 38206 has no F444W point-source threshold crossing that repeats at the expected sky position in both 2025 epochs.** The blind search covers 2.5–5.5 arcsec (46–101 au and 177–389 au in projected separation, respectively). It has no sensitivity statement for the inner working angles of the instruments, obscured regions, or the whole surrounding sky.

## Measurements

| Star and observation | Blind F444W crossings | 95% empirical recovery, µJy, at 2.65, 2.95, 3.55, 4.15, 4.75, 5.35 arcsec |
|---|---:|---|
| HD 53143, 2024 Apr 4, o001 | 8 | 11.94, 9.02, 19.23, 8.99, 1.15, 0.92 |
| HD 53143, 2024 Apr 4, o008 | 10 | 16.82, 5.06, 21.11, 9.97, 1.75, 1.11 |
| HD 38206, 2025 Jan 23, o141 | 5 | 16.29, 5.81, 5.30, 4.12, 5.68, 4.81 |
| HD 38206, 2025 Oct 24, o150 | 3 | 6.84, 5.16, 9.04, 8.12, 5.28, 4.57 |

Each row is an independent invocation of the merged F444W pipeline. Five contemporaneous same-program PSF dithers constitute each visit's reference library. The archived off-axis G191-B2B images, CALSPEC spectrum, F444W transmission curve and coronagraph transmission map have the original hashes. Each curve uses 36 angles and 11 fluxes at six radii (2,376 injections). A full raw-cube injection through refitting agrees with the accelerated calculation to <4.0e-13 in the reducer's amplitude units for each exposure. Annular finite-score coverage is 99.55–99.64%. The search and injection definition, including the treatment of pre-existing peaks, are inherited from PR #142. The empirical 95% value is the lowest tested level with recovery at or above 95% also at all brighter tested levels. The tables retain angle-level results and flux contrasts. Stellar F444W normalizations are approximate Ks-anchored Planck values, 2.323 and 0.937 Jy. The absolute µJy limits do not depend on these assumed host fluxes.

**HD 53143.** Independent empirical-PSF fits give 91.74 and 91.47 µJy with core correlations 0.998 and 0.997. The fitted star-relative offsets are (−1.922, −2.772) and (−1.929, −2.789) arcsec east/north. At the 2009.03 Gemini/NICI position, Gaia DR3 stellar proper motion and observer parallax predict (−1.788, −2.838) arcsec for an inertial background source in April 2024. The residual is 0.149 arcsec in each roll, versus 4.67–4.68 arcsec for a constant star-relative position. Differential astrometric calibration, background-object motion and rounding of the NICI epoch limit the first comparison. The 15-year displacement and the published independent NICI/SPHERE background classification exclude a co-moving interpretation of this re-detection. Separately centering the short-wave PSF gives F200W fluxes of 268.5 and 268.3 µJy, SW fit correlations of 0.994, and F200W–F444W AB colors of −1.166 and −1.168. These are much bluer than the +3.77 to +8.55 range of the 67 flux-matched Bobcat nodes for o001. Independent centering was necessary because the SW core lies 0.18 arcsec from the nominal LW-to-SW projection in each roll. Both April rolls are on the same date. After translating the detector pattern with the fitted bright source, each of the eight o001 crossings matches an o008 crossing within 1.01 pixels; eight of ten o008 crossings match back. A single independently calibrated PSF at the bright source predicts 81% to 132% of the matched-filter amplitude at all 18 threshold crossings, including the two weaker unmatched lobes. Some residuals remain above threshold after that subtraction. No additional source is securely validated in these crowded residuals.

**HD 38206.** All eight crossings fail the other epoch under both predictions. At the co-moving predictions the highest score anywhere in a 3-by-3-pixel patch is 2.35. Two maxima repeat at identical detector pixels (98,215 and 166,235) while the sky rotates by about 119 degrees. Both central pixels carry DO_NOT_USE throughout each observation. The strongest October trigger (21.9 annular S/N) likewise sits on a DO_NOT_USE plus JUMP_DET pixel, and does not repeat on the sky. At the January location corresponding to the October trigger the forced flux is 0.63 µJy, compared with 28.35 µJy in October. Background versus co-moving predictions differ by 35.4 mas across these epochs, including observer parallax. The two-epoch absence of a repeat can reject a persistent source at the measured brightness where the other epoch is sensitive, but cannot classify a faint source below the sensitivity of the other epoch.

## Checks and limits

- The paired-roll table records each crossing under co-moving, stationary-background and detector-fixed predictions and supplies distances to the nearest other-epoch threshold crossing and non-host Gaia neighbor. Absolute astrometry uses the nominal coronagraph center. Gaia's source list is incomplete in these bright halos.
- The paired-color table is forced small-aperture F200W or F210M photometry at the nominal F444W location. Treat aperture sums as diagnostics. Saturated or misregistered positions, sky residuals and the short camera's smaller footprint can invalidate colors. The known HD 53143 source requires separate centering on the SW camera.
- The conditional control table reports JADES GOODS-S F444W source counts at or above each trigger's brightness over the 61.5-square-arcminute quality-selected control area. It scales counts to the geometrical search annulus (75.4 square arcsec) and gives Poisson intervals. The field provides extragalactic counts only. As a separate Galactic control, Gaia DR3 resolves 3 and 9 non-host objects in the 10–40 arcsec local annulus in HD 53143 and HD 38206, respectively. Scaling those optical counts gives 0.048 and 0.144 objects per search annulus, subject to glare, incompleteness and dissimilar bandpasses. The same table brackets SW–F444W AB colors over cloud-free Sonora Bobcat models at three metallicities with F444W within 25% of the measured flux at the stellar distance. Colors alone cannot distinguish a chance projection. A grid comparison at an instrumental peak does not identify an astrophysical source.
- Historical tests include both November 2005 HST/ACS F606W rolls for HD 38206, and both June 2024 unocculted F444W mosaics for HD 53143, with forced apertures at both hypotheses. The host core saturates in June and absolute WCS and structured halos limit these tests. Values are local halo scores, not Gaussian significances. The HST roll depths and previously published VLT/SPHERE and Gemini data supplement the two-epoch test. The ground-based images of HD 38206 have not been re-reduced here. HD 53143 MIRI and additional HST checks are documented separately.
- An independent fully distortion-calibrated analysis of the faint PSF residuals, a full-frame blind search in the June exposure, and an inner-disk search fall outside these measurements. No occurrence constraint or mass-completeness limit is inferred.

## Replay

Use Python with NumPy, SciPy and Astropy. Set `P` to the merged `results/jwst-nearby-companions` directory, `C` to a scratch cache and `O` to the directory containing the archived results. Run these commands in the directory containing the adapters.

```sh
python literature_query.py "$O"
python acquire.py --parent "$P" --cache "$C" --out "$O"
for v in hd53143_o001 hd53143_o008 hd38206_o141 hd38206_o150; do
  python reduce_fields.py "$v" --parent "$P" --cache "$C" --out "$O"
done
python check_peaks.py "$P" "$C" "$O"
python controls.py "$P" "$C" "$O"
python gaia_counts.py "$O"
python known_source.py "$P" "$C" "$O"
python get_wide.py "$C" "$O"
python arch382.py "$O"
python fetch_hst.py "$O" "$C"
python archive_tracks.py "$P" "$C" "$O"
```

The `inputs.csv` manifest records the exact bytes used. To reproduce a frozen revision, download through the parent `nircam/nircam_fetch.py --manifest "$O/inputs.csv" --inputs "$C"` (verifies SHA-256); check the provenance tables for the auxiliary products. Archive revision changes require a fresh comparison. The main reducer's binary intermediate maps regenerate locally. Published priors and program search coverage are in `publication.md`.
