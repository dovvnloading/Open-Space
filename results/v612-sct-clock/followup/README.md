# V612 Scuti deeper follow-up: harmonic, independent reduction, matched-exposure tests

This is **new work added to PR #24**, not a revision of its original 16 files. All times refer to BJD(TDB) where properly barycentrically corrected. Numbers are a **retrospective targeted analysis** and do not establish global discovery false-alarm probabilities or an orbital period.

## A newly detected first harmonic, recovered from raw TESS images

The baseline ZTF period is 0.15956200155 day = 3.829488 hours; the additional finding is a coherent **2f Fourier harmonic at 1.914744 hours**. This component was not included in the original PR. We independently processed the archived *raw* 5×5 TESSCut FFIs in sectors 80 and 92 (different camera/CCD geometry), without using the original PCA or rolling-median filter. We summed four Gaia-centered pixels, subtracted four times the same-cadence **mean of 11 outer pixels**, fitted a separate linear intercept/slope for each observing day, then fitted f/2, f, 2f and 3f simultaneously.

| Raw TESS FFI test | Sector 80 (2024) | Sector 92 (2025) |
|---|---:|---:|
| Valid exposures | 7,563 | 9,175 |
| Occupied days | 22 | 25 |
| f semiamplitude (e−/s) | 2.903 | 3.192 |
| f phase relative to ZTF maximum, cycles | −0.0053 | −0.0099 |
| **2f semiamplitude (e−/s)** | **0.653** | **0.528** |
| 2f 95% day-resampling amplitude distribution | 0.458–0.865 | 0.371–0.729 |
| **2f phase, cycles at 2f** | **−0.4645** | **−0.4780** |
| Day-block sign-flip null (10k trials) | 0 exceedances | 3 exceedances |

The 2f phases disagree by just **0.0135 cycle**. Fixing the 2f phase from sector 80, one-sided sector 92 testing found **0/10,000** flipped-day null realizations as aligned as observed. These probabilities condition on a previously selected star and frequency, extraction, and daily-noise model. The two independent observational seasons were investigated retrospectively, not under a true blind preregistered analysis. The original PR PCA-corrected files give 2f amplitudes 0.641 and 0.512 e−/s with respective phases −0.466 and −0.453; the independent linear analysis therefore does not depend on a PCA artifact.

The source's dominant f waveform in the red TESS band is **not a pure sinusoid**: 2f/f amplitude ratios are approximately **0.225 and 0.165**. The near-antiphased second harmonic gives broadened optical brightness maxima and sharper minima. This can arise from nonlinear orbital geometry or irradiation/eclipses, but the data **do not demonstrate an eclipse**; 3f is not securely detected. The harmonic's source position from simple Gaussian pixel maps is less secure than the f localization, and a calibrated TESS PRF fit would be preferable.

## Stricter ZTF ground-based systematics control and temporal holdout

The original PR tested a nearby comparison star independently. Here we **pair 246 exposures by their exact exposure IDs**, subtract the simultaneous nearby g-band comparison star (OID 333114400172086) from the nova-position source (OID 333114400007906), and detrend the differential magnitude with a ±25-day median. This removes same-exposure zeropoint variations rather than relying on separate unmatched light curves.

At the already-frozen frequency, the differential series has **0.04715-mag semiamplitude** and a **magnitude phase of +0.4941 cycles** against the original brightness maximum (magnitude minima should occur at ±0.5). We fit the sinusoid coefficients on the first **152** selected observations, before HJD 2460300, and predicted the remaining **94** 2024–25 observations without refitting amplitude or phase. The prediction has **Pearson r = 0.5653** and a conditional one-sided **14-month-block sign-flip p = 0.00034** (33/100,000 randomized trials at least as aligned). This holdout split was chosen after identifying the source, so it is an internally held-out *coefficient test*, not a blind independent survey discovery.

A simultaneous target-minus-control astrometry fit using 150-day nuisance offsets returns a **first-moment-assumption** variable-component displacement **(+0.019, −0.153) arcsec**, east/north, with marginal 95% 30-day bootstrap intervals **[−0.117,+0.153]** and **[−0.284,+0.013] arcsec**. Gaia's nearest separate star is **2.71 arcsec** away. This strongly disfavors that neighbor under an *ideal first-moment blend model*, **but ZTF's PSF-fit catalog positions do not obey that model exactly**, and the interval does not include calibration of astrometric systematics. Neither a subarcsecond unresolved neighbor nor pipeline PSF bias is definitively excluded.

Seven yearly paired-exposure records can be regenerated with the additive reanalysis script. The high-S/N 2022–25 optical modulation retains its phase while its amplitude changes. This is new independent-exposure subtraction, not a new target-selection event.

## Orbital period: both 3.8295 h and 7.6590 h are still viable

| Interpretation | Model | What the evidence says |
|---|---|---|
| **P_orb = 3.829488 h** | One brightness variation per orbital cycle (accretion geometry, irradiation, eclipse, beaming), with a nonlinear 2f harmonic | Compatible; **not measured by radial velocities** |
| **P_orb = 7.658976 h** | Nearly equal double-wave ellipsoidal variations at twice the orbital frequency; observed 2f is then fourth orbital harmonic | Also compatible; perfectly symmetric double-wave need not show f/2 |
| White-dwarf spin, nonorbital beat or superhump | An optical clock not exactly at the orbital period | Cannot be eliminated without spectroscopy / polarized or X-ray timing |

The *subharmonic* f/2, corresponding to a 7.658976-hour odd-even modulation, is **not robust across observations or extraction choices**. In a daily-linear FFI reduction f/2 semiamplitude is **0.219** e−/s (sector 80, day-sign p≈0.45) and **0.830** e−/s (sector 92, p≈0.054), but the respective phases are **−0.385 and −0.059 cycle** and do not form a coherent two-sector alternating signal. Bootstrap amplitude distributions reach 0.60 and 2.04 e−/s; these are **not rigorous intrinsic-amplitude upper limits**. No f/2 detection is not evidence that P_orb must be 3.8295 h.

## Attempted third external photometric validation: inconclusive, not fabricated

MAST's sky-location sector query confirms **only TESS sectors 80 and 92**. ASAS-SN Sky Patrol contains 1,117 measurements toward the proposed Gaia counterpart. Restricting to **711 accepted post-eruption g measurements**, the free-phase fit at the ZTF frequency has **0.0108-mag semiamplitude**, against **0.187-mag residual rms**; month-block bootstrap of fitted amplitude spans about 0.003–0.040 mag, and the ZTF-predicted *sign* has **30-day block sign-flip p≈0.23**. This is **not a third survey detection**. The three nearest ASAS-SN master-list entries have paired-epoch g magnitudes correlated **0.995–0.998** because their apertures cover the same crowded region; they are not independent neighboring-star controls. ASAS-SN provides unbarycentered JD timestamps, so a precise absolute phase test also requires observatory correction absent from those rows.

Gaia DR3 source **4104113350446549888** has the metadata flag **has_epoch_photometry=false**, so no Gaia DR3 epoch curve is available to phase-test. The original Pan-STARRS 69 detections (27 g, 10 r, 13 i, 14 z, 5 y after cuts) predate the 2017 nova and cannot test late-decline phase coherence or prove pre-eruption presence. ATLAS forced photometry was not accessed via an authenticated account; no ATLAS detection is claimed. Additional catalogs are not replaced with inferred data.

The image-level ZTF difference-imaging avenue is also incomplete: of 120 phase-stratified public original-image cutout requests, **112 returned IRSA HTTP 503** and only eight FITS cutouts were retrieved. Eight images spanning years do not support a calibrated differential-PSF centroid fit; we **do not** assert such a localization from them. Full-calibration image differencing remains a decisive next test.

## Falsifiers and decisive predictions

1. **Source association**: high-resolution (≲0.3-arcsec) time-domain PSF or difference photometry that places the f variation on a *different* Gaia object, while the candidate progenitor is nonvariable at sufficient sensitivity, would falsify the proposed nova-source identification.
2. **Real astrophysical harmonic**: recovery of the original f but **failure of 2f in calibrated independent reference-difference imaging** at a precision capable of detecting a 0.17–0.23 relative-amplitude harmonic, or convincing localization of 2f on another source, would falsify its shared astrophysical waveform interpretation.
3. **Clock stability**: a high-S/N, deblended time series several years later that rejects the original ephemeris even after its accumulated uncertainty and phase-correction errors are propagated would falsify the “phase stable across years” interpretation.
4. **Binary orbital period**: time-resolved absorption-line spectroscopy of one orbital component distinguishes the two periods without relying on the absolute optical phase. If P_orb=3.829488 hours, **v(t+3.829488 h)=v(t)**. If P_orb=7.658976 hours, **v(t+3.829488 h)=2γ−v(t)**, where γ is the systemic velocity. The full 7.658976-hour cycle restores the velocity. In either case a velocity line that traces only variable disk emission is insufficiently reliable, and lack of absorption does not prove a spin clock.
5. **Unequal eclipses/ellipsoidal parity**: ≲2-minute cadence, high-S/N simultaneous two-band photometry over at least several dozen optical cycles should explicitly fit f/2 and even/odd minima. A coherent f/2 difference across seasons favors 7.659 h; **no f/2 cannot reject perfectly symmetric ellipsoidal modulation**.

**Requested observing program:** at least **16 hours total** covering two candidate 7.66-hour cycles across suitable nights, resolving donor absorption features at approximately R ≥ 5,000–10,000 with **≤10–15-minute exposures** and enough signal-to-noise for line cross-correlation; simultaneous photometry to update phase. Period-vs-period RV model comparison must allow systemic velocity and phases to float. Multi-line donor absorption, where detectable, is preferable to accretion-line emission.

The original **brightness** ephemeris is T_max(BJD_TDB) = 2459499.94732823 + 0.15956200155 E, integer E. As a concrete test, **E=12,896** predicts maximum at **BJD(TDB) 2461557.658900** near 1 June 2027. Propagating only the original conditional 95% period range gives roughly **−30/+25 minutes** of timing uncertainty by then, before adding epoch or systematic errors. It would be misleading to supply a seconds-accurate 2027 UTC eclipse schedule. Re-phase with contemporaneous photometry.

## How to reproduce and audit

Run with NumPy, pandas, SciPy and Astropy; all new scripts read the original PR files without rewriting them.

- Exposure-matched, time-held-out ZTF analysis and original archived TESS harmonic checks: run **python results/v612-sct-clock/followup/reanalysis.py --check --output /tmp/v612-photometry.json** from repository root.
- Independently reduce two original raw TESS FFIs: first run **python results/v612-sct-clock/followup/harvest.py** (defaults to /tmp/v612-harvest), then **python results/v612-sct-clock/followup/linear_pixels.py --fits-dir /tmp/v612-harvest --check --output /tmp/v612-linear.json**.
- ASAS-SN supplemental test: add **--asassn-dir /tmp/v612-harvest** to the photometric replay.
- The new [measurement ledger](new-results.json) captures all claim-relevant values, source checksums and time-window limitations. An automated [GitHub Actions raw-pixel artifact](https://github.com/dovvnloading/Open-Space/actions/runs/43668969680) is available under GitHub's retention policy.

**Claim boundary:** new replicated **2f waveform harmonic** and stricter matched-exposure light/centroid analysis; *not* a spectroscopic binary period, *not* a third-survey confirmation, *not* calibrated source separation at 2.7 arcsec, and *not* a global false discovery probability.
