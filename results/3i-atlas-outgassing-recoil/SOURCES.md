# Source ledger, prior overlap and physical interpretation

This is a new **gas-recoil versus independently refitted RTN orbit** investigation, not an independent new detection of 3I or an original interstellar-object discovery. The date-resolved CSV is [production.csv](production.csv); full runner network receipts, MPC 80-column rows and JPL current model parameters are committed under `generated/` with digests. Dates are UT and `Q` units are molecules or atoms per second, except the separate Hubble dust estimates.

## Primary measured gas and dust sources

| Facility and epoch | Input and interpretative caveat |
| --- | --- |
| [JWST/NIRSpec 2025 Aug 6](https://doi.org/10.3847/2041-8213/ae2fff) | Cordiner et al. total/terminal CO2 1.70e27, CO 3.7e26, H2O 2.23e26 s^-1; near-nucleus extraction and radial Q-curve smaller. The extended H2O component is likely released in coma grains. |
| [SPHEREx/IRTF 2025 Aug 8-12](https://arxiv.org/abs/2508.15469) | Lisse et al. detected CO2 9.4e26; 3-sigma **upper limits** H2O 1.5e26 and CO 2.8e26, not detections. |
| [Swift/UVOT 2025 Jul 31-Aug 1 and Aug 18-20](https://arxiv.org/abs/2508.04675) | Xing et al. OH-based water rates are reddening- and aperture-sensitive; published analyses of the same epochs vary (e.g., [Combi et al.](https://doi.org/10.3847/2041-8213/ae3bd8)). CSV retains original Jul model and later-compiled Aug rate, explicitly labeled. |
| [Lowell Observatory 2025 Aug 19](https://www.astronomerstelegram.org/?read=17352) | D. Schleicher Q(CN)~7.23e24 with asymmetric ~1.1-1.3e24 uncertainty in a large aperture; CN traces composition, not the total gas momentum. |
| [VLT/X-shooter and UVES 2025 July-Sep](https://www.research.ed.ac.uk/en/publications/very-large-telescope-observations-of-interstellar-comet-3iatlas-i/) | Ni and CN and [homogeneous OH/CN monitoring through Feb 2026](https://doi.org/10.1051/0004-6361/202662346). Trace Ni/CN excluded from bulk recoil. Telescope-dependent outflow prescriptions cause nontrivial Q changes. |
| [JUICE/MAJIS 2025 Nov 2, 12, 19, 25](https://doi.org/10.3847/2041-8213/ae8a47) | Bockelée-Morvan et al. near-perihelion H2O and CO2 rates; Nov 25 H2O marginal. Low gas expansion speeds suggested but not directly measured as a heliocentric series. |
| [SOHO/SWAN 2025 Nov onward](https://doi.org/10.3847/2041-8213/ae3bd8) | Combi et al. full-coma H-Lyman-alpha water. Large-field Nov 6 rate reported as 3.17e29 s^-1 in [independent IR spectroscopy](https://www.nature.com/articles/s41550-026-02850-5); aperture/source differences with JUICE are substantial. |
| [TRAPPIST-North 2025 Nov 29 and Dec 7](https://www.astronomerstelegram.org/?read=17515) | Jehin et al. OH and CN, using a Haser 1 km/s coma velocity; [Dec 7 follow-up](https://www.astronomerstelegram.org/?read=17538). In calculations OH is **approximate water proxy** Q(H2O)=1.1Q(OH) — not a direct H2O measurement. |
| [SPHEREx 2025 Dec 8-15](https://doi.org/10.3847/1538-3881/ae578d) | H2O ~(1.4±0.28)e28, CO2/H2O~0.21 and CO/H2O~0.55; CSV ratio-derived CO2/CO entries are tagged as derived. |
| [JWST/NIRSpec 2025 Dec 22-23](https://arxiv.org/abs/2603.20460) | Roth et al. CO 3.82e27, CO2 1.59e27, H2O 1.55e27 s^-1; different activity mix on outbound arc. |

HST constraints that set the size grid (not independent gas-production measurements):

- [Jewitt et al. HST July 2025](https://arxiv.org/abs/2508.02934): red-albedo-0.04 effective **upper radius 2.8 km**; model-dependent lower radius from sustaining dust coma ~0.16 km (CO). A revised presentation gives dust production roughly 12–120 kg/s depending on grain size ([paper](https://faculty.epss.ucla.edu/~jewitt/papers/2025/JHM25.pdf)); dust is entrained by gas and **must not** be added twice as separate independent thrust.
- [Hui et al. HST Dec 2025–Jan 2026](https://doi.org/10.3847/2041-8213/ae471c): photometric effective radius **1.3±0.2 km** for assumed pV=0.04; pV×cross-section=0.22±0.07 km². This is model-dependent and in tension with some gas-recoil interpretations, not a hard known nucleus mass.
- [BHTOM 2025 Aug–Sep photometry](https://doi.org/10.1051/0004-6361/202659913): model-dependent **dust mass-loss upper limits** ~217–328 kg/s for 100-µm grains at 5 m/s. This is not directly measured gas nor a measured recoil force.

## Prior papers already making the orbit-to-gas comparison

- [Spada, Królikowska & Dones (2026), A&A 714 A79](https://doi.org/10.1051/0004-6361/202659645): independent MPC fit with explicit astrometric weighting and outlier tests; NG detected; pre/post arcs incompatible in coefficients; JPL #54 3-axis + delayed outgassing law broadly reproduced. Their solution makes clear formal uncertainties understate systematics.
- [Thoss, Loeb & Burkert (2026)](https://arxiv.org/abs/2603.15735): CO2-dominated recoil gives R~0.42 km at rho=500kg/m³ and eta=0.5 in their assumptions; compares with Hubble photometric radius and distributed coma water.
- [Ahuja et al. (2026)](https://doi.org/10.1093/mnras/stag1800): uses published JPL A1/A2/A3 and gas production to estimate R~0.60±0.04 km at rho~538kg/m³; not an independently fitted orbit.
- [Hui et al. (2026)](https://doi.org/10.3847/2041-8213/ae471c): HST nucleus extraction and gas-recoil size estimate; uses an earlier JPL orbit release with numerically different coefficients.
- [Spada et al. extended-discussion and newer long-arc analysis](https://arxiv.org/abs/2610.38485): distant motion detected but sensitivity to activity law and astrometric systematics remains.
- [Micheli et al. (2018)](https://www.nature.com/articles/s41586-018-0254-4): 1I/'Oumuamua's nongravitational acceleration without a measured gas coma motivates the comparison, not a benchmark of 3I's physical behavior.

## Interpretation safeguards

A 3-vector outgassing acceleration equals the **vector sum** of all nucleus-coupled gas/dust momentum, divided by nucleus mass; spectra measure scalar molecular production rates over finite apertures. The recoil-grid quantity is `eta_required=M|a|/(dotM v)` and is a **conditional sufficiency check**, not a jet-direction fit, and missing molecules cannot be set to zero. Values based on upper limits, aperture-integrated coma ice sublimation and line-of-sight dust all tend to **overstate** nucleus-coupled thrust. No single epoch's failed scenario excludes all plausible nuclei or gas flows. No interval around perihelion has complete contemporaneous direct H2O, CO2, CO *and* 3D exhaust-vector constraints.

The independent astrometric runner deliberately does not implement JPL's delayed-peak DT, stellar aberration, Gaia catalog de-biasing or a nongravitational mass-loss law tied continuously to the measured Q time series. Those differences may dominate a raw formal-sigma comparison.
