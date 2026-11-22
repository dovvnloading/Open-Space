# Survey definition and reproduction

## Archive population

The unit of selection is a Chandra ObsID. The CXC TAP inventory uses the public availability date, inclusive from 2021-01-01 to 2026-11-21. An older pointing can enter the survey when its proprietary period ends. The query reads disjoint numeric ObsID ranges and reconciles both native rows and rows joined to proposal metadata against independent counts. When an ObsID has multiple archived imaging records, its earliest known public release controls inclusion. Gratings and continuous-clock exposures fall outside this two-dimensional imaging search. Missing release dates, ambiguous configurations, inaccessible files, failed scans and partial follow-up have separate receipts. Target identity and Galactic latitude play no part in the event search. The inventory preserves target and proposal information for reviewing the selection of any surviving source. [CXC TAP schema](https://cda.harvard.edu/cxctap/tables)

Each observation uses the latest archive Level 2 event file for each distinct imaging stream at inventory execution. The worker records the source URL, archive filename, FITS revision, byte length and SHA-256 of both the transport object and the decompressed FITS bytes. A repeat dispatch uses the original pinned input after acquisition. Exact-byte replay fails visibly when an archive revision is unavailable. Unclassified surviving real events also cause the workflow to retain the original EVT2 as a separate run artifact for up to 90 days.

## Trigger construction and vetoes

The search keeps 0.5–7 keV, zero-status events in ACIS grades 0, 2, 3, 4 and 6. On each CCD and uninterrupted good-time interval it scans 8- and 24-pixel sky cells in four offset tilings. At each photon it tries windows of 60, 180, 600, 1,800 and 3,600 seconds. A seed requires eight photons, four times the average cell rate and a local binomial upper tail at most 1e-7. The local statistic ranks correlated trials. The survey-wide significance estimate uses the scrambles described below.

A source aperture has radius 0.7 times the cell width. A concentric annulus runs from two to four aperture radii. Each trial requires at least one trial-length of contiguous coverage on each side. The code tests the flank rates, their adjacent counts, local spatial contrast, annulus and chip-wide background variation, edge truncation, detector-pixel concentration, frame occupancy, fixed columns, S4 streaks, and alignment with the readout direction. It rejects subminute photon concentrations. The ledger records each deduplicated trigger, its measured diagnostics and the veto reasons. Data-driven instrumental cuts have limited power for sparse or unusual artifacts.

Sky coordinates come from the binary-table WCS. The reported duration is the 5th-to-95th percentile span of aperture photon times in the selected window. At low counts this interval depends strongly on individual photons and residual background. The reported peak is the brightest observed fully enclosed 60-second aperture interval. Real leads receive background-subtracted energy flux through the CIAO per-event mirror and detector response. Counting uncertainties are conditional on the selected peak. Aperture losses, absorption, redistribution, pile-up and the bias from selecting a maximum remain relevant. Surviving real events receive nine-panel diagnostics with WCS sky apertures, detector coordinates, source and background light curves, energy versus time, and an independent field-derived dither trace. The cutout records the zero-based photon rows in the original event table. The photometry step preserves individual weighted photons and the CIAO/CalDB build receipts. [CIAO event-flux method](https://cxc.cfa.harvard.edu/ciao/threads/eff2evt/)

## Catalogs, false alarms and completeness

Every real and scrambled trigger follows the same detector and catalog selection. Gaia DR3 positions propagate to the observation epoch when proper motions are available. SIMBAD and NED supply classifications and possible host galaxies. Local object densities quantify chance projections under a locally uniform surface-density model. Published Chandra variability lists flag prior detections, including earlier work overlapping this release window. Unclassified counterparts, missing catalogs, crowded fields, repeaters and a blank field leave the event's origin open. A positional galaxy redshift belongs to the galaxy. Human review of the observing proposal, discovery notices, repeat visits, spatial cuts and the exact photons is required before declaring a previously unreported extragalactic event. [Gaia DR3](https://www.cosmos.esa.int/web/gaia/dr3) · [SIMBAD](https://simbad.cds.unistra.fr/simbad/) · [NED](https://ned.ipac.caltech.edu/) · [Quirola-Vásquez et al.](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/A+A/675/A44) · [Dillmann et al.](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/MNRAS/537/931) · [Saathoff and Larsson](https://doi.org/10.1051/0004-6361/202557261)

A time-scramble permutes timestamps within each CCD and good-time interval, keeping each event's sky position, detector position and energy together. Each of 24 survey replicates uses a separate deterministic seed for each product and segment. The controls preserve chip-wide time histories, integrated images and the frame-time distribution. They destroy correlations between event time and position, including residual dither artifacts, afterglows and true variability. The reported false-alarm probability is conditional on that exchangeability assumption. The catalog identifiers and response bodies can be replayed offline inside the workflow.

For each real lead take the maximum surviving seed score over each scrambled survey. If k maxima equal or exceed the lead score, the permutation estimate is (k+1)/25. The output also gives the two-sided exact binomial interval for the unknown exceedance probability. The smallest reportable estimate with this design is 0.04. An incomplete catalog stage switches the reference to the more inclusive instrument-passing control maxima. The empirical mean number of control events and the rate per million searched seconds describe this particular mixture of fields.

Each event product receives 12 pre-seeded draws from the trial durations, 8/16/32/64 detected photons and Gaussian widths of 1.2/3/6 native sky pixels. Each top-hat in time is inserted individually into the recorded background. The simulation uses observed frame times, chip orientation, energy sampling and dither offsets. Recovery requires a spatial and temporal match, an instrument-passing trigger and at least four injected photons among its selected events. The table gives the scheduled, attempted and recovered counts and exact binomial intervals in each bin. These are product-weighted instrumental recovery fractions conditional on the simulated positions, shape and observed counts. They do not measure sensitivity to incident flux over the sky or completeness of host identification.

## Denominator and interpretation

A trial of length D in a good-time interval of length L has max(0, L−3D) seconds of possible start times owing to the quiet flanks. Intervals are weighted by the live fraction. Simultaneous CCDs and streams contribute the greatest available live fraction once within an ObsID. This is conservative for complementary interleaved frames. A partial multi-product ObsID contributes no exposure to inference. Time support does not measure effective sky area. A zero-survivor bound uses −ln(0.05) divided by the stated searched start-time exposure for the specified duration. It applies to detectable events in the surveyed field mixture. Stronger population limits require a luminosity, duration, host and survey-area model.

## Run and replay

The maintainer enabled pushes that change the Chandra workflow file on the research branch. The Actions run link in each deposited summary identifies the exact revision. After default-branch registration, the same workflow also accepts the manual commands below. All archive retrieval, catalog queries, calibration, control generation, injections and repository checks execute in that workflow. The first three moderate-volume observations are selected using event counts and exposure only. Their complete evidence is checkpointed before the remaining inventory is assigned. The primary search uses at most eight matrix workers. One independent worker examines the fourteen ambiguous release dates. Run artifacts preserve partial jobs. A subsequent manual dispatch resumes from deposited ObsID ledgers and pins previously acquired data. Repository writes are confined to the research branch.

```sh
gh workflow run chandra-acis-fxt.yml --repo dovvnloading/Open-Space --ref research/chandra-acis-fast-transients
gh run list --repo dovvnloading/Open-Space --workflow chandra-acis-fxt.yml --branch research/chandra-acis-fast-transients
```

A second checkpoint selects up to 32 ObsIDs in each of four archive-count intervals, 0–8,000, 8,001–24,000, 24,001–72,000 and 72,001–200,000. Within each interval it ranks the SHA-256 of the ObsID and fixed cohort name before opening events. This stage runs the same blind search on up to 128 observations and deposits the individual ledgers in a single batch. The full census then schedules every remaining observation, including all higher-count and unknown-count fields. A separate worker processes every CXC release-date exception still unresolved in HEASARC and records its inclusion evidence independently of the primary census.

The preflight queries Gaia, SIMBAD, NED, and three published Chandra transient lists at the deposited position of an instrument-vetoed detection in ObsID 31579. Evidence from the initial run lives in `evidence/history/44981778618/`. The catalog correction preserves its original decisions and search arrays. The migration checks the event-product hashes, control and injection indices, inventory checksum and implementation hashes before reusing a scan. Workers apply the corrected catalog checks to every archived trigger, including all scrambled triggers. The HEASARC release-date cross-check for incomplete CXC metadata has its own input and outcome receipts under `evidence/release_supplement/`.

For retained leads the report job checks other visits in the release window and searches the GCN Circular archive for potentially related notices. Inspect the source and detector diagnostics, the photon time series, catalog matches, nearby galaxies, and target/proposal metadata for each event. A non-detection in another visit constrains the mean aperture count rate only for the time its source position was covered. The alert comparison records the notices and search coverage used for a novelty assessment.

After a scientific deposit, substitute an integer ObsID from the coverage ledger for OBSID to repeat its event scan, controls, synthetic injections, catalog choices and measured peak photometry.

```sh
gh workflow run chandra-acis-fxt.yml --repo dovvnloading/Open-Space --ref research/chandra-acis-fast-transients -f replay_obsid=OBSID
```

The verification job and the evidence-commit job both require these commands to pass.

```sh
python -m open_space.verify
python -m unittest discover -s tests
python -m unittest discover -s results/chandra_acis_fxt -p 'test_*.py'
```
