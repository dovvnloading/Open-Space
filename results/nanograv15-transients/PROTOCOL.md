# Decisions made before examining residuals

## Discovery and validation samples

Select the complete narrowband timing release, including short baselines. Keep selection and fit failures in the coverage ledger. Use the wideband timing solution for the same pulsar as a reduction check. Compare fixed epochs and templates against European Pulsar Timing Array DR2 where its observations bracket the event. The two NANOGrav reductions share telescope data and cannot count as independent discoveries.

Use PINT timing predictions, the timing-parameter design matrix and multi-frequency arrival times. Refit or project timing-model errors. Represent both achromatic red fluctuations and frequency-dependent dispersion. Retain observing-system labels and timing uncertainties. Profile and dispersion-measure changes can have similar signatures in arrival times. A profile change requires direct pulse-profile evidence to identify its cause.

## Tests and error budget

Freeze the event families, epoch grids, widths and vetoes in code before opening the primary residuals. Scan positive and negative amplitudes. The family-wise false-alarm budget is 0.005 for individual pulsars and 0.005 for the array. The combined budget is 0.01 under the specified noise nulls. Keep primary threshold crossings, including failures of the noise and cross-release checks. Label angularly unresolved array responses and ambiguous single-pulsar waveforms as such.

Sensitivity is a detection probability for the specified sampling, noise model, event shape and search threshold. Record injection amplitudes, number of trials and recoveries, the sampled epoch/width distribution, and what the recovery rule tests. Report any interval of sensitivity with its units, its detection stage and its coverage. A limit on timing amplitude does not by itself give a population limit on compact-object abundance.

Record the input archive versions, checksums of downloaded bytes, per-file provenance, software versions, run revision and the commands that produced the tables. Archive selection failures, missing comparison epochs and rejected candidates. A failed fetch or failed search does not establish a nondetection.

## Publication history and interpretation

Keep published timing-event epochs out of the template grid construction and candidate ranking. After running the blind search, cross-match each survivor against prior work on that pulsar and epoch. Document the cited source, its event epoch and the overlap or mismatch. Mark unresolved literature checks explicitly. A previously reported event may validate the search. Independent detection and publication history are separate questions.

- The NANOGrav 15-year release paper describes how the arrival times and timing solutions were constructed. [Agazie et al.](https://arxiv.org/abs/2306.16217)
- The NANOGrav new-physics search tested Doppler and Shapiro signals from compact substructure. Recovered waveforms require comparison with its tested regime. [Agazie et al.](https://doi.org/10.3847/2041-8213/acdc91)
- NANOGrav has published a memory search using this release. Any comparison of strain sensitivity must specify sky position, epoch and confidence construction. [Agazie et al.](https://authors.library.caltech.edu/records/t8d0v-hvz44)
- Earlier chromatic-noise studies include the two events in J1713+0747 near MJD 54750 and 57510. [Larsen et al.](https://doi.org/10.3847/1538-4357/ad5291)
- Pulse-profile work includes known changes in J1713+0747 and J1643−1224 and additional ranked epochs. The customized noise study covers 67 pulsars. [Jacobson-Bell et al.](https://doi.org/10.3847/1538-4357/ae6db6) · [Larsen et al.](https://arxiv.org/abs/2606.28571)
