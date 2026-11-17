# Secondary continuous-period model

This model was specified after the first Sector 109 measurements and partial extension results. It is a secondary analysis. Use the union of all distinct crowding-qualified Sector 109 and extension white dwarfs, one FFI visit per star. Give every host one circular, opaque companion with radius 11.2 Earth radii and a period drawn uniformly in log P between 0.25 and 3 days. Use the catalogue stellar parameters, the nominal aperture fraction and the existing blind detector. The short-period end assumes that the body can exist on the specified orbit. Companion density and dynamical survival are outside this model.

For isotropic orbits the transit probability is proportional to P^(−2/3). Integrate it analytically across the intrinsic log-uniform distribution for each star. In 50 real-flux injections per star, draw the period from the transit-conditioned density proportional to P^(−5/3), draw uniform transit impact and orbital phase, and require the original detector's one-for-one period and endpoint phase match. Use seed 909563 and the same exposure integration. Count a star once as a possible host if it has any original periodic trigger in the full period interval, with a conservative timing buffer at both edges. Keep triggers of every prior classification in that count. Use the original finite-trial confidence allocation. Report the transiting-host ceiling, the isotropic-host ceiling if the data constrain it, and the count of distinct compatible stars.

## Recovered signals and host bounds

The union contains 42 distinct white dwarfs. The detector recovered 2,035 of 2,100 injections. The two Gaia IDs 4672306015773211008 and 5063539946887524864 had one or more compatible blind periods, so both count as possible hosts for the upper limit regardless of their classification. The pointwise 95% ceiling on the fraction hosting an *already-transiting* companion under the specified period and radius model is 20.13%. With isotropic orientations the sum of the estimated detection probabilities is 1.68 effective stars, and even 100% intrinsic occurrence remains compatible. [Per-star outcomes](extension/population_stars.json) · [Individual blind matches](extension/population_trials.csv.gz) · [Confidence calculation](extension/population_limit.json)

The audit recomputes all 2,100 random periods, phases, impact parameters and transit durations. It checks the matched ephemerides at both measured endpoints, reconstructs the compatible-host count from the unmodified detector outputs, integrates each transit probability numerically, and verifies the survey confidence calculation. From the result directory run

```bash
python audit_population.py extension record/s109_ffi_visits.csv extension/extension_visits.csv --survey-record .
```
