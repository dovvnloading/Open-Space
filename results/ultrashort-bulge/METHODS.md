# Reproduction and definitions

Start at the repository root. The scripts import [lens.py](../ogle-dark-lenses/lens.py) for OGLE cleaning, magnitude conversion and point-lens magnification, and [galaxy.py](../ogle-dark-lenses/galaxy.py) for Galactic density, velocity and quantiles. Python, NumPy, SciPy, pandas, Astropy, pyerfa, lxml, Matplotlib, Poppler, VBMicrolensing 5.6.1 and a C++17 compiler were used. Preserve independent caches for different dates. Live alert files can change. The survey manifests give the upstream SHA-256 hashes for comparison.

## Main survey census

From the repository root, the commands below rebuild the common screen. The directories under `staging` hold per-event checkpoints. On a clean machine the downloads and computations can take hours. KMT's release status can be checked in the annual [2026 event index](https://kmtnet.kasi.re.kr/ulens/event/2026/). This snapshot's KMT numerical inputs cover 2024 and 2025.

```sh
export OPENBLAS_NUM_THREADS=1
r=results/ultrashort-bulge
mkdir -p staging
g++ -O3 -std=c++17 -fPIC -shared "$r/grid.cpp" -o "$r/grid.so"
python "$r/ogle/freeze.py" .ogle-cache/raw --jobs 6
python "$r/kmt/freeze.py" --cache staging/kmt --download --jobs 4
python "$r/kmt/fallback.py" --cache staging/kmt
python "$r/kmt/manifest.py" --cache staging/kmt
python "$r/other/fetch.py" --cache .survey-cache --workers 4
python "$r/other/gaia.py"
python "$r/run_grid.py" --checkpoint staging/screen --workers 3
python "$r/run_grid.py" --kmt-bound --checkpoint staging/bound --workers 3
python "$r/cadence.py" --out staging/all-cadence.csv --workers 3
python "$r/consolidate.py"
python "$r/associations.py"
python "$r/kmt/field_coverage.py"
```

Each event ID is the unit of processing. Official cross-identifications are preserved separately from positional neighbors. The model uses all valid epochs for each stream. OGLE uses HJD and magnitudes, KMT uses its HJD and differential flux, MOA/PRIME use the native JD convention and difference flux, and Gaia uses JD(TCB) and G magnitudes. The merged `first_hjd` and `last_hjd` columns are native time values. See `time_system` before comparing epochs between surveys. Gaia supplies no per-transit photometric uncertainty, so the analysis estimates each source's uncertainty from its baseline. Some archives span years before the discovery season.

The point-source and uniform finite-source point-lens models share an event time, impact parameter and Einstein time across streams. Finite-source fits add the angular source-to-Einstein-radius ratio ρ. Each stream has an independent additive baseline and nonnegative source flux. The primary difference-flux screen permits arbitrary positive amplitude. For KMT, a second screen caps the source at 1.15 times the reported total baseline light where its magnitude-to-flux calibration succeeds. Some streams and alternate DIAPL-only reductions have no comparable cap. Total-flux survey streams use their measured median baseline. Difference-flux data have an error floor at 0.15 times the median quoted uncertainty per stream.

The common bank uses four point-source impact parameters, ten extended-source profiles (ρ = 0.25, 0.75, 1.25, 2.5, 5) and four genuine finite-source binary-lens trajectories. Its time grid spans 0.03–250 d plus half, once and twice the catalog timescale where usable. Peak trials include the catalog epoch, four nearby offsets and data-driven peaks. This samples the whole archive per reported event, and can miss short anomalies at unsampled times. Primary and capped flags include any fitted point or finite-source time <0.5 d or finite-source ρ>1, without a significance requirement.

## Follow-up and validation

The KMT audit refits each flagged capped geometry by detector field and site. It measures rise/decay sampling, seeing dependence, off-peak excursions and V-band support. Qualifying short fits or excess finite-source gains trigger continuous per-stream fits. A second selection revisits free fits near the half-day boundary when the flux cap would remove them. This recovers a 0.451 d fit for published KMT-2024-BLG-0816, whose capped grid time is exactly 0.50 d. The public per-alert rows retain cases below the follow-up thresholds.

```sh
python "$r/kmt/audit.py" --cache staging/kmt --grid staging/bound
python "$r/queue_refits.py" | xargs -r python "$r/refine.py" --bound --grid staging/bound --out staging/continuous
python "$r/queue_refits.py" --crossings | xargs -r python "$r/refine.py" --bound --grid staging/screen --out staging/crossings
python "$r/compare_leads.py" --workers 3
python "$r/alternative_checks.py" --workers 3
python "$r/inject.py" --per-survey 24 --seasons 12 --out staging/injections-final.csv
python "$r/summarize_leads.py" --injections staging/injections-final.csv
python "$r/earth_control.py" "$r/controls/earth-figure1.tsv.gz"
python "$r/mass_controls.py"
python "$r/conditional_geometry.py"
python "$r/validate.py" --models
```

For manually explored peak fits, use the same `window_halfwidth_d`, `bin_width_d`, reference epoch and source cap as the appropriate [stored record](tables/peak-comparisons.json.gz). Earlier individual files without those window fields use ±4 d and native epochs. The full-curve and peak-window fits use different baseline summaries and error inflations. Compare chi-square only across models evaluated in the same experiment. `peak_models.py` uses separate per-stream baselines and compares point lens, finite-source point lens, an asymmetric exponential with free amplitude, and multiple binary-lens starts. The binary search samples a limited range of mass ratios, separations, orientations and starting geometries. It can reject a proposed isolated interpretation when another model predicts the data comparably well. The finite grids and local optima cannot exclude arbitrary hosts.

In [peak-comparisons.csv](tables/peak-comparisons.csv), the reference single lens is the lowest chi-square of the unrestricted fit, restricted restarts and point-source boundary. The four profile costs restrict tE ≥ 0.5, ρ ≤ 1, tE ≤ 0.5 and ρ ≥ 1, respectively. Blank costs denote missing profiles. A positive binary-minus-single value is an optimizer result for the sampled companion models. A binary with negligible effect can approach the single lens. Photometric errors can be correlated within nights and detector fields. The tabulated costs are fit comparisons, without a conversion to Gaussian detection significance. [Non-KMT procedures](other/README.md) and [OGLE-specific procedures](ogle/README.md) give their own thresholds, reduction conventions and input records.

## Cadence, injections and population limits

The [exact window table](tables/cadence.csv.gz) integrates the portion of 1 March through 1 November that has two distinct observations on each side of a trial peak within a specified half-width. Records before or after these dates can support a peak close to a season boundary. Each field has its actual interruptions and nightly cadence. [OGLE's separate five-exposure calculation](ogle/window_by_field.csv) uses 1 February to 1 November and includes source brightness, baseline scatter and blending.

The [injection experiment](tables/injection-summary.csv) uses the lowest 24 deterministically hashed, available numerical arrays per survey, including both released KMT years. Each curve has 12 off-event trial times uniform in the March–November window after removing a neighborhood of the catalog peak. The simulations use u₀ = 0.35, tE = 0.03, 0.10 or 0.30 d, and ρ = 0 or 2. An injection peaks at 12 times the median off-event *reported* error in each stream. Its paired control is the same real photometry without injected flux. Recovery requires two distinct exposures on each side inside tE × max(1,ρ), a localized short or large-ρ single-lens peak competitive with the trial binary shapes, and at least 80 units of gain above the control's best single lens. Off-peak baseline information is included per stream. The search is supplied the trial peak, which makes it easier than a blind alert search. Trials at different geometries share the same observing windows. The trial fractions characterize this experiment on previously alerted sightlines and have no monitored-star-time denominator.

## Physical interpretation

A measured Einstein radius needs both the lensed source's angular radius and ρ. For a finite source, θE = θ★/ρ and μrel = θE/tE. Mass also depends on lens distance through θE² = κ M πrel, unless parallax or another mass constraint is measured. The [conditional source hypotheses](kmt/extended/angular.json) propagate candidate IR source identity, extinction and the surface-brightness relation. The [mass sensitivity table](tables/conditional-source-masses.csv) gives conditional 16th, 50th and 84th mass percentiles separately for each assigned star and each geometric model. These rows use fixed angular size and timescale within a case, two slopes for dN/dlogM (0 and −0.96), and a mass range of 0.0001–1,000,000 Earth masses.

The Galactic calculation assumes bar sources at 5–12 kpc and disk/bar lenses from 0.15 kpc up to the source. It uses the repository's disk/bar density and kinematics, corrects to the geocentric frame, and weights by microlensing event rate. At fixed θE and tE, the integration samples direction at the required proper-motion magnitude in the two-dimensional velocity density. It remains numerically stable for very slow projected motions. This calculation assigns no probability to the star identification, the single-lens hypothesis, a stellar host, or free-floating status. The [published-event validation](controls/conditional-masses.json) propagates reported measurement errors and uses its stated 0.03–3,300 Earth-mass prior range.
