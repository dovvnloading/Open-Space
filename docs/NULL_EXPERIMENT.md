# How much recurrence survives a redshift shuffle?

On the first SDSS plate visit, one same-sign spectral feature group contained eight of twelve galaxies, and 25 groups contained at least four. We held the measured wavelengths and the 12 measured redshifts fixed, then reassigned an entire redshift to each spectrum. Across 9,999 unrestricted shuffles the largest group reached six galaxies. Restricting exchanges to galaxies with neighboring redshifts brought the observed result within reach of eight of 1,296 possible assignments. This experiment uses only the original observing visit. [Original data and selection](../results/first-visit.json) · [Full results, input hashes and histograms](../results/redshift-null.json)

| Statistic at 150 km/s | Measured | Unrestricted shuffles, median / 95th percentile | Redshift-neighbor assignments, median / 95th percentile |
| --- | ---: | ---: | ---: |
| Size of largest group in distinct galaxies, primary | 8 | 3 / 4 | 4 / 6 |
| Number of groups containing at least four galaxies, secondary | 25 | 0 / 1 | 4 / 16 |
| Number of groups containing at least two galaxies, descriptive | 47 | 32 / 48 | 44 / 54 |

Each count uses the original algorithm for selecting nonoverlapping groups separately in emission and absorption. The table includes all 333 detections before catalog screening. Percentiles use the nearest-rank convention. The total of 47 groups was reached or exceeded in 635 of 9,999 unrestricted draws. A group shared by eight targets tests a different property of the data from a count dominated by two-target coincidences.

## What is exchanged

1. Use one completed, compatible detection report per distinct sky position on the same plate and MJD. Keep galaxies with zero unassigned features, since their redshifts are part of the assignment pool. Each galaxy retains its observed wavelengths, signs, scores, pixel spacings and feature identifiers.
2. Apply one reassigned redshift to every detected feature in that spectrum. Compute its proposed rest wavelength as observed vacuum wavelength divided by `1 + assigned redshift`. Reapply the same cosmic and airglow line catalogs and matching tolerances used in the discovery scan. Re-running the original association reproduces all 47 group memberships exactly. All simulated outputs stay separate from the evidence packets.
3. Form groups using the original maximum span of 150 km/s, one vote per distinct galaxy, sign separation, scoring and overlap resolution. Record the strongest multiplicity anywhere across both signs and the count of groups supported by four or more targets. We selected these global summaries after the discovery results were known and before drawing the random assignments. The secondary statistic supplies an additional description of this sample.
4. Draw 9,999 independent assignments uniformly with replacement from the 12! possible index assignments. Fixed points and the original assignment are permitted. Separately sort galaxies by measured redshift, split into four consecutive blocks of three, and enumerate every within-block permutation. This retains redshift neighborhoods and permits all `(3!)^4 = 1,296` assignments, including the original one. The software records the ordered block membership and a digest of each assignment sequence.

For each of the two named statistics, the unrestricted run gave zero shuffles reaching the measured value. With the conventional finite-simulation correction, its upper-tail fraction is `(0 + 1) / (9999 + 1) = 0.0001`. The restricted calculation gave the exact fraction `8 / 1296 = 0.00617` for each. Equal fractions do not imply independent evidence. Their interpretation as randomization probabilities depends on the posited exchangeability of the galaxy redshifts. [Phipson & Smyth on Monte Carlo permutation probabilities](https://doi.org/10.2202/1544-6115.1585)

A second unrestricted calculation retains the 241 detections that originally missed the catalog, even if a proposed redshift would change a match. This gives the same maximum group range of two to six galaxies and zero draws reaching either high-multiplicity statistic. The primary calculation repeats the catalog screen because eligibility for a rest-frame line match depends on the assigned redshift. It yields a median of 313 unassigned detections after shuffling, compared with the observed 241. The additional detections include physically familiar lines displaced from their catalog wavelengths by the shuffle.

## Physical interpretation

The measured redshifts organize the spectra in a way a random reassignment rarely reproduces. Ordinary host-galaxy absorption, blends and features missing from the small source catalog can all generate this pattern. For example, the original calcium investigation found three centroids consistent with Ca II K. The SDSS redshift fit itself uses spectral structure from the same coadds used in the recurrence screen. [Ca II K investigation](FIRST_INVESTIGATION.md) · [SDSS redshift estimation](https://www.sdss4.org/dr17/algorithms/redshifts/)

All twelve spectra share a plate, night and reduction. Their searchable center counts span 3,274–4,357 and their total detected features span 6–72. Each shuffle inherits all gaps and noise effects present in the detected observer-frame features. Assigning a different redshift moves those gaps and the sensitivity function to different rest wavelengths. Luminosity selection, spectral class, line strength and coverage may vary with actual redshift. Restricting the shuffles to redshift neighbors addresses only part of that dependence. No flux is simulated or remeasured. Blank-sky fibers, exposure differences and held-out observations test other mechanisms. These results supply neither a false-discovery rate for the catalog nor a probability for any individual transition. [SDSS spectroscopic caveats](https://www.sdss4.org/dr17/spectro/caveats/)

## Reproduce

Use the current checkout to analyze an original campaign directory containing its twelve completed reports. The command runs offline using the Python standard library.

```sh
python -m open_space.redshift_null runs/replicate-first-visit \
  --draws 9999 --seed 6497 --tolerance-kms 150 > runs/redshift-null.json
```

To regenerate the original evidence from SDSS, use the exact code revision recorded as `software_commit` in the [first-visit manifest](../results/first-visit.json) and the [full visit recipe](../results/README.md). For example, run these commands from the current repository checkout and use a local POSIX filesystem for both directories.

```sh
git worktree add ../open-space-first-visit 3b04f60cda92dc7add47ea2e4ba171c67a332f58
(cd ../open-space-first-visit && python -m venv .venv && .venv/bin/pip install -e '.[archive]' && \
  .venv/bin/python -m open_space.campaign --root runs/replicate-first-visit \
    --plate 6497 --mjd 56329 --limit 3 --pages 5 --max-calls 0)
mkdir -p runs
python -m open_space.redshift_null ../open-space-first-visit/runs/replicate-first-visit \
  --draws 9999 --seed 6497 > runs/redshift-null.json
```

The output records the input source, derivative and feature-vector digests, the detector fingerprint, the null-analysis implementation hashes, the random seed, the exact assignment sequence hashes and the full histogram of each statistic. Sorting observation IDs before generating assignments makes the result independent of report traversal order.
