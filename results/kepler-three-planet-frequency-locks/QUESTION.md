# How Kepler-359 was selected

We tested whether a slow combination of three transit frequencies reveals a persistent orbital resonance. We first checked the main-branch result directories and open PRs. Their planet work covered transit crossings and individual timing targets. We found no equivalent three-frequency search.

An exploratory look at NASA's composite planetary parameters flagged Kepler-1530. Its periods vary enough across published solutions that the putative long clock was sensitive to the choice of catalog. We therefore repeated the selection using the periods in a single published table. [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/) · [Lissauer et al. 2024, Table 1](https://doi.org/10.3847/PSJ/ad0e6e)

We used the positive-period entries with the authors' P disposition. For each three-planet combination around the same catalog target, both period ratios had to be at most four. For inner and outer positive integers a and b up to six, with a+b at most ten and no common factor, we tested frequency combinations with coefficients (a,−a−b,b). We sorted by the absolute frequency offset multiplied by the middle planet's period and divided by a+b.

The resulting search covered **201 host stars, 500 triplets and 10,500 combinations**. Kepler-60 ranked first and served as the known resonant control. Kepler-359 b/c/d ranked 19th. Its nearly century-long nominal circulation period and available joint transit-timing posterior made a direct test possible. The interactive selection and availability of the posterior preclude treating this rank as a statistical discovery probability. [Candidate ranking through Kepler-359](period_screen.csv)

Run `python selection.py` to repeat the catalog screen. It verifies the original source checksum. The [timing and dynamical analysis](README.md) tests the candidate against individual transit measurements and gives the replay instructions for the N-body experiment.
