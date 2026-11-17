# DESI DR1 catalogue census

Of the 63,968 objects visually inspected by Swan et al., 63,152 have distinct Gaia DR3 designations and 816 have only a DR2 designation. All 63,152 DR3 IDs were located in Gaia. This parent includes stars, extragalactic contaminants and ambiguous spectra as well as classified white dwarfs. Keeping them preserves the chance of recovering an unusual spectrum misclassified by the programme.

The catalogue-first search has no raw speed, colour or magnitude preselection. We required a usable five- or six-parameter solution, finite two-component proper motion, an applicable individual Lindegren zero point and a corrected parallax at least five times its reported uncertainty. We then required 0 < RUWE < 1.4, at least nine visibility periods, a clear duplicate flag, IPD multi-peak fraction ≤2% and harmonic amplitude <0.1. Of 33,322 sources that pass, zero have corrected point-estimate transverse speed above 1,000 km/s. With a three-sigma parallax threshold the count is 0/37,832. In both selections the maximum is 704.56 km/s, for the published hydrogen-atmosphere white dwarf SDSS J125834.93−005946.1. The formal uncertainty follows the distance as well as proper motion. A cut on this point estimate has different coverage from a lower confidence bound. [Igoshev et al.](https://arxiv.org/abs/2209.09915)

Only two of the twenty published comparison sources have exact DR3 matches anywhere in this particular catalogue. Both are classified DZ. J1109+0001 has corrected parallax S/N 1.32 and J1240+6710 has transverse speed 423 km/s. Both enter the separate [composition-selected branch](../metals/README.md). The precise-parallax result bounds the observed catalogue tail and cannot measure the intrinsic fraction of thermonuclear remnants.

The data source is the curated v1 main catalogue. The [Zenodo record](https://doi.org/10.5281/zenodo.23453342) and [Warwick distribution](https://cygnus.astro.warwick.ac.uk/phsdaj/WhiteDwarfs/) have the same compressed MD5. The per-batch SHA-256 hashes in [receipt.json](receipt.json) cover all Gaia replies. The numeric sort, 220-ID batch size and exact query are in [fetch.py](fetch.py), allowing reconstruction without committing the bulk downloads. From the study directory, run the following commands with one writable cache directory.

```sh
python desi/fetch.py /path/to/survivor-cache
python desi/screen.py /path/to/survivor-cache
```

[Catalogue paper](https://doi.org/10.1093/mnras/stag1701) · [DESI DR1 data license and acknowledgments](https://data.desi.lbl.gov/doc/acknowledgments/) · [Catalogue counts](summary.json) · [Per-source tail](tail.csv) · [Matched controls](control_matches.csv)
