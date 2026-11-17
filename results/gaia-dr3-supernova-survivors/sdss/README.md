# SDSS-V DR19 and DR20 spectroscopic classifications

We read every row in each release's Astra AllStar SnowWhite FITS product and retained every nonblank classification with a positive Gaia DR3 identifier, including uncertain or mixed classifications and CVs. DR20 includes many MWM rows on which SnowWhite provided no classification. Those rows are outside this white-dwarf-programme census. A classification is a model output. All 51,379 distinct identifiers in the union matched the live Gaia DR3 archive. Of these, 11,196 had already been retrieved for the DESI search. We retained one record per star per release and used the highest-S/N spectrum for its display label. There is no S/N cut on the spectra in this census. [SnowWhite classifier and programme targeting](https://www.sdss.org/dr20/mwm/astra/pipelines-in-astra/snowwhite/)

| Release | Rows read | Classified rows | Distinct classified DR3 IDs | Five-sigma and astrometric-quality eligible | Measured transverse maximum |
|---|---:|---:|---:|---:|---:|
| DR19, Astra 0.6.0 | 10,364 | 10,364 | 10,364 | 7,843 | 486.90 km/s |
| DR20, Astra 0.8.1 | 1,803,883 | 50,960 | 50,862 | 36,400 | 601.97 km/s |

We used the same source-specific zero point and five astrometric quality conditions as the [DESI census](../desi/README.md). Both releases have zero eligible sources with a point estimate above 1,000 km/s, including when the corrected parallax S/N threshold is lowered to three. The DR20 three-sigma maximum is 681.18 km/s. No exact-ID match to the twenty class controls has a nonblank SnowWhite classification in these products. Other SDSS programs and earlier BOSS spectra give coverage for some controls. The two releases share 9,847 classified Gaia IDs. Survey counts cannot be added.

From the study directory, first run the DESI exact-ID download, then run `python sdss/screen.py /path/to/survivor-cache`. The script checks the FITS file hashes and uses the compressed source files to recreate its local memory maps. [Release files, denominators and hashes](summary.json) · [Additional Gaia batch receipts](gaia_receipt.json) · [Fast-tail coadds](tail.csv) · [SDSS acknowledgement](https://www.sdss.org/collaboration/citing-sdss/)
