# Separate aperture photometry of the discovery exposures

McNeill et al. (2023) extracted asteroid photometry independently from TESS Cycle 1 full-frame images using forced circular apertures. We selected the two asteroid numbers from their public sigma-clipped epoch table. We kept their flux values, fixed the Pál-derived spin frequency, dip spacing, epoch and event width, and fitted only a quadratic baseline and the spin harmonics. [McNeill et al. (2023)](https://doi.org/10.3847/1538-3881/acf194) · [Original data and metadata](https://rcdata.nau.edu/snaps/published_datasets/TESSC1/) · [Selected original rows](photometry.tsv)

| Using the archive dates and six spin harmonics | (29489) 1997 WQ | (3177) Chillicothe |
|---|---:|---:|
| Usable measurements | 821 | 798 |
| Individual events dimmer in the predicted window | 32/33 | 22/23 |
| Median depth [cycle-bootstrap 95% interval] | 0.115 [0.095, 0.129] mag | 0.043 [0.037, 0.048] mag |
| Median depth after 28 spin harmonics | 0.114 mag | 0.017 mag |
| Adjacent pairs with the even event deeper | 6/12 | 5/8 |

The two catalogs timestamp the same half-hour image sequences differently. Our primary check uses McNeill's reported modified Julian dates plus 2,400,000.5 without a fitted lag. A second check joins exposures to the nearest Pál timestamp within 15 minutes, using only times. This shifts the median timestamp by +12.39 minutes in 1997 WQ and +12.93 minutes in Chillicothe. The matched-cadence depths are 0.120 and 0.044 mag, respectively. The dimmings survive both timestamp conventions and extraction pipelines. The sparse parity comparisons in the second reduction do not reproduce a consistent alternation. The telescopes, exposures, field backgrounds and observing geometry are common to both reductions. [All timestamp and spin-order tests](test.csv)

Run `OPENBLAS_NUM_THREADS=1 python validate.py` from this directory to reproduce the CSV offline. The tab-separated file is an exact selection from the release's five-column `3sigma_Cycle1.txt` table, using the first field as the asteroid number. Verify the full download and extract the rows with the following commands.

```sh
curl -fL https://rcdata.nau.edu/snaps/published_datasets/TESSC1/3sigma_Cycle1.txt -o 3sigma_Cycle1.txt
echo 'e8e619c45d82b061c327224bed0d3a26b1df9ba9e2f92730f3837f28143612ab  3sigma_Cycle1.txt' | sha256sum -c -
awk '$1 == 29489 || $1 == 3177' 3sigma_Cycle1.txt > retrieved.tsv
cmp photometry.tsv retrieved.tsv
```
