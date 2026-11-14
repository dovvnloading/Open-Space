# EROS-2 eclipse measurements before the OGLE-IV survey

Public EROS-2 observations from 1996–2003 provide an earlier test of the three OGLE-IV eclipse-depth leads. We froze the OGLE catalog coordinates, orbital periods and primary epochs. We retrieved all 39 EROS source records and 28,787 paired exposure rows within the three 8-arcsecond cones. We verified the exposure count for every record. Four records lie within 2 arcseconds of the three OGLE positions. OGLE-SMC-ECL-07730 falls in two distinct EROS pointings. [EROS data release](https://doi.org/10.1051/0004-6361/202557321) · [OGLE-IV binary catalog](https://ogle.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/) · [Results and receipts](../results/smc-eclipse-eros.json)

We pooled the primary and secondary eclipse cores, within 0.05 cycle of phases 0 and 0.5, and measured their mean magnitude excess over points at least 0.13 cycle from either minimum. The robust fit accounts for year-bin offsets. Standard errors cluster residuals in seven-day blocks. Both passbands share observing times, and both records of 07730 sample the same sky. The 2-arcsecond matching rule and the EROS resolution allow flux from a neighbor.

| OGLE-SMC-ECL | EROS record | Offset | Red contrast ± SE (mag) | Blue contrast ± SE (mag) |
| --- | --- | ---: | ---: | ---: |
| 07730 | sm0030k7615 | 0.14″ | +0.095 ± 0.033 | +0.041 ± 0.011 |
| 07730 | sm0046l6504 | 0.51″ | +0.019 ± 0.030 | +0.026 ± 0.011 |
| 06703 | sm0010n23802 | 1.84″ | −0.001 ± 0.011 | +0.019 ± 0.010 |
| 09564 | sm0063k22812 | 0.32″ | +0.024 ± 0.008 | +0.021 ± 0.004 |

At the OGLE orbital ephemerides, the EROS observations show modest dimming toward 07730 and 09564. For 07730 a scan of 100 phase offsets using only each record's early red and blue observations selected +0.030 and +0.015 cycles. With those offsets locked, the held-out late blue contrasts were 0.052 ± 0.015 and 0.033 ± 0.018 mag. The corresponding late red measurements were 0.056 ± 0.060 and 0.042 ± 0.037 mag. For 09564 the early data selected +0.005 cycles and the later blue contrast was 0.016 ± 0.006 mag. Toward 06703 the older time series and 1.84″ association are inconclusive. The JSON retains all cone neighbors, phase/width/error-cut sensitivities and the full holdout results. Photometric selection could suppress deep or faint measurements.

The same phase-window estimator measured 07730 OGLE-I contrasts of 0.181, 0.062, 0.090 and −0.004 mag across the four OGLE-IV epochs. For 09564 they were 0.044, 0.095, 0.110 and 0.148 mag. EROS-2 has different filters, a different angular response and gaps before OGLE-IV. The physical depth of either eclipse in the intervening years remains unconstrained. A shared cycle of disappearance and reappearance would require further evidence.

### Time standard and reproduction

The EROS paper's Appendix A.4 defines the date as heliocentrically corrected MJD minus 49999.5, giving **HJD = EROS date + 2450000.0**. The VizieR image-table column description instead gives 2449999.5. To resolve this, we matched a light-curve timestamp to its actual image. For sm00300trr9e17119.fits the image archive gives MJD 51316.40930 and EROS 1316.90988. Adding 2450000.0 gives a heliocentric correction of +0.00058 d (50 seconds). Using the column label would imply −0.49942 d. The script repeats the check before writing the result. A half-day offset changes the conclusion for these short-period systems. [Image archive query](https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=II%2F390%2Fxtinfo&HJD=1316.9098..1316.9100&-out=FileName%2CMJD%2CHJD%2CDate) · [EROS release paper, Appendix A.4](https://arxiv.org/abs/2512.12745)

Run from the repository root with NumPy installed, passing the OGLE-IV ecl.dat next to the extracted phot/I directory. The script downloads EROS only on cache miss. The JSON records hashes of the OGLE inputs, the EROS raw responses and date-independent sorted data rows.

```sh
python examples/smc_depths/eros.py runs/smc-eros/cache runs/smc/raw/ecl.dat runs/smc-eros/result.json
```
