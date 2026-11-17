# A fixed field source reproduces the strongest marginal W2 associations

**All five wide-screen cases with at least six clean predicted-track matches spanning at least three held-out visits, but falling short of the full motion criterion, encounter a previously present field source.** This additional test allows the fixed alternative to be centered on the field source near the late track rather than at the proposed object's starting position. All the tabulated moving-window hits are within 3.5 arcsec of that fixed field position.

| CatWISE ID | Held-out clean hits within 2″ of moving model | Hits within 2″ of free fixed model | Moving-model hits within 3.5″ of fixed model |
|---|---:|---:|---:|
| J102015.59−593344.0 | 15 | 109 | 15 |
| J102036.54−594326.6 | 16 | 117 | 16 |
| J103528.91−591635.2 | 11 | 55 | 11 |
| J152202.65−510741.7 | 12 | 89 | 12 |
| J152353.76−512725.9 | 12 | 56 | 12 |

The fixed peak was independently present before 2019. Each pair of counts uses the same quality-filtered exposures beginning on 2019-01-01. The broader 3.5-arcsec diagnostic estimates overlap of the two hypotheses' association regions. A close passage by a fainter moving object remains possible; counting these detections cannot establish one. The [per-source records](results.json) contain early counterpart counts, fitted coordinates and URLs, and [fixed-centroid measurements](fixed-centroids.csv) report available half-year bins.

After running the [wide exposure analysis](../widebatch/README.md) to regenerate the loose receipts, invoke `python ../fast-alternatives/static.py --l1b ../widebatch --out . --names J102015.59-593344.0 J102036.54-594326.6 J103528.91-591635.2 J152202.65-510741.7 J152353.76-512725.9` from this directory. The [method description](../fast-alternatives/README.md) gives the exact peak choice and aperture sizes.
