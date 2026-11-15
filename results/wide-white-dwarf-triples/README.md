# Correcting Gaia identifiers for ten wide white-dwarf triple systems

Eleven component identifiers in Heintz et al. (2022), machine-readable double-white-dwarf pair table, point to a different member of the same triple. The errors affect 11 of 26 pair rows representing 10 of the 24 distinct triple systems. Ten rows assign the same ID to both stars. The eleventh row substitutes the wrong edge in the three-white-dwarf system. A coordinate-and-photometry check against Gaia EDR3 resolves all eleven fields. The separate printed triple table in the source paper contains the correct identifiers. The demonstrated impact is on joins by source ID. This audit does not re-evaluate the ages that the authors inferred from the component photometry.

## A false double-helium classification

Joining row 1582 by its two published IDs retrieves Gaia 5874024769842933760 twice. Vincent et al. (2024) classify this star as DB with a score of 0.99. The second position and G magnitude resolve to 5874024774146311808, classified DA with a score of 1.00. The original paper's spectral types, DBA and DA, corroborate this distinction. The two catalog classifications are preserved in `xp-classifications.csv`. The two stars are 2.57 arcsec apart. An ID-only join would both fabricate a DB+DB pair and remove the real hydrogen-atmosphere counterpart. The correction recovers known objects.

## Scope and falsification

The input is the original IOP fixed-width table (SHA-256 `a741e22aac89e2bdf28c5376a44e380cc48b7de2bac686b5e983a93f74dc9128`). Preserve the ID fields as decimal strings. The first data row is row 1. All coordinates and proper-motion reference epochs are Gaia J2016.0. For each component in every triple, look for all EDR3 sources within 0.54 arcsec of the published position, without filtering on source ID or G magnitude. Exactly one result meets the position criterion in each of 52 component records, representing 49 distinct white dwarfs. The largest positional difference is 0.000217 arcsec and the largest independent G-magnitude difference is 0.000048 mag. The angular precision reflects matching the same astrometric measurements across tables.

For a negative control, look up the published component IDs in all 1,565 ordinary binary rows. All 3,130 positions agree within 0.000246 arcsec and all G magnitudes agree within 0.000050 mag. The defect is confined to the triple subset in this release under this coordinate test. We also checked the ten distinct replacement IDs against the printed triple table. In all ten cases it lists the replacement at the matching component position. The publication already includes those correct identifications. We found no separately published correction in a targeted search for the authors, title, erratum and EDR3 component identifiers.

| Machine table row | Field | Published Gaia EDR3 ID | Corrected Gaia EDR3 ID |
| ---: | :--- | ---: | ---: |
| 1574 | Secondary | 338979301700437120 | 338979297405358848 |
| 1575 | Secondary | 3681220906102862720 | 3681220906103896064 |
| 1576 | Secondary | 4438147422454430464 | 4438147418158270848 |
| 1578 | Secondary | 5041268200697001600 | 5041268204992049664 |
| 1579 | Secondary | 5250998119556563200 | 5251022416193428096 |
| 1581 | Primary | 5560518000048274304 | 5560518000052459520 |
| 1582 | Secondary | 5874024769842933760 | 5874024774146311808 |
| 1583 | Primary | 6087695613244135936 | 6087695617537911040 |
| 1586 | Primary | 763543308609191936 | 763543312899793664 |
| 1589 | Primary | 4190500054845023488 | 4190499986125543296 |
| 1590 | Primary | 4190500054845023488 | 4190499986125543296 |

## Reproduction

Run `python audit.py` to summarize the frozen triple projection in `triple-components.csv`. Run `python audit.py --online` to download the checksummed publisher table, reproduce the complete triple sky search, check the frozen rows, and verify all ordinary binary components against Gaia EDR3. Both commands use the Python standard library. The Gaia archive requests use the Heidelberg TAP service. The CSV stores the original projected fields, the matched Gaia source, its coordinates and G magnitude, and the number of Gaia sources found inside the matching radius.

Heintz et al., *Testing White Dwarf Age Estimates Using Wide Double White Dwarf Binaries from Gaia EDR3* ([paper and printed triple table](https://arxiv.org/pdf/2206.00025), [original machine-readable table](https://content.cld.iop.org/journals/0004-637X/934/2/148/revision1/apjac78d9t2_mrt.txt), [VizieR mirror](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/934/148)). Gaia Collaboration, *Gaia Early Data Release 3* ([release](https://doi.org/10.1051/0004-6361/202039657), [Heidelberg Gaia TAP](https://gaia.ari.uni-heidelberg.de/tap)). Vincent et al., *Classification and parameterization of a large Gaia sample of white dwarfs using XP spectra* ([spectral catalog](https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/682/A5)).
