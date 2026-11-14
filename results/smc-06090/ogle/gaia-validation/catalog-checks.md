# Prior reports and catalog neighbors for OGLE-SMC-ECL-06090

OGLE-IV gives J2000 coordinates 01h24m15.06s, −73°19′33.4″ and the field identifier SMC738.16.542. The OGLE-III catalog lists OGLE-SMC-ECL-6090 and field star SMC121.3.2046 at the same coordinates. No OGLE-II counterpart is listed. Searches performed on 2026-11-13 US Eastern included both binary identifiers, both field aliases, the sky position and the Gaia DR3 identifier with eclipse timing and apsidal motion terms.

| Record checked | Result |
| --- | --- |
| [Pawlak et al. 2013, OGLE-III release](https://ogle.astrouw.edu.pl/cont/4_main/var/ogleiii/ebssmc/) | All 32 objects labeled with eccentric orbits and visible apsidal motion in [the release remarks](https://www.astrouw.edu.pl/ogle/ogle3/OIII-CVS/smc/ecl/remarks.txt) were checked. Source 6090 has no remarks entry. |
| [OGLE-IV SMC binary release](https://www.astrouw.edu.pl/ogle/ogle4/OCVS/smc/ecl/remarks.txt) | No 06090 remarks entry. |
| [Hong et al. 2015](https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/AJ/150/1/stars&-out=OGLE-III&-out.max=100) | None of the 27 source IDs is 6090. |
| [Hong et al. 2016](https://arxiv.org/abs/1604.02200) | No 6090 in the 90-system sample. Its I magnitude is fainter than the search's 18.0 mag cutoff. |
| [Zasche et al. 2014](https://arxiv.org/abs/1412.0905) and [Zasche et al. 2019](https://arxiv.org/abs/1902.04961) | 6090 occurs in neither apsidal study's target list, containing 18 and 21 systems respectively. |
| [Juryšek et al. 2018](https://arxiv.org/abs/1709.08087) | 6090 is absent from the 21 reported SMC cases of changing eclipse amplitude in Table 3. |
| SIMBAD and VSX | A 5″ cone in SIMBAD returned one object, OGLE SMC121.3 2046, linked to Gaia DR3 4686449235706841984 and five references. They are the OGLE and Gaia catalogs and [Bódi & Hajdu's morphology study](https://arxiv.org/abs/2106.01039). VSX returned one eclipsing object with the same Gaia source ID. Its listed period is 1.16968 d. Raw query replies are provided. |

The catalog and paper checks found no earlier report of this timing trend. They leave open additional catalogs, less common aliases, unpublished analyses and papers that do not identify individual targets in searchable text.

The [ESA Gaia archive](https://gea.esac.esa.int/archive/) query below returned the three rows included in gaia-neighbors.csv. Distances are relative to the OGLE catalog coordinate. The target lies 0.083″ away. The two remaining Gaia sources lie 3.65″ and 7.51″ away and are 2.78 and 2.75 G mag fainter. These are detected neighbors only. Unresolved companions and differences in OGLE I-band blending remain possible.

```sql
SELECT source_id, ra, dec, phot_g_mean_mag, phot_variable_flag
FROM gaiadr3.gaia_source
WHERE 1 = CONTAINS(
  POINT('ICRS', ra, dec),
  CIRCLE('ICRS', 21.06275, -73.3259444444, 0.002222222)
)
```
