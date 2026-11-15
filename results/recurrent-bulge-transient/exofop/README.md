# Prior identification of the 2016 event in the NASA ExoFOP K2 Campaign 9 archive

NASA's March 2022 database export identifies **K2C9-R-0367** as **OGLE-2016-BLG-0935**. Event-parameter row 443 gives J2000 coordinates 271.2724583330°, −27.8226666667° and the ARTEMiS status `pspl` (point-source, point-lens). The campaign archive lists two ARTEMiS model products, two pyLIMA products and the OGLE finding chart. Its changelog contains zero rows matching either alias. The NASA archive establishes prior recognition and lens modeling of the 2016 event. Neither the candidate's presence in Campaign 9 nor this model establishes an actual lens. [NASA release](https://exoplanet.ipac.caltech.edu/docs/exofop_k2c9.html) · [Table definitions](https://exoplanet.ipac.caltech.edu/docs/README_k2c9_tables.txt)

The KMTNet Campaign 9 table independently identifies KMT-2016-BLG-0228 with OB160935. KMTNet calls it a “possible” microlensing event. Henderson's 2017 Keck Campaign 9 presentation shows the May 26 OB160935 alert and OGLE light curve. McDonald et al. (2021) report a rejected variable at JD 2457535.478, but its position 271.3306°, −27.2857° lies 0.54° away. Their published paper and tables contain no entry matching the target coordinates or the 0935 alias. [KMTNet event table](https://kmtnet.kasi.re.kr/~ulens/event/2016k2/) · [Keck presentation](https://nexsci.caltech.edu/conferences/2017/microlensing/k2c9_keck_too_cbhenderson_share.pdf) · [McDonald et al.](https://arxiv.org/abs/2107.02746)

## Catalog checks

A cone of radius 0.00556° (20.016 arcseconds) centered on 271.272458°, −27.82268° returned no rows in either the OGLE Collections of Variable Stars time-series table or SIMBAD. The nearest of the 1,091 objects in the 2015 OGLE dwarf-nova catalog is OGLE-BLG-DN-0658 at 271.288875°, −27.79947222°, 98.55 arcseconds away. A source can be absent from these catalogs through selection or crowding. [OGLE service](https://skvo.science.upjs.sk/ogle/q/ts-web/form) · [SIMBAD TAP](https://simbad.cfa.harvard.edu/simbad/sim-tap) · [Mróz et al. catalog](https://vizier.cds.unistra.fr/viz-bin/VizieR?-source=J/AcA/65/313/table1)

Send these respective ADQL queries to `https://skvo.science.upjs.sk/tap/sync` and `https://simbad.cfa.harvard.edu/simbad/sim-tap/sync` with `REQUEST=doQuery`, `LANG=ADQL` and `FORMAT=csv`.

```sql
SELECT TOP 30 ssa_targname, ssa_location, mean_mag FROM ogle.ts_ssa
WHERE 1=CONTAINS(ssa_location, CIRCLE('ICRS',271.272458,-27.82268,0.00556))

SELECT TOP 30 main_id,otype,ra,dec FROM basic
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',271.272458,-27.82268,0.00556))
```

## Replay

`exofop_rows.json` preserves the full selected rows, including the empty changelog result. The standard-library Python script checks the SHA-256 of the compressed release and of each extracted table. The archive is 597,396,749 bytes. Its SHA-256 is `69ebcfdd40d85a1c911e314b359fe2f21444da0deac64bbe98d7d48c21a0d2eb`. [Source archive](https://exoplanet.ipac.caltech.edu/data/Contributed/K2C9/K2C9.tar.gz)

```sh
python extract_exofop.py > replay.json
diff -u exofop_rows.json replay.json
```

Use `--archive /path/to/K2C9.tar.gz` to replay from an existing local copy.
