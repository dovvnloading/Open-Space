# Earlier investigations and source identities

The review on 2026-11-16 covered the main tree at `a62807165dd86b1c405751042dbd0a98be00a777`, the main README and the titles and bodies of the 106 pull requests then present. Independent open and closed queries recovered four open PRs (#102, #104, #105, #106) and 102 closed PRs. Ascending and descending creation-order queries covered both ends of the 100-row search limit. The untruncated main tree contained 4,029 entries and 182 top-level research directories under `results/`. The frozen PR inventory records the title, state and URL of every PR. PR #107 opened during the search. Its radio turnover study of NVSS J100005+812702 was checked separately and added to the inventory before closing this audit, bringing the total to 107. [PR #107](https://github.com/dovvnloading/Open-Space/pull/107)

No PR described a 5XMM thermal isolated-neutron-star search. Existing work does include source classes this search can select. PR #12 examined the M31 supersoft transient XMMU J004347.8+413310; #15 examined nova sightlines in M31. PRs #13, #16, #21, #23, #31, #47, #50, #61, #67, #86, #88 and #89 include pulsars, magnetars or compact-source timing. PRs #17, #18, #20, #42, #44, #46 and #59 include X-ray variables, GRBs, supernovae or ULXs. The main tree also includes two eROSITA-selected cataclysmic-variable candidates, the M83 X039 flare, 4U 1538−522, GX 1+4 and SMC X-1. The new selection uses the full DR15 source table rather than seeding from those campaigns. Positional checks at thirteen named repository targets supplement this investigation-level review. Many timing and population studies have large comparison sets; those sets were not exhaustively converted to an all-object positional veto. [Repository snapshot](https://github.com/dovvnloading/Open-Space/tree/a62807165dd86b1c405751042dbd0a98be00a777) · [PR inventory](pr-inventory.csv)

## Reproducible position checks

`prior-positions.csv` includes the following position records. Aliases and revised coordinates are retained across papers, so the number of rows exceeds the number of astrophysical sources.

| Earlier work | Position records | Treatment |
| --- | ---: | --- |
| Pires et al. 2009, Table 2 | 8 | Positions and 90% errors of all eight optical follow-up targets. [Paper](https://doi.org/10.1051/0004-6361/200912180) |
| Rigoselli et al. 2022, Table 3 | 47 | All soft-source designations, including known INSs, contaminants and ambiguous objects. Coordinates decoded from designations are approximate. [Paper](https://doi.org/10.1093/mnras/stab2974) |
| Kurpas et al. 2024, Table 2 | 33 | All 33 eRASSU sources with the paper's positions and uncertainty convention. The two early follow-up targets are already included. [Paper](https://doi.org/10.1051/0004-6361/202450248) |
| Demasi et al. 2024, Tables 2 and 4 | 55 | Full ROSAT blank-field shortlist joined to published positions and errors. Includes prior candidates with subsequently identified contaminants. [Paper](https://doi.org/10.3847/1538-4357/ad03e7) |
| Pires et al. 2026, Tables 1 and A.1 | 10 | All ten investigated sources, including the ambiguous and rejected cases, with corrected coordinates and 90% errors. [Paper](https://doi.org/10.1051/0004-6361/202659695) |
| Open-Space | 13 | Selected overlapping named repository targets with per-row references. |

Each DR15 queue was compared with all 166 rows by great-circle separation. The review envelope is 90 arcsec for ROSAT 1/2RXS and 30 arcsec for everything else. These intentionally broad radii flag aliases for inspection. A positional hit alone cannot establish physical identity. Each output includes all in-radius hits or the nearest reference and distance when none fall inside. An additional independent check used all 147 ROSAT bright-source study entries in the HEASARC mirror of Turner et al. 2010, including later rejected sources. Use `turner-prior.csv` in place of the core table for that replay. The published samples and the input catalogue share selection biases. Clearing the list does not establish universal priority. [ROSAT BSC study table](https://heasarc.gsfc.nasa.gov/W3Browse/all/rasscndins.html)

| DR15 queue | Rows | Distinct DR15 positions within core envelope | Within Turner envelope |
| --- | ---: | ---: | ---: |
| Measured HR2 | 185 | 6 | 0 |
| Boundary | 8 | 0 | 0 |
| Rescue | 12 | 0 | 0 |
| Bright sensitivity check | 116 | 1 | 0 |

The six measured-queue hits are J072608.1−261238, J081952.3−131932, J134725.6−363413, J022141.4−735634, J225139.6−162750 and J063045.4−603113. The bright-screen hit is J114051.9−641848. J022141.4−735634 matches two papers. J225139.6−162750 illustrates why a stable soft spectrum warrants a positional and extragalactic check. Pires et al. found a faint galaxy at the refined XMM position. [Pires et al. 2026](https://doi.org/10.1051/0004-6361/202659695)

To replay one queue with the standard library, pass its native CSV, the position list and an output filename.

```sh
python match_prior.py measured.csv prior-positions.csv measured-prior.csv
```

## Four fields checked individually

The full catalogues, query URLs, search areas and selected returned columns are in `target-history.json`. A successful empty catalogue cone is a statement about listed detections. Flux limits require exposure and sensitivity maps. The four-class Tranin et al. classifier provides AGN, star, X-ray binary and CV scores and has no isolated-neutron-star outcome. [Tranin et al. 2022](https://doi.org/10.1051/0004-6361/202141259)

- **J085618.5+380711.** 4XMM-DR12 and DR13 contain 4XMM J085618.5+380711, SRCID 203025818010014, at 0.499 arcsec. The DR12 band-2 flux of 2.295 × 10⁻¹⁵ erg cm⁻² s⁻¹ lies below the 10⁻¹⁴ erg cm⁻² s⁻¹ nominal 0.5–1 keV threshold in the 2026 Pires search. Tranin et al. list no prior class or optical/IR counterpart and give four-class model scores of 0.724 for a star and 0.264 for an AGN. The source had no SIMBAD object inside 90 arcsec, no CSC2, ATNF, Swift 2SXPS or LSXPS object inside 60 arcsec, and no Chandra aimpoint inside 20 arcmin. The nearest GALEX AIS entry is at 43.71 arcsec. The nearest prior in the core list is over 5.75 degrees away. Searches on the observation ID lead to the targeted galaxy cluster RXC J0856.1+3756 and studies of other field sources. The DR15 object lies 10.32 arcmin from the cluster position in the SSC target table. Searches on the XMM designations and abbreviated coordinates returned no earlier source-specific thermal INS identification. The prior catalogue detections use the same 2005 photons. [SSC pointing and target](https://xmmssc.irap.omp.eu/2XMMi/xcat_targets.html) · [DR13](https://cdsarc.cds.unistra.fr/viz-bin/cat/IX/69)
- **J180019.4+080504.** DR12 and DR13 list 4XMM J180019.5+080503 (1.255 arcsec). Its Tranin star score was 0.585 without an optical/IR association or input identification. No SIMBAD objects lie inside 90 arcsec. The V2301 Oph eclipsing polar targeted by this programme is 6.50 arcmin away. It remains important to check the actual field and relative fluxes of all contributing epochs. [ESA proposal](https://doi.org/10.5270/esa-mzgmbki)
- **J092212.2+032038.** DR12 lists 4XMM J092212.0+032040 (2.713 arcsec). Tranin et al. reported a Pan-STARRS association and an AGN model score of 0.959. An independently queried Legacy Survey source with resolved REX morphology, r_AB = 21.20, lies 0.77 arcsec from the DR15 position. This supports an extragalactic competing interpretation. No source-specific spectrum or time history establishing a tidal disruption event was identified. [DR12](https://cdsarc.cds.unistra.fr/viz-bin/cat/IX/68)
- **J052015.5−654420.** This detection in 0841320101 is 6.64 arcsec from its DR12 counterpart 4XMM J052015.1−654426, 6.95 arcsec from the published optical position and 9.38 arcsec from the SIMBAD eRASSU position. SIMBAD links both earlier X-ray designations with [HP99] 159. Greiner et al. analysed the very same exposure. They measured a ~45 eV soft spectrum, identified a ~16 mag optical counterpart, and located the helium-line source in the LMC. The inferred thermal radius is around 3,700 km at 50 kpc. Subsequent spectroscopy and photometry support a 2.327-day binary containing an evolved helium donor. This field is a direct example of a known supersoft white-dwarf system appearing without a close DR15 catalogue counterpart. [Greiner et al. 2023](https://doi.org/10.1038/s41586-023-05714-4) · [Follow-up](https://doi.org/10.1093/mnras/stag189)
