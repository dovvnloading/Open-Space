# Why the J1539+3855 pair does not establish two polluted white dwarfs

Gaia EDR3 1377401818623167104 and 1377401822919131008 are separated by 3.35″. They appear as a WD+WD pair in the wide-binary catalogues, with a projected separation of roughly 1,100 au and chance-alignment statistic 0.000956. The Planetary Enriched White Dwarf Database carries calcium entries for the brighter star, and the Montreal White Dwarf Database (MWDD) lists the fainter star as DZ. The 2026 DESI DR1 catalogue classifies the brighter star as DC at blue S/N 8.1. Its DESI spectrum yields no persuasive Ca K absorption. Crowding, heterogeneous reductions and the original low-S/N spectra make the catalogue labels insufficient evidence.

The two separately named SDSS spectral exports share fine-scale fluctuations. We converted the brighter export from the per-frequency spectral shape to per-wavelength, interpolated the fainter export to that wavelength grid and subtracted a 21-pixel median in each spectrum. Across 871 common pixels from 4500 to 5500 Å the residual correlation is **0.892**. Shifting one spectrum by 5, 10, 20 or 40 pixels gives correlations of −0.094, −0.092, −0.005 and +0.003. Strong aligned fine structure is consistent with shared input data or extraction contamination. Without distinct exposures and extraction provenance for both objects, these records cannot support two independent detections. A resolved spectrum of each component, including Ca K, would test the hypothesis.

An additional positional match pairs the polluted star WD 1425+540 with G200−40. Published resolved spectroscopy identifies the companion as a K dwarf. A wider angular match lies in the NGC 6633 field and does not establish a bound binary. No double-pollution discovery is claimed. These outcomes close the specific search preserved in [the initial checkpoint](SEARCH.md).

Replay the decisive common-spectrum check with Python, NumPy and SciPy. The command downloads two SHA-256-verified spectra and prints [the recorded result](spectrum-check.json).

```sh
python replay.py
```

The investigation used the [Montreal White Dwarf Database](https://www.montrealwhitedwarfdatabase.org/), [PEWDD](https://github.com/jamietwilliams/PEWDD/tree/023dba3176ac367c6cf2b5968847a80d79e3c594), [Heyl's wide systems](https://zenodo.org/records/12803292), [Heintz et al. pair sample](https://arxiv.org/abs/2206.00025), [DESI WD release](https://zenodo.org/records/23453342) and the [G200−40 companion study](https://academic.oup.com/mnras/article/503/2/1877/6132259). The two spectral input URLs and checksums are in [sources.json](sources.json).
