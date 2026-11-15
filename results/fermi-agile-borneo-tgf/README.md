# One Borneo terrestrial gamma-ray flash viewed from directions 93° apart

On 2016-03-20 at 07:55:50 UTC, Fermi GBM recorded a submillisecond gamma-ray pulse while flying south of Java. AGILE MCAL recorded a pulse north of the equator 7.465 ms earlier on its spacecraft clock. The two photon-selected catalogs independently recover the event. Their radio association tables point to the same WWLLN lightning location in Indonesian Borneo, 2.6516° S, 111.6029° E. At an assumed source height of 15 km the spacecraft lay **93.06° apart as seen from the flash**. The angle remains above 91.9° across the source-location and height grid below.

The new result here is the cross-identification and its viewing geometry. Both observatories had already cataloged their individual detections. Stereo TGFs have precedents. Gjesteland et al. reported three RHESSI/Fermi events using observations ending in 2014. Marisaldi et al. reported ASIM/Fermi and ASIM/AGILE events in 2018–2019. Searches by both identifiers, event time/date, the instrument pair and “stereo”/“simultaneous” did not turn up a prior analysis of this particular pairing. Publication priority remains open. [Gjesteland et al. 2016](https://meetingorganizer.copernicus.org/EGU2016/EGU2016-11968-1.pdf) · [Østgaard et al. 2019](https://doi.org/10.1029/2019JD031214) · [Marisaldi et al. 2020](https://doi.org/10.5194/egusphere-egu2020-12804)

## What the instruments measured

| Observation | Measurement |
|---|---|
| Fermi GBM offline detection | oTGF160320330, discovery window starts 07:55:50.783900 UTC and spans 200 μs. Catalog counts are 5 and 4 in the two BGO detectors and 15 summed over NaI. GBM reports a corrected joint Poisson probability of 2.11 × 10⁻¹³. This tests the GBM count burst. |
| Fermi hourly event files | The all-channel extraction gives 5 and 6 events in BGO 0 and 1, against 0.297 and 0.298 expected from adjacent exposures. Nine events have channel lower boundaries above 300 keV, five in BGO 0 and four in BGO 1. Three channels have lower boundaries above 1 MeV. The all-channel stream includes two additional low-energy events. Energy intervals describe recorded detector channels. |
| AGILE MCAL photon selection | 20160320.075550.776435, orbit 46084. The independent Part II event selection reports 11.0 ± 3.3 fitted counts, a single pulse, and a middle-50% duration of 27 ± 5 μs. The individual event file contains 11 deposits within ±100 μs of the fit epoch, spanning 69 μs (−25 to +44 μs). Three exceed 10 MeV, with a maximum of 17.270 MeV. A twelfth deposit at +610 μs falls outside the prompt window. |
| WWLLN field | The GBM event file contains 556 radio strokes in the surrounding 20 minutes. The matched stroke has a travel-corrected absolute offset of 208 μs. The next closest in time is 1.044 s away. AGILE Part I lists the same coordinates rounded to 0.01°. Both catalog checks use one radio network. |

AGILE classifies this event as a single gamma-ray pulse. The GBM terrestrial-electron-beam list and special-comments file have no entry for it. The official MCAL time/energy plot labels the event table’s units and records the fitted epoch and orbit. [AGILE event-level counts](https://www.ssdc.asi.it/mcal3tgfcat/lc2/20160320.075550.776435.dat) · [AGILE light curve and energy scatter](https://www.ssdc.asi.it/mcal3tgfcat/lc2/20160320.075550.776435.pdf)

The Fermi radio offset meets GBM’s 3.5 ms association criterion and narrowly exceeds its stricter 200 μs “simultaneous” convention. A common lightning flash is well supported. The radio data alone cannot identify the microscopic emission episode within the flash. [GBM field definitions and selection](https://fermi.gsfc.nasa.gov/ssc/data/access/gbm/tgf/) · [AGILE catalog and samples](https://www.ssdc.asi.it/mcal3tgfcat/) · [AGILE photon selection](https://doi.org/10.1029/2019JD031986)

## Cross-catalog search and false matches

We compared all 4,135 GBM offline events against the independently photon-selected AGILE Part II samples. The latter contain 5,344 entries and 5,249 distinct photon times. In the overlapping timestamp range there are 980 Fermi events and 1,193 distinct AGILE photon times. Only the Borneo event has a match within 100 ms. The second-closest match is 543.498 ms apart, with the spacecraft nadirs 19,631 km apart. AGILE’s radio-selected entries were excluded from the discovery search and used only after matching.

Shifting the entire timing comparison by every nonzero integer second from −1,000 to +1,000 preserves short-term event clustering. With a fixed 1,000.1-second edge margin, 17 of 2,000 shifted catalogs have at least one match within 100 ms anywhere in the catalog. The unshifted, edge-matched search has one. This measures the frequency of timing-only false matches on these shifts. No spatial or radio evidence enters that count. The shifts share data and do not provide an exact independent-trial significance.

## Geometry and timing

WGS84 Earth-centered positions give spacecraft separation of 1,052 km, a 979 km surface distance between nadirs, and source-to-spacecraft distances of 836 km (Fermi) and 596 km (AGILE). The respective viewing angles from the source’s local zenith are 54.05° and 39.41°. For source heights from 10–20 km in 1 km steps, and horizontal displacements of 0, 5, 10 and 15 km around the radio fix in 5° bearing steps, the separation spans 91.93° to 94.15°. This is a sensitivity range, conditional on the stroke localization. AGILE’s 0.05° longitude difference between its photon and radio catalog entries has negligible impact at the precision of the headline.

The photon path to Fermi is 0.799 ms longer in this model. After subtracting that difference, the catalog epoch difference is 6.666 ms. AGILE was in its documented GPS drift period, when systematic errors persisted for days and absolute times could differ by tens of milliseconds. AGILE reports a 6.572 ms travel-corrected radio offset for this event. The catalog epochs refer to different pulse estimators. They cannot support submillisecond measurements of the intrinsic lag. [AGILE timing calibration and radio localization](https://doi.org/10.1029/2019JD031985)

The observed directions span more than a right angle. A single circular cone containing both direct rays would require a half-angle of at least 46.53° at the nominal location. Atmospheric scattering, tilted or multiple emission regions and multiple pulses within one lightning flash can change the relation between this viewing angle and the electron or primary-photon beam. The archive measures the viewing geometry. Inferring an intrinsic beam profile would require photon transport and instrument response modeling.

## Reproduce

From this directory, run the following with Python, NumPy, and Astropy. Allow about 65 MB for nine public input files. The script writes the field values, source URLs and SHA-256 hashes to [evidence.json](evidence.json) and the eleven all-channel BGO events to [gbm_photons.csv](gbm_photons.csv). The script reads the two hourly Fermi FITS streams and the AGILE event-level time/energy file directly, and preserves every row of the latter in the evidence record.

```sh
python analyze.py --cache /tmp/fermi-agile-borneo --output /tmp/fermi-agile-result
diff evidence.json /tmp/fermi-agile-result/evidence.json
diff gbm_photons.csv /tmp/fermi-agile-result/gbm_photons.csv
```
