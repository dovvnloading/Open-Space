# Research checkpoint: new SPHEREx moving-object spectral lane

**Bottom line:** a reproducible, target-centered SWAN C/2025 R2 spectral-excess lead repeats on December22, December25 and held-out January2. This is not yet a validated molecular identification, new discovery or unusual-composition result. Calibration-aware testing has exposed limitations that prevent promoting it as one.

## Completed
- Frozen untouched-in-repository instrument lane and early branch checkpoint, then draft PR.
- 29P known positive control:90 attempted D4–D6 records,82 retrieved,8HTTP503;15 strict/75 inclusive-valid2.5px apertures.
- Echeclus sensitivity pilot:8selectedD5,5retrieved3HTTP503; insufficient continuum signal, no gas-absence inference.
- Wierzchos known-chemistry comparison:15availableD5,9retrieved6HTTP503. CO2 previously published byJWST, so not novel.
- SWAN completeD5 extension:57attempted,56retrieved1HTTP503;27strict/46inclusive-valid2.5px apertures. First9exploratory;remaining48selected before inspection. FullJSON, finalCSV, exactURLs/hashes preserved.
- Official spectral calibration tar and hashedresponse/map products; device/sky coordinate convention checked; approximateembeddedWCS replaced for response modeling.
- Aperture/background sensitivity; midpoint agreement44microseconds; targetpositions and ePSF profiles; AllWISE contamination cross-check; empiricalshape forwardmodel and heldouttest.
- Four focused unit tests pass; no browser suite/fullrepoCI invoked.

## Evidence and falsifications
December25 excess survives three background annuli and is centered near the ephemeris. January2 independently sampled the rising feature again. A simple continuum cannot fit these formal measurements. But the six-point December Gaussian fit has six parameters and zeroDOF, and its frozen shape does not transfer toJanuary2 under formalnoise (chi-square44.8/2DOF; largeraperture91.6/2). An assumed10%floor changes the answer but is not a measuredsystematic. A single4.260micron delta-line does not reproduce red-wing behavior. Thus proper molecularband/extendedsource modeling is required, not a confidentCO2 claim.

Point-source morphology is inadequate in brightfeaturechannels; Gaussian broadening helps, but continuum also broadens. Do not call this resolvedgas. October has roughly3pixelmotiontrailing and strongcontamination, so exclude unmodeledseasoncomparisons. September band-adjacent measurements are weak; largerheliocentric/geocentricdistances and extractiongeometry preclude a simple physicaldeclineclaim.

All Decemberpilot featureapertures require explicitly permittingOUTLIERflags. Anothercoreflag rejects theDec22bluefeatureimage even afterrecovery. No suchflags were overridden. Nearbyhistoricalcatalogsources and brightstarwings require deblending; catalogbroadbandflux is not contemporarynarrowbandphotometry.

## Next necessary work
Physical fluorescentbandtemplate throughaperture-weightedresponse, commoncoma/trail+star model, inject/recover uncertaintycalibration, and independentcross-instrument/epochvalidation. Literature search didnotfind indexedSWANCO2/SPHERExpaper, but absenceofsearchhit isnotpriorityproof. Do not report unusualCO2/H2O, productionrate, or firstdetection.

Start with SWAN_CANDIDATE.md, SWAN_SYSTEMATICS.md, RESPONSE_VALIDATION.md, EMPIRICAL_FIT_REPORT.md, and independent spatialreviews. Files ending in pilot are deliberately retained frozencheckpoints, not thefinal57sample.
