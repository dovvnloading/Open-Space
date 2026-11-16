# Nearby white-dwarf transit screen

A search for deep recurring dips found apparent near-total flux losses that fail field and detector-position checks. No astrophysical occultation is established.

## Sample and query receipts

The southern Gaia EDR3 catalogue parent comprises 7,372 objects with declination below zero, parallax ≥10 mas, white-dwarf probability ≥0.90, and 16 ≤ G ≤ 20.5. The original Harvard VizieR response with query URL, column metadata, catalogue DOI and retrieval timestamp is frozen in `parent.tsv.gz`. The astrometric parent comes from [Gentile Fusillo et al. (2021)](https://doi.org/10.1093/mnras/stab2672).

Propagate each Gaia position from epoch 2016.0 by 10.6 years using both proper-motion components and search a 1.512″ radius in ALeRCE's Rubin object table (`sid=1`). Across 154 batches, 176 DIA objects match 169 white-dwarf positions. These associations include split detections and dipoles of high proper motion sources; a match alone has no astrophysical interpretation. The 176 objects returned 2,151 detections and 3,435 forced measurements. An absence of alert history does not establish sky coverage or constant flux. The known debris transiter ZTF J0328−1219 has no ≥12-forced-epoch ALeRCE object within the separate control cone and supplies no recovery efficiency. [Known debris transit characterization](https://doi.org/10.1093/mnras/stab3475)

## Why the dips fail

WDJ040653.64−485428.82 (Rubin 313853497407176716) and WDJ041617.75−463632.06 (313897383915683844) have six and two forced epochs with direct-image flux below half the object's same-band median and difference flux below −0.35 times that median. None has a same-object DIA detection on the corresponding visit/detector. Multiple other white-dwarf positions in the same field show similar dropouts.

For each of these eight epochs, compare against the nearest 59 other Rubin DIA objects with at least 12 forced measurements within 0.10°, retaining sources on the same detector in the same visit with ≥5 band measurements and median science flux >15 times that epoch's formal flux error. On visit 2026021900078, the first target fades 99.5%; 10/19 neighbours also lose more than half their flux, with an 84.9% median loss. Other visits show more local effects. A local affine fit to sky coordinates and measured CCD coordinates of up to 16 nearest DIA sources from each visit maps two events for the second target to CCD rows 1 and 0 (fit residuals 0.60 and 0.44 pixels). Two events for the first target fall at CCD columns 5 and 6. One of its six epochs lacks three astrometric references. These repeated low measurements cannot establish eclipses without calibrated pixel/mask confirmation. The classification of the remaining individual drops is unresolved. No occurrence rate follows from this alert-conditioned sample.

Run `python check.py` to regenerate the sample counts and local-control summaries from frozen measurements. The field-control extract contains each object/measurement ID, band baseline, measured flux, uncertainty and fitted detector location where available. In `tap_receipts.json.gz` each query has the literal ADQL, raw returned CSV digest and record count. Replay against [ALeRCE TAP](https://tap.alerce.online/tap/sync) with REQUEST=doQuery, LANG=ADQL, FORMAT=csv and MAXREC=200000. The live database can evolve or return a different row order.
