# JWST companion search around nearby stars

**No new bound companion is established.** The observing unit is a target, epoch, detector, filter and mask. Repeat observations and previously unsearched fields receive priority within verified public data.

## Completed

- Reviewed the repository's [inclusion criteria](../../docs/REMARKABLE_FINDINGS.md) and [repeat-imaging campaign](../jwst-repeat-imaging-transients/README.md). That campaign measured extragalactic variability. This search uses its provenance conventions and covers stellar fields.
- Opened a separate branch for direct images of nearby stars. SPHEREx excesses and small-body astrometry are outside the selection.

## Pending

- Freeze the MAST inventory with explicit public status, release date, instrument and calibration coverage. Link each target to prior searches.
- Reduce a defined subset, record input hashes, measure throughput and injection recovery, and report sensitivity only inside the searchable footprint.
- For any point source, record detector and roll checks, astrometry against the background track, colours, contamination and prior reports. A single image cannot establish a bound companion.

## Observational checkpoint

The two-epoch F1500W experiment in the 4–12 arcsec annulus around WD 1202−232 and WD 2105−82 yielded 16 first-epoch peaks. Nine matched the stationary frame and seven lacked a qualifying second-epoch detection. None satisfied the common-motion selection. [Source-by-source measurements](miri/measured/sources.csv) · [Injected-source recovery](miri/measured/contrast.csv)

An F444W source 3.20 arcsec from L 788−34 in newly public program 8826 survives a change of reference-star library. Its association and prior-report checks are in progress. No bound companion has been established.
