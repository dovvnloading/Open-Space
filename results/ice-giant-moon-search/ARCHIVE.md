# Public imaging inventory

MAST returned 1,492 candidate observation groups and 3,262 distinct file URIs for 2,067 detector exposures. The inventory covers observations starting on or after 2022-01-01 UTC and before 2026-11-17 UTC. Its 179.524 GB total includes different reductions of the same exposures. Public access and archive membership refer to the retrieval time in [inventory_query.json](manifests/inventory_query.json).

| Planet | Observation groups |
|---|---:|
| Uranus | 610 |
| Neptune | 882 |

| Instrument / product type | Groups or files | Stored GB |
|---|---:|---:|
| HST WFC3/UVIS groups | 1,427 | |
| JWST NIRCam image groups | 48 | |
| JWST MIRI image groups | 17 | |
| JWST NIRISS image groups | 0 | |
| JWST calibrated CAL files | 640 | 55.013 |
| JWST count-rate cubes RATEINTS | 640 | 93.084 |
| HST calibrated FLT files | 1,427 | 16.959 |
| HST charge-transfer-corrected FLC files | 555 | 14.468 |

The catalog query takes the union of two selections. One retrieves target names containing URAN or NEPTUN, case-insensitive. The other retrieves all public HST/JWST groups with catalog positions in the swept sky boxes and compares their midpoints to daily interpolated geocentric planet ephemerides. A 10 arcmin center-to-center tolerance includes the JWST observing baseline and instrument fields. Science reduction uses spacecraft-specific coordinates. Eight MIRI Triton and Triton-background groups entered through sky matching. Calibration program 8975 remains included because its science detector imaged Uranus. The selection removed 101 lamp, dark and bias groups carrying coincident telescope pointings.

The instrument filter accepts image-type groups for WFC3, NIRCam, NIRISS and MIRI. Planet-labeled backgrounds, masks, chip gaps and off-detector targets require footprint verification. Slitless spectra, IFU reconstructions and disk-projected OPAL maps are outside the photometric imaging sample. A serendipitous group missing a relevant target label and an accurate pointing, or an association whose reported observation interval is too broad to locate the planet at its midpoint, can escape this inventory. The observation table records target name, selection route, temporal midpoint separation, instrument, group ID, proposal and filter.

Every selected group maps to at least one retained exposure file. [Observations](manifests/observations.csv) and [products](manifests/products.csv) join through observation IDs in the pipe-separated parent list. Matching file names need no proprietary credentials. [Downloaded JWST inputs](manifests/downloaded_jwst.csv) gives complete-file SHA-256 hashes for the 82 early program 6379 and 8975 inputs acquired during this investigation. Subset and byte-range receipts in the individual reductions distinguish the actual measurements from catalog membership.

From this directory, with Python 3, NumPy and Astropy installed, rerun the inventory and retrieve any chosen data subset.

```sh
python archive.py inventory manifests --start 2022-01-01 --end 2026-11-17
python archive.py products manifests
python archive.py fetch manifests/products.csv /data/moons --pattern 'jw06379001001_*_nrcb1_cal.fits'
```

MAST can reprocess files. An unchanged URI does not guarantee unchanged bytes. The download script checks the cataloged length and records SHA-256 per file. The astronomical observation time uses the time recorded in the FITS file.

[MAST API services](https://mast.stsci.edu/api/v0/_services.html) · [JWST archive products](https://archive.stsci.edu/missions-and-data/jwst) · [OPAL program and source HST proposals](https://archive.stsci.edu/hlsp/opal)
