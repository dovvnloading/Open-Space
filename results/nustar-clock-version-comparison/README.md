# NuSTAR clock-version comparison: calibration feasibility

For the as-flown schedule envelope of Crab observation **11202302010**, three public NuSTAR fine-clock solutions differ by microseconds. The observation is an identifying calendar window for a calibration comparison; no photon data, pulse phase, clock accuracy, coherence improvement or astrophysical finding was measured. The original selection and calculation were frozen in [PROTOCOL.md](PROTOCOL.md) and [CALIBRATION_FREEZE.md](CALIBRATION_FREEZE.md). The results below are differences between solutions evaluated on the same dates, not independent observations.

| Later minus earlier version | Minimum (µs) | Maximum (µs) | Largest change from envelope start (µs) |
| --- | ---: | ---: | ---: |
| v235 − v234 | −4.776 | +3.854 | 8.630 |
| v236 − v235 | +6.728 | +13.236 | 6.508 |
| v236 − v234 | +8.460 | +10.582 | 2.122 |

The interval is 2026-10-10 02:00:01–07:50:00 UTC, 20,999 seconds. It is the SOC as-flown schedule envelope, **not** measured live time: it includes Earth occultation and other gaps. NUMASTER reports archived, public on 2026-10-16, with effective exposures 4,979.0385 s (A) and 5,205.2320 s (B); [selected.xml](selected.xml) preserves the query response and [obs_metadata.json](obs_metadata.json) its decoded values. The archive directory was accessible, and its [event listing](event_cl-index.html) shows A and B cleaned files. We did not download or inspect those events. A public orbit file was downloaded only as an access test; it was not used in this calculation.

## Two calculations

The [sampled implementation](sampled.py) evaluates cubic Hermite interpolation directly at 60 s spacing and the endpoint, converting UTC epochs to TT relative to each FITS reference. All 351 sample times have bracketing rows in all three files. The [independent implementation](extrema.py) uses SciPy's `CubicHermiteSpline`, forms differences on the union of knots, and checks all interval endpoints and real derivative roots. For all three pairs it found 22 polynomial segments and no interior derivative roots; the extrema are envelope endpoints. The two numerical results agree to better than 1e-9 µs. Both methods use the native `NU_FINE_CLOCK` header's convention: interpolate offset and fractional-frequency derivative and add the resulting correction to spacecraft time. They do not independently establish the accuracy of this convention's upstream calibration.

All versions bracket the full envelope and pass the header `TSTART`/`TSTOP` validity gate. Their 23 selected knots have a maximum 1,000 s separation. The bracketing-knot `CLOCK_ERR_CORR` columns range approximately 92.04–93.41 µs for v234, 90.45–91.95 µs for v235 and 90.48–92.07 µs for v236. These are table fields; no statistical independence or phase-uncertainty interpretation is assumed. The updated files are revisions of the same mission clock data, not independent timing experiments. In particular, v236 was generated on November 16 and announced November 17; it supersedes the v235 file generated November 2.

## Inputs and access

Each clock file below was downloaded from the official SOC `Clockfile_v2` directory. They are retained outside Git at `D:/Open-Space/campaigns/nustar-clock-access-20261117/` as `v234.fits.gz`, `v235.fits.gz`, and `v236.fits.gz`; the three exact files total approximately 38 MB. The durable links, lengths and SHA-256 hashes are:

| Version | Official input | Bytes | SHA-256 |
| --- | --- | ---: | --- |
| v234 | [nuCclock20100101v234.fits.gz](https://nustarsoc.caltech.edu/NuSTAR_Public/NuSTAROperationSite/Clockfile_v2/nuCclock20100101v234.fits.gz) | 12,651,240 | `fe872bf396be44c6e4ee20551e4fde89c9cff0809e16057b5a35bcc2b0015e8e` |
| v235 | [nuCclock20100101v235.fits.gz](https://nustarsoc.caltech.edu/NuSTAR_Public/NuSTAROperationSite/Clockfile_v2/nuCclock20100101v235.fits.gz) | 12,679,426 | `24a440096719c2f8c992d6c03e3c3ac769bc7f7638826b1a422d64a1a9d75d19` |
| v236 | [nuCclock20100101v236.fits.gz](https://nustarsoc.caltech.edu/NuSTAR_Public/NuSTAROperationSite/Clockfile_v2/nuCclock20100101v236.fits.gz) | 12,709,629 | `f30744159e4dcff43a525d15318bd4fb88e561c1fda9fcc508e0db752bb30f96` |

[NUMASTER](https://heasarc.gsfc.nasa.gov/W3Browse/nustar/numaster.html) query: `SELECT obsid,name,time,end_time,public_date,exposure_a,exposure_b,status FROM numaster WHERE obsid = '11202302010'` via HEASARC TAP `https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync`, `REQUEST=doQuery`, `LANG=ADQL`. [Archive directory](https://heasarc.gsfc.nasa.gov/FTP/nustar/data/obs/12/1/11202302010/) was HTTP 200; FPMA event HEAD returned 200 and 571,562,998 bytes, not a download. The [orbit product](https://heasarc.gsfc.nasa.gov/FTP/nustar/data/obs/12/1/11202302010/auxil/nu11202302010_orb.fits.gz) was 1,002,372 bytes with SHA-256 `afd3eb6ee1cdea3cdf792e13754d4bfccd9255c144f4de5b08d9d0472b1f14b9`. Preserve the local bytes or reacquire exact hash-matching bytes.

The independent [Jodrell Bank monthly Crab ephemeris](https://www.jb.man.ac.uk/pulsar/crab/crab2.txt) was retrieved and checked again on **2026-11-18 01:50:56 UTC**. Its last nonblank row is September 15, 2026 (MJD 61298); no October solution is present in the checked file. The [preserved bytes](crab2.txt) have SHA-256 `70dbc14986defcb00421d341e9554e369d0610bcb8fd83f0ca49a430d19685dd`. This is a negative lookup in one public product, not proof that no contemporaneous ephemeris exists. The September row's 100 µs accuracy is not assigned to October. An independently measured, time-valid phase reference and a demonstrated photon-level sensitivity are needed before testing coherence or absolute timing. The observed model drift of 8.63 µs is about 0.00026 rotations at a 33 ms period, but this scale is not a detection limit.

## Reproduce

Install Python with Astropy, NumPy and SciPy into an isolated environment. Place the three hash-matching files in one directory, named `v234.fits.gz`, `v235.fits.gz`, `v236.fits.gz`; set `NUSTAR_CLOCK_INPUTS` to that directory and `NUSTAR_CLOCK_OUTPUT` to a new empty directory, then run:

```powershell
python results/nustar-clock-version-comparison/sampled.py
python results/nustar-clock-version-comparison/extrema.py
```

The extrema script verifies each input SHA-256, time references, monotonicity and stated coverage, and refuses to overwrite its output. The sampled script also verifies all three input hashes before reading FITS tables. The reproducibility output [sampled.json](sampled.json) and [extrema.json](extrema.json) used Python 3.14, NumPy 2.5.5, Astropy 8.0.3 and SciPy 1.18.1; [original independently generated extrema](extrema-original.json) used Astropy 8.0.2 and SciPy 1.17.1 and agrees. The revised scripts and packet files are frozen by [SHA256SUMS](SHA256SUMS) and the Git commit; editing a script changes that identity.

Sources: [NuSTAR clock notes](https://nustarsoc.caltech.edu/NuSTAR_Public/NuSTAROperationSite/clockfile.php), [NASA NuSTAR updates](https://heasarc.gsfc.nasa.gov/docs/nustar/news/nustar_whatsnew.html), [NuSTAR as-flown timeline](https://nustarsoc.caltech.edu/NuSTAR_Public/NuSTAROperationSite/AFT_Public.php), and [barycorr documentation](https://heasarc.gsfc.nasa.gov/docs/software/lheasoft/help/barycorr.html). For observations on or after 2026-09-24, NASA specifies NuSTARDAS 2.1.6 / HEASoft 6.37.1; no event analysis or barycentering was performed here.
