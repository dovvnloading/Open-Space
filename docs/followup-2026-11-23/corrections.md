# Catalog corrections: evidence and follow-up packets

Checked on **23 November 2026**, against Open-Space commit
`972cd5169eea2b83d7407df84eb60c167b3c084e`.

The strongest immediate actions are a narrowly scoped NICER timestamp report,
a replacement coordinate index for the SMGPS images, and a warning about the
unsupported NGC 1261b radial-velocity track. The numerical replays support
these reports. The more ambitious interpretations—identifying the physical
cause of the NICER flashes or deciding whether NGC 1261b exists—remain open.

## Work completed in this follow-up

| Check | Result | Evidence and limits |
|---|---|---|
| NICER detector packet, geometry, catalog count summary, blank-sky summary, year strata and frame controls | All six compact replays completed successfully | Reprocessing the committed tables; the original archive photons were not downloaded again in this follow-up. |
| SMGPS coordinate reconstruction | All 693 centers reproduced; output CSV byte-identical | Reused the native FITS header measurements collected on 15 November. |
| Independent SMGPS WCSLIB calculation | Maximum difference 0.000025 arcsec | Agreement between coordinate implementations on the same headers; not an absolute astrometric calibration of the sky. |
| Fresh SMGPS image-index retrieval | All 693 filenames, indexed coordinates and image sizes match the frozen header table | The Harvard VizieR mirror still served the affected values on this date. |
| NGC 1261b primary-source crossmatch | 24,540 atlas rows, 87 streams, 129 Ji targets, 59 Simon targets; 12 branch-b matches | Fresh primary downloads passed the embedded SHA-256 checks. Regenerated `matches.csv` and `coverage.csv` are byte-identical. |
| Current upstream Galstreams track | Same bytes as the pinned audit; 690 of 1,044 velocity samples exceed the speed of light in absolute value | Fresh GitHub retrieval, independently parsed and hashed; all 690 are outside the measured stars' RA interval. |
| NGC independent coverage | 70 Gaia IDs with zero DR3 RVS velocities; zero matches in 31,777 S⁵ DR1 rows; one match in the later 115-star Awad table | Fresh archive retrievals reproduced the earlier coverage result. Other archives have not been exhausted. |
| Spectroscopy candidate manifest | 58 retained atlas candidates; one later Awad spectrum explicitly flagged | Complete selection after removing the 12 exact Reticulum II matches. These are not 58 established stream members. |

The [replay log](corrections-replay.json) contains commands and outputs.
The [live-check receipt](corrections-live-checks.json) records source hashes,
versions and comparison counts; the
[fresh SMGPS index response](corrections-live-smgps-index.tsv) preserves the
server response. Source-level reports remain in
[NICER](../../results/nicer-magnetar-microbursts/README.md),
[SMGPS](../../results/smgps-image-astrometry/README.md), and
[NGC 1261b](../../results/ngc1261b-reticulum-ii/README.md).

## 1. NICER: request a correction to the published absolute timestamps

**Question and observing unit.** Do the published calendar strings identify
the same detector events and independent spacecraft detections when read in
their stated time system? The units are individual catalog rows within named
NICER observations, with shared burst storms treated as dependent samples.

The compact-table replay recovers 1,132 electronic-table rows in 159
observations. Of 944 pairs with coverage at both tested times, shifting the
predicted local window forward by 69.184 seconds increases the count in
843, ties in one and decreases it in 100. The median difference is positive
in 120 of 121 eligible observations. Excluding **all SGR 1935+2154 rows**
leaves 145 wins and one loss among 146 pairs. That exclusion differs from
excluding only its dense April 2020 observation.

The existing independent-spacecraft comparison is the strongest additional
timing evidence: after excluding the April 2020 storm, the +69.184-second
prediction matches seven Fermi bursts in seven visits within 0.25 seconds;
each of the three control offsets has zero such matches. This follow-up
replayed the compact NICER summaries, not the original Fermi extraction.
See the [all-catalog comparison](../../results/nicer-magnetar-microbursts/clock-validation/README.md#all-catalogue-gbm-comparison).

**Exact scope for the correction report.** The affected material is the
absolute calendar labeling of Chu et al. Table A1 / the machine-readable
`TDB Start` values and the two Figure 3 zero points. For the two reconstructed
microsecond events, the published and recovered barycentric labels are:

| Event | Published label, TDB | Recovered label, TDB, from the existing precision reconstruction |
|---|---|---|
| 2018 flash | 2018-11-22 11:58:14.425441 | 2018-11-22 11:59:23.609441 |
| 2022 flash | 2022-07-24 23:59:07.273212 | 2022-07-25 00:00:16.457212 |

The displayed digits preserve the reconstruction's numerical output; they
are not a new claim of absolute microsecond accuracy. A constant correction
does not change a burst duration. It is also not a conversion directly to
spacecraft UTC. NASA distinguishes the local calibrated TT event clock from
solar-system-barycentric TDB, and requires the applicable reference keywords
and `TIMEZERO` in conversion. [NASA timing definitions](https://heasarc.gsfc.nasa.gov/docs/nicer/analysis_threads/time/)

**Concrete report text for author review.**

> We recover the two Figure 3 event clusters in the public NICER products,
> but the existing independent barycentric reconstruction places their
> absolute TDB epochs 69.184 seconds after the printed labels. The same fixed
> shift improves 843 of 944 fully covered catalog windows, and independent
> Fermi comparisons outside the April 2020 storm favor the shifted times.
> Please check the serialization of Table A1 and Figure 3 from the
> post-barycorr products, including the time system and reference location.
> We can supply the row-level window table, source hashes and comparison
> code. The evidence concerns published absolute times; it does not identify
> the exact author-pipeline step or imply an error in NICER's hardware clock.

The report should include [the count table](../../results/nicer-magnetar-microbursts/catalog-clock/windows.csv),
[the independent timing comparisons](../../results/nicer-magnetar-microbursts/clock-validation/README.md),
and [the competing-convention analysis](../../results/nicer-magnetar-microbursts/timescale-forensics/README.md).
Request the actual post-`barycorr` EVENTS/GTI headers, source direction,
ephemeris and the calendar-formatting code for a named row. A same-event,
same-frame reconstruction at the unshifted epoch would challenge the
proposed correction. Uniformly adding a constant to every table row before
that review would conceal rows with ambiguous coverage or event identity.

### Keep the flash-origin question separate

The detector replay confirms zero fast-chain triggers among seven retained
2018 events and 32 retained 2022 events, versus 131/140 and 555/640 matched
controls. Five simultaneous 2018 hits satisfy the particle-filter criterion.
The 2022 proposed magnetar is 13.17 arcmin off axis. These are substantial
reasons to revisit the source assignments. They do not establish a particular
physical disturbance. The blank-sky summary contains **no** events at the
seven- or 32-hit multiplicities in 3,231,475.848 seconds; that negative control
must accompany an instrumental interpretation.

The existing Shao et al. ten-event comparison also needs to travel with the
packet. Its approximately 61.3-second discrepancy compares potentially
different spatial references. Replaying the frame bookkeeping reproduces
the discrepancy, but supplies no third clock and no diagnosis of an error in
Shao et al. The decisive additional work is an author-assisted, same-frame
reconstruction and, separately, a detector response/pile-up investigation
for the flashes. [Forensics and falsifiers](../../results/nicer-magnetar-microbursts/timescale-forensics/README.md)

## 2. SMGPS: repair the image index and resolve the duplicate row

**Question and observing unit.** Does each index position describe the
geometric center of the FITS file named in that same row? The unit is one
image, joined by its exact filename.

The standard-library reconstruction and independent WCSLIB calculation
agree for all 693 images. The median displacement of the old index from the
native WCS center is 88.776 arcmin; the maximum is 179.97642 arcmin. The old
indexed center lies outside its own image for 649 files. Only 20 of the 693
old centers fall within five arcminutes of the corresponding native center.
The index is consistent with a reflected linear longitude-offset formula
to 0.02545 arcsec, identifying the numerical pattern without identifying the
software that generated it.

The fresh Harvard VizieR response still contains exactly the same coordinates
and sizes as the frozen table. This establishes persistence at that mirror;
it does not establish the update state of every CDS service or mirror.
The native FITS headers and full images were not all re-downloaded today.

**Concrete repair request for the catalog maintainers.**

> Please regenerate `J/A+A/693/A247/list` GLON and GLAT from the native FITS
> WCS at each file's geometric pixel center, joining by exact filename.
> Our replacement `coordinates.csv` covers all 693 files and agrees with
> WCSLIB within 0.000025 arcsec after CSV rounding. In the indexed values,
> 649 centers fall outside their own images. This request concerns the image
> list; the native radio images and the five source-coordinate tables have
> separate scope.

Attach [headers.csv](../../results/smgps-image-astrometry/headers.csv),
[coordinates.csv](../../results/smgps-image-astrometry/coordinates.csv),
[audit.py](../../results/smgps-image-astrometry/audit.py),
[check_wcslib.py](../../results/smgps-image-astrometry/check_wcslib.py),
and this follow-up's live-index receipt. Use the full spherical native WCS,
not a replacement sign alone. The geometric center in FITS one-based pixels
is `((NAXIS1+1)/2, (NAXIS2+1)/2)`; software using zero-based centers must
adjust the convention explicitly.

Acceptance checks are all 693 filename joins, agreement with WCSLIB at the
chosen rounding tolerance, and successful self-retrieval at the native
centers after the corrected index is published. The existing three cone
examples supply fixed regression positions. No new live cone search was
needed to establish that today's full index values match the audited ones.

Keep the candidate-count question in a separate paragraph to the authors:
Table A.4 records 152 and 153 duplicate `G354.353−00.018` in the original
manuscript as well as the released table. The existing audit finds 236
distinct scientific records among 237 rows. Ask whether a duplicate should
be removed or a different candidate was omitted; do not silently invent a
replacement row or revise the authors' intended total. The original paper
and source remain [Anderson et al.](https://doi.org/10.1051/0004-6361/202451038)
and [arXiv:2409.16607](https://arxiv.org/abs/2409.16607).

## 3. NGC 1261b: separate the contaminated velocity sample from the stream test

**Question and observing unit.** Are the 12 original branch-b velocity stars
independent evidence of NGC 1261 debris, and does the distributed velocity
track stay within its data support? The units are exact Gaia source IDs and
individual tabulated track samples.

Fresh downloads of the Ibata and Ji machine-readable tables and Simon author
PDF passed their pinned hashes. The full reconstruction yields the same 12
matches and all-stream coverage CSVs byte for byte. All 12 original branch-b
velocity stars are firm Reticulum II members in at least one of the two
spectroscopic assessments. Removing those identifications leaves 58 atlas
candidates and zero original-atlas velocities for branch b. This does not
change the original Galactic mass fit: Ibata et al. excluded the NGC 1261
family from its 29-stream fitting sample. [Original atlas and scope](../../results/ngc1261b-reticulum-ii/README.md)

### A live upstream defect remains reproducible

The current `cmateu/galstreams` default-branch file retrieved today has Git
blob `870458a736dea9213b0e63d0837c2d8d7e016fb9` and SHA-256
`03a116fa782fc7c9a59ec31daff0c0bcc8758e5a65461755a89c27ffdf7f3a31`.
It is byte-identical to the pinned audit. The path's latest commit returned
by GitHub is `5f3125fec686c1aa199c374135010e2595ba56ed`, dated 6 August 2024.

Independent CSV parsing and the original audit both find 1,044 samples,
only 28 within the measured stars' RA interval, and 690 samples with
`abs(radial_velocity) > 299792.458 km/s`, all outside that interval.
The extrema are −32,929,072.1608 and +396,486.0309 km/s. This is a derived
track failure; it is not an observation of superluminal stars.
[Upstream file](https://github.com/cmateu/galstreams/blob/main/galstreams/tracks/track.st.NGC1261b.ibata2024.ecsv)

**Concrete report text for Galstreams and the atlas authors.**

> The NGC1261b Ibata2024 track contains 690 radial-velocity samples with
> magnitudes above the speed of light, outside the RA interval of its 12
> supporting velocity stars. Those 12 stars are exact Gaia-ID matches to
> published Reticulum II members or candidate members. Please flag the
> branch's distributed velocity track as unsuitable for an NGC1261 debris
> velocity constraint pending a member-aware reconstruction. Preserve the
> sky/proper-motion track and the original stellar measurements with their
> provenance; a numerical clip of the extreme velocities would not repair
> the underlying support or membership issue.

The supporting attachments are [matches.csv](../../results/ngc1261b-reticulum-ii/matches.csv),
[coverage.csv](../../results/ngc1261b-reticulum-ii/coverage.csv),
[track.json](../../results/ngc1261b-reticulum-ii/track.json), the source
reproduction scripts, and the current-file receipt in this packet. No
third-party issue or correction request has been sent in this follow-up.

### Deliverable for spectroscopy follow-up

The new [58-target manifest](corrections-ngc1261b-targets.csv) keeps the exact
Gaia IDs, published J2016.5 coordinates, proper motions, extinction-corrected
photometry, and every atlas stream assignment for each retained star. Import
Gaia IDs as text; floating-point spreadsheet imports can change their digits.
Propagate positions to an acquisition epoch before constructing a telescope
target file. The supplied magnitudes are the atlas's extinction-corrected
values, not measured acquisition-band magnitudes.

Gaia DR3 `4734256791553917952` is flagged with its later Awad et al. spectrum:
85.706 ± 3.365 km/s, [Fe/H] −1.114 ± 0.173, and a reported NGC 1261 membership
probability of 0.946. It appears in atlas groups 12, 13 and 14, so this one
star cannot distinguish the branches. The fresh independent-coverage replay
reproduces this exact match. The other 57 targets lack a match in that checked
Awad table; this is not a claim that no other spectrum exists anywhere.
[Awad et al.](https://doi.org/10.1051/0004-6361/202451930)

Select follow-up positions across the proposed branch outside the dwarf and
across regions that distinguish the overlapping branches. Joint velocity,
metallicity and distance consistency with the cluster's tidal model is the
test. Removing the 12 known matches alone does not make the retained set
free of contamination. Preserve targets that fail the model and repeat
velocities where binaries could dominate a membership decision.

## Reproduction and external-search receipt

Compact replay commands are in [corrections-replay.json](corrections-replay.json).
The primary NGC replay used the existing commands below, with a new cache;
the cache now contains the hashed input bytes. The target-list builder
requires the same verified atlas bytes and introduces no network queries.

```sh
cd results/ngc1261b-reticulum-ii
python reproduce.py --cache /path/to/cache --output /path/to/rebuilt
python track_audit.py --cache /path/to/cache --output /path/to/rebuilt/track.json
python independent_coverage.py --cache /path/to/cache
# From the repository root:
python docs/followup-2026-11-23/corrections-build-targets.py --atlas /path/to/cache/ibata.txt
```

A bounded public-web search on 23 November 2026 used the literal queries
`"NICER" "ae5655" correction erratum`,
`"SMGPS" "693" "A247" correction coordinates`, and
`"NGC 1261b" "Reticulum" correction`, with a 365-day recency filter.
It returned original papers and catalog documentation and did not surface
a corresponding correction. This is not proof that no correction or private
author discussion exists. The fresh VizieR and GitHub data comparisons above
provide the more specific checks of the currently served products. Direct
web-tool opens of the Chu and Awad DOI pages failed; their primary data
products were accessed through the successful archive replays described
above. NASA's timing guide and the Anderson arXiv record were accessible.

The completed work is a reproducible evidence package and a concrete
spectroscopy candidate list. Author correspondence, publication and changes
to external catalogs remain subsequent actions.
