# Fast extragalactic transients in the November 2026 public alerts

**Opening status.** The selection below was fixed before examining candidate photometry. Survey coverage, candidate identities and an event rate are unmeasured at this commit.

## Frozen search

The observing window is [2026-11-01 00:00, 2026-11-17 00:00) UTC. Search the public Rubin alerts through ALeRCE (Fink/ANTARES as fallbacks) and the public ZTF alerts through ALeRCE, Fink or the University of Washington packet archive. One astrophysical position is one candidate across streams. Record each survey's time scale, broker watermarks, retrieved rows, paging/coverage and hashes. For any necessary subset, fix its boundaries before inspecting its candidates and report the actual selection exposure.

Select objects with their first available positive detection in the window. Require two independent, positive detections in one passband, each at signal-to-noise ≥5, separated by 0.02–1.00 observer-frame days inclusive, with an absolute change ≥0.30 mag, >1.00 mag/day and ≥3 combined magnitude errors. Either direction qualifies. Use all eligible pairs. A single steep pair does not establish that the event both rose and faded in a day. Record meaningful upper limits separately. A prior outburst more than 30 days before discovery excludes an isolated-event interpretation where history can be checked.

For survivors, verify image stamps, motion, repeated stationary detections, prior activity, Gaia/Legacy/PS1 stellar or host context, source association and distance. Check public multi-survey coverage, TNS, GCN and other high-energy/multimessenger notices. Unavailable checks remain open. A redshift or a gamma-ray nondetection alone cannot establish an orphan jet. Calculate recovery separately from observational counts. A sky or volumetric rate requires measured effective exposure.

## Completed

Read the [remarkable-findings inclusion bar](../../docs/REMARKABLE_FINDINGS.md) and the [existing Rubin investigation](../rubin-alerce-2026/README.md). This search concerns extragalactic rapid transients. The earlier binary minor-planet, blue stellar variability and white-dwarf analyses address different populations.

## Pending

Survey retrieval, image/host/history checks, novelty checks, and recovery measurement.
