# Juno UVS impact control: independent new investigation

Status: published flash recovered in photon time/geometry and longwave counts. See [measured control and limits](VERIFIED_CONTROL.md). Not a new impact discovery or full spectral-model reproduction.

Question: can photon-level Juno UVS data recover the published 2020-04-10 bolide before a frozen search of previously unexamined post-2020 observations?

Observing unit: one calibrated photon-list acquisition and its spacecraft scan; neighboring time bins are not independent events.

Positive control: Giles et al., https://doi.org/10.1029/2020GL091797. Instrument guide: https://pds-atmospheres.nmsu.edu/data_and_services/atmospheres_data/JUNO/uvs.html. H2 transient confounder: https://doi.org/10.1029/2020JE006659.

The selected manifest pins public product identifiers, sizes and the source index hash. Source byte provenance is in acquisition.json; raw photon bytes are fetched by the replay script rather than committed. Files on the control day are approximately 0.7–1.3 GB each; the frozen four-second replay downloads only 9.43 MB. Archive release dates are not observation dates.

Falsification gates:
1. Recover the published flash time and location from real calibrated photons.
2. Reject detector-coordinate fixed artifacts, particle contamination, masked events and scan/exposure changes using housekeeping and off-disc controls.
3. Distinguish thermal continuum from H2 auroral/lightning-related transients.
4. Freeze detector thresholds before screening new passes; measure injection recovery and shuffled/null global false alarms.
5. Review literature and independent evidence before calling any candidate new.

A failed reproduction blocks population claims. A candidate count without exposure and completeness cannot establish an impact rate. No new impact or new physics is claimed here.

