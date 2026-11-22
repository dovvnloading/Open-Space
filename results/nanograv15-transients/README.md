# Blind pulsar timing transient search

**Awaiting observations.** This branch contains a search protocol. It has produced no event measurements, exclusion limits, or evidence for a nondetection.

The discovery sample consists of the 68 pulsars in the NANOGrav 15-year timing release. Narrowband arrival times define the search. Wideband timing of the same observations tests dependence on how arrival times were measured. European Pulsar Timing Array DR2 supplies a separate telescope check where the pulsar and epoch overlap. A lack of overlapping observations is recorded as missing coverage.

The search covers individual frequency changes, compact-object timing shapes, chromatic events, and signals shared across the array. A timing waveform alone cannot establish which physical mechanism produced it. The preregistration records the discovery threshold, how to distinguish previously published events, and what counts as a sensitivity limit. [Search protocol](PROTOCOL.md)

The NANOGrav release is v2.1.0 and the EPTA release is DR2 v2.1. The release checksums below identify the expected inputs. The workflow records SHA-256 hashes of the bytes it actually downloads.

| Archive | Publisher MD5 | Release record |
| --- | --- | --- |
| NANOGrav timing files | `557d42dd8486a5f8272d90dec9b228a8` | [Zenodo 16051178](https://doi.org/10.5281/zenodo.16051178) |
| EPTA timing files | `307b2b99cd352876b07409899448af8b` | [Zenodo 8300645](https://doi.org/10.5281/zenodo.8300645) |

[Individual search and sensitivity](SINGLE_SEARCH.md) · [Common-signal search and sensitivity](ARRAY_SEARCH.md) · [Reproduction commands](RUN.md)
