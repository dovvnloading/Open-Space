# Contributing to Open-Space

The current research program investigates unidentified narrow optical features in SDSS DR17 spectra. Contributions may add a reproducible observation, test a known explanation, measure a failure mode or improve the tools an agent can use. A rejected hypothesis with inspectable evidence is a research result.

## Submit a research result

Open a pull request containing the question, dataset and outcome. For observational work, include the observation ID, archive release, source URL, SHA-256 digests, wavelength convention, code revision, command and eligibility criteria. Report how many positions or targets were searched and how many had adequate coverage. Include a plot or machine-readable evidence packet, the strongest competing explanation and a falsifiable next test. For a literature check, record queries, dates, databases and primary references. Cite original data providers and retain their attribution and terms.

Keep generated ledgers, model credentials, downloaded archive batches and unrelated raw logs out of Git. Check in compact derived results and the recipe or pinned inputs needed to regenerate them. Explicitly identify any synthetic or model-authored content. Reviewers should be able to disagree with an interpretation while reproducing its measurements.

## Research projects ready to take on

| Project | First experiment | A reviewable result |
| --- | --- | --- |
| Independent sky controls | Select sky-fiber controls and an exposure-difference sample on held-out nights. Apply the frozen detector to each. | Eligibility denominators, per-spectrum false-trigger counts and bounds, with observing conditions recorded. |
| Astrophysical identification | Search an unassigned rest wavelength against a versioned line list using consistent vacuum units. Predict a companion before inspecting that window. | Query, provenance, blend alternatives and a held-out measurement or upper limit. |
| Repeated observations | Crossmatch one flagged target to another visit. Inspect exposure IDs and record the second spectrum's sensitivity. | Frozen prediction and a confirmation, rejection or inconclusive test on independent photons. |
| Improve the detector | Compare a resolution-aware template family with the shipped baseline on disjoint targets. | Curves for recovery and false alarms at matched operating points, including artifacts and missed signals. |

For code changes, run the unit tests and bundled pilot. Changes to masking, reference tables, continuum estimation or thresholds should explain how their search coverage and selection effects change. Keep the prior measurement and compare it with the new version on an identical input. The local ledger keys experiments by input, model configuration and code fingerprint.

For additional archives, choose one instrument and release. Document wavelength and velocity conventions, pixel masks, flux units, exposure relationships, data rights and how to recreate the selected inputs. Demonstrate a measured known feature and a usable negative control before interpreting unidentified entries. [Observational contract](docs/DATA.md) · [Science protocol](docs/SCIENCE.md)
