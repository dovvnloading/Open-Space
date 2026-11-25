# Formal mathematical verification (Lean 4)

This is an initial, deliberately narrow complement to the [Python verifier](../docs/VERIFICATION.md).
It does not replace observational checks, archive provenance, numerical replay or scientific review.

## Run

Install Lean 4 via Elan and Python 3.10 or later. From the repository root:

```sh
bash formal/check.sh
```

The pinned toolchain is in [lean-toolchain](lean-toolchain).
Only Lean's bundled library is required; no mathlib or astronomical downloads.
The script verifies the pinned Lean version, runs adapter unit tests, generates an ephemeral
Lean certificate from the checkout's CSVs, compiles [Foundations.lean](OpenSpace/Foundations.lean),
and checks the certificate with warnings as errors. An [axiom audit](audit_axioms.py) enforces
that all exported claim reports are present, that exact certificates use no axioms and that
general theorems use only the standard Lean foundations. No `sorry`, custom axioms or
`native_decide` are used. [Negative controls](test_negative.py) alter a temporary copy of
each of the three input tables and require Lean to reject the false proposition.
A successful run is the criterion for calling the stated Lean proofs checked.

## Claim-to-proof map

| Source | Formal statement | Boundary |
| --- | --- | --- |
| [NICER catalog clock](../results/nicer-magnetar-microbursts/catalog-clock/README.md), [window counts](../results/nicer-magnetar-microbursts/catalog-clock/windows.csv) | `nicerClockOffset`, `nicer2018PrintedLabels`, `nicer2022PrintedLabels`, `commonShiftPreservesDuration`, `commonShiftPreservesWindow`, `pairedOutcomePartition`, generated `nicerCatalogExact` | The TT−TAI and TAI−UTC constants and rounded barycentric labels are inputs; neither the orbit correction nor timing pipeline is verified. The generated certificate counts the committed narrow-window rows, treating blanks as ineligible, including the storm-excluded subset and event totals. It does not establish independence, a p-value, an instrumental origin or a source association. |
| [Venus occultation forecast](../results/solar-orbiter-venus-jovian-occultations/README.md), [contact table](../results/solar-orbiter-venus-jovian-occultations/contacts.csv) | `orderedOverlap`, `orderedCoveragePartition`, `largerLimbPreservesCoverage`, generated `decemberContactLabels`, `decemberContactOrder`, `decemberContactDurations`, `decemberDurationRounding` | Exact arithmetic on rounded same-day UTC labels, assuming continuous intervals in the selected model. No claim of certified ephemeris interpolation, root completeness, angular geometry, navigation errors, future UTC or radio detectability. The larger-limb lemma assumes the modeled angular radius has increased. |
| [OEC orbital-scale audit](../results/oec-orbital-scale-audit/README.md), [primary inputs](../results/oec-orbital-scale-audit/primary_checks.csv) | `stellarRadiusScaling`, generated `oecAxis0` through `oecAxis5` | Under the adopted central-value stellar solution and nominal solar radius, exact rational inequalities (via scaled naturals) check the 1.5% agreement of the erroneous transform, rounding of the six corrected axes to five decimals, factors strictly between 4 and 5 and fixed-luminosity flux ratios strictly between 17 and 23. No new orbital fit, uncertainty propagation, catalog priority or physical model validation. |
| [Spectral injection recovery](../open_space/calibrate.py) | `calibrationAccounting`, `noAssessableIff` | Conditional bookkeeping when recovered, missed and confounded exhaust attempted trials. Does not verify that the Python classifier implements the partition or calibrate significance. |

The two time label proofs use microseconds from midnight (adding one day for the second
2022 label); the December certificate uses integer milliseconds from the same midnight.
Duration agreement allows 6 ms for the two endpoint roundings and separately rounded
centisecond duration. This is an arithmetic rounding bound, not timing uncertainty.

## Evidence trust and extension

The Python [generator](generate.py) is an explicitly trusted input adapter. It parses CSV
and emits only numeric literals and fixed Lean proof templates. For NICER, each small
block's proposed summary is checked in Lean against its encoded rows; the aggregate
target is the published summary, not a Python recomputation used as the target.
The Python adapter's row selection, units, transcription and the underlying data remain
outside the kernel. Decimal scaling uses integer arithmetic, rejecting unrepresentable
values rather than rounding; the adapter rejects malformed tables and duplicate
(observation, burst number) keys. A change to a selected CSV regenerates the certificate.
The unit tests and three negative controls help catch transcription and verification regressions;
they do not formally prove the adapter or guarantee detection of every conceivable data change.
The initial source revision was `d69fae651c458ca0c22e7f9fb4a85dbc0aa24d70`.
The three input Git blob IDs at that revision were, respectively,
`0d4e703553eac3535400773c0a5a0c64a4443302`,
`c7ef0757862abbb1f7bbde2352a965208160fa02`, and
`b7a8aefcb989483a07a2b15f0789c51c967d3fbb`.

For new work, write explicit units, domain assumptions and theorem statements linked to
specific source files. Keep exact arithmetic, model assumptions, interval-certified
numerical bounds, statistical inference and observations separate. Prefer kernel-checked
proofs and constructive counterexamples; do not mark conjectures proved or bypass missing
proofs with axioms. Next candidates include barycentric and spherical-vector geometry,
redshift and unit transformations, least-squares identities, multiple-testing bounds,
and interval-enclosed ephemeris roots. They require additional modeling or mathlib work
and are not part of this initial certificate.

The [initial Lean 4.19.0 run](https://github.com/dovvnloading/Open-Space/actions/runs/45648185212)
passed on the first implementation commit. Consult the latest pull-request or main-branch
workflow for the current revision's status; no whole discovery has been certified by Lean.
