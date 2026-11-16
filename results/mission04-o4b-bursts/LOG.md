# Mission 04 investigation log

## 2026-11-16 · Investigation opened

The target is previously unrecognized coherent burst power in validated public O4b strain. A null search can constrain explicitly defined waveform and exposure classes only after measuring efficiency in observed noise. GWOSC released the full O4b data on May 26, 2026, including Hanford, Livingston and Virgo. The observing run starts at GPS 1396796418 and ends at GPS 1422118818.

The isolated branch was opened before inspecting any foreground strain. Work in progress includes the overlap audit, checksummed data selection, network statistic, controls, time shifts, sensitivity and independent checks. No on-source result has been established at this checkpoint.

[GWOSC O4b release](https://gwosc.org/O4/O4b/)

## 2026-11-16 · Before the first foreground scan

The initial protocol fixed ten DATA-selected frames, three durations, four bands, two polarities, 41 delays, 128 file-bounded time shifts, an empirical escalation criterion, the separate count threshold of 20 and the primary injection population. GWTC-5 public candidate rows and all public GraceDB production superevents supplemented the confident catalogue. The primary positive control was a separate, published event. [Prespecified decision record](https://github.com/dovvnloading/Open-Space/pull/114#issuecomment-6923893538) · [Input and control checkpoint](https://github.com/dovvnloading/Open-Space/pull/114#issuecomment-6924015857)

## 2026-11-16 · First ten epochs

All 32 input files passed their GWOSC checksums. The searched exposure was 40,315 s. Of 128 shift campaigns, 24 reached the highest unshifted score. No observed second exceeded the rate threshold. The 150 Hz, Q=9 primary signal had 1,792 recoveries in 2,048 injections and a pointwise 90% upper rate of 7.45/day. The loudest observed tile, the held-out GW250114 event and the ten monthly maxima were checked separately. [Foreground checkpoint](https://github.com/dovvnloading/Open-Space/pull/114#issuecomment-6924082866)

## 2026-11-16 · Disjoint extension and durable checkpoints

The first-stage exposure motivated the extension. Before inspecting the new foreground we fixed the calendar slots, 8,192 injection draws, coverage convention, and two-search multiplicity treatment. The 125 selected blocks and 43 empty slots were committed in d6a739f. Each finished month received its own reviewable archive containing every live foreground second, the delay-shifted tails, masks, and injections. Transient GWOSC download failures were retried, and the official OSDF origin could supply the same MD5-verified bytes. [Extension freeze](https://github.com/dovvnloading/Open-Space/pull/114#issuecomment-6924145706) · [Coverage checkpoint](https://github.com/dovvnloading/Open-Space/pull/114#issuecomment-6924313872)

## 2026-11-16 · Extension complete

All 125 scheduled frame pairs passed file validation. The final 499,071 live seconds contain no foreground crossing at the fixed threshold of 20. Of the 128 shifted campaigns, 118 reach the loudest foreground score. On exactly equal Hanford exposure, 106 shifted campaigns still reach it. Independent raw-strain and detector-context tests, complete shifted and unshifted ledgers, and all assigned injection draws are included. With 6,950/8,192 recoveries, the pointwise 90% rate upper bound is 0.601/day. Accounting jointly for the two primary campaign statements gives 0.741/day for the extension. The repository's 150 tests, five burst-specific regressions, bundled pilot and 12-seed benchmark passed. [Final numbers in the PR](https://github.com/dovvnloading/Open-Space/pull/114#issuecomment-6924626786) · [Report](README.md) · [Independent audit](results/audit_results.json)
