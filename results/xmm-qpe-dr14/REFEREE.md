# Independent referee: 4XMM-DR14 QPE pilot (19 November 2026)

**Verdict: NOT READY FOR SCIENTIFIC MERGE.** The Actions execution succeeded, the catalogue receipt's published MD5 matches, and the sample was frozen before PPS retrieval. The proposed QPE result **fails the positive-control and completeness bars**. This is an inconclusive engineering pilot only.

## Rechecked facts

1. The [committed summary](outputs/summary.json) reports 1,035,832 DR14 detections, 372,313 TSERIES, 33,488 soft/long prefiltered detections, 12 PPS requests, 9 with parsed curves, zero unvetted repeated-peak flags. The [frozen manifest](outputs/frozen_selection.json) gives eight ranked targets and four control requests. The data were downloaded from release-tagged HEASARC/NXSA links as shown in [receipts](outputs/receipts.json).
2. **Known positive control failure:** GSN 069, RX J1301.9+2747 and eRO-QPE2 were rejected by a hard 120 MiB PPS archive cap; eRO-QPE1 was parsed for 25.264 h but the flagger reported zero peaks. eRO-QPE3 was matched to the catalogue but not requested. AT2019qiz had no match to the supplied coordinate; Ansky and eRO-QPE4 were absent from the control table. Thus **no published QPE flare was recovered**.
3. **No candidate nucleus:** none of the eight ranked light curves triggered a repeat flag, and no host-galaxy central coincidence or stellar rejection was attempted. It would be misleading to call these eight galactic nuclei or to infer a DR14-wide non-detection.
4. **No QPE-physics measurement:** the output lacks flare intervals, amplitudes and thermal evolution for any positive flag. The band-averaged catalogue hardness cuts do not establish spectral hardening or fitted blackbody kT. Background median `BACKV=0` does not constitute a calibrated, scaled source/background flare comparison; pile-up, GTIs, edge effects and burst significance were not validated.
5. **No calibrated false-positive or recovery fraction:** no synthetic injections into observed light curves, and no matched stellar/field/background controls. A zero flag count from a positive-control-failing detector is not a null result.
6. **Limited literature check:** targeted papers establish that GSN 069, RX J1301.9+2747, eRO-QPE1–4, AT2019qiz and Ansky are previously reported; there was no novel flagged object for which a targeted ATel/TNS/Swift/Chandra/eRASS1 check could be completed. No novelty was established.
7. **Reproducibility caution:** `frozen_selection.json` embeds non-finite `NaN` tokens (not valid standards-compliant JSON); `frozen_detections.csv`'s SHA-256 is recorded in the manifest but was not independently redownloaded/rehash-verified by this review. PPS tarballs are hashed in receipts but do not have independently supplied publisher digests. The selected sample is biased to extremely low chi-square variability probabilities, with many values floored at 1e-15, and is not a statistically characterized parent-sample search.

## Acceptance conditions not satisfied

Demonstrate recovered **actual** published QPE flare trains with source and scaled-background significance and energy/temperature evolution; fix control retrieval and nonstandard JSON; validate nuclear hosts, assess stellar/background false positives and pile-up, perform injection recovery over a stated denominator, and screen any surviving sources against prior QPE/TDE/AGN reports and other X-ray epochs. Until then, **keep PR #205 draft; do not label [REMARKABLE], claim a new QPE, claim an archive-wide null, or merge as a validated scientific result.**

Evidence: [Actions run #44470633598](https://github.com/dovvnloading/Open-Space/actions/runs/44470633598) · [Outputs](outputs/) · [Remarkable-findings policy](../../docs/REMARKABLE_FINDINGS.md).
