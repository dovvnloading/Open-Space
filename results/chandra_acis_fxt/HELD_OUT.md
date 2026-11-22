# Prospective secondary test for isolated high-latitude flares

The protocol was fixed at 06:35 UTC on 22 November 2026 for primary shard groups 96 through 191 inclusive. Those groups were queued and their outcomes had not been examined in designing this test. Earlier groups informed the hypothesis that a high-latitude flare on a quiet baseline warrants separate evaluation. This is the only secondary selection specified here.

**Sample.** Use zero-based shard assignments from commit `0655e5b54215a3d7bbac1f6b2561379bb50c81cf`. Remove every ObsID assigned to the moderate-count early cohort, the three-observation pilot, or the earlier source-ledger migration, including any ObsID assigned again for retry. Freeze the original assignment and exclusion bytes, hashes and resulting ObsID list. The 14 release-date-ambiguous observations are outside the dated census. [Fixed assignments](https://github.com/dovvnloading/Open-Space/blob/0655e5b54215a3d7bbac1f6b2561379bb50c81cf/results/chandra_acis_fxt/evidence/inventory/shards.json)

**Selection.** Apply all four rules to each raw real trigger and to each raw trigger in all 24 primary time-scramble replicas, before using the primary overlap-deduplication rule independently within each stream.

1. The original primary instrumental decision is a pass.
2. The recorded primary catalog disposition is `no_counterpart` or `host_candidate`.
3. The absolute catalog-vetting Galactic latitude is at least 30 degrees.
4. Both flank exposures are positive and the expected number of aperture photons in the trial at the measured outside-window rate is at most one. Use `trial_duration_s * (before_counts + after_counts) / (before_s + after_s)` and retain equality.

**Missing data.** Record every assigned ObsID and all exclusions. An ObsID contributes to inference only when every acquired product completed the frozen search, all 24 replicas are present, and every real and scrambled trigger can be classified by these four rules. A known failed condition resolves a trigger as a failure even when another field is missing. A missing required field or incomplete catalog query that could change membership leaves the trigger unresolved and excludes that ObsID from all 25 streams. Keep available rows and reasons from excluded observations in the audit. Missing real-only flux or diagnostics is reported separately and does not change the statistical selection. Preserve the primary ledgers intact.

**Statistic.** Rank survivors by the original seed score. For each real survivor, count the number `k` of the 24 held-out control-catalog maxima greater than or equal to its score. An empty control catalog never exceeds a real score. Report `k/24`, the empirical survey FAP `(k+1)/25`, the two-sided exact 95% Clopper–Pearson interval for the underlying control exceedance probability, and `min(1, 2*(k+1)/25)`. The factor of two covers the primary analysis and this one secondary analysis. The smallest attainable adjusted value is 0.08. These conditional results describe only the evaluable fixed held-out subset. The primary injections characterize instrumental recovery. This protocol makes no population-rate bound or automatic discovery claim.
