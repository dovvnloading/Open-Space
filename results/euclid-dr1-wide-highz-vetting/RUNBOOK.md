# Full-release execution and validation runbook

This is one branch-scoped GitHub Actions job, not a repository-wide PR check.

- Source: IRSA S3 HATS wide MER `partition_info.csv`, frozen SHA-256 `f7940bcaed19b6f67295c5592cc81783d798c0f858062d7293e5e3dfbfa57c2a`; validated index has **768 HATS partitions and 888,152,165 source rows**.
- Each run streams full Parquet row groups with projected photometry, processes 24 partitions per checkpoint, saves object-key, ETag, row count, row group count, selection counts, source matches and a complete strict-selected gzip JSONL register.
- Resume is based on per-batch SHA-256, not unverified progress counters. No skipped or failed partitions count toward a full-survey result.
- Only `research/euclid-dr1-wide-highz-vetting` pushes modifying this result folder (excluding generated `output/**`) trigger the workflow. Commits of generated outputs use `github-actions[bot]` and never trigger the job again.
- At scan completion, `aggregate.py` verifies 888,152,165 read rows and exact partition coverage; `vet.py` records conditional quasar-versus-dwarf colour likelihoods (only if control support is sufficient), motion, IR, published object types, spectrum metadata and image *metadata*, leaving pixel-level artifact inspection explicitly pending.
- No footprint-wide recovery fraction, calibrated QSO posterior, image-artifact clearance or spectroscopic discovery is asserted without supporting evidence. Inspect `output/status.json`, `output/summary.json`, `output/vetting_status.json` and the referee report.
