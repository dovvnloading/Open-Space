# Prospective held-out secondary analysis

**3 distinct real triggers met the frozen rule in 3,345 evaluable ObsIDs.**

Groups 96–191 assigned 3,351 ObsIDs. Prior exposure excluded 4. 0 lacked ledgers and 2 had incomplete statistical coverage. The 24 equal-scope control catalogs retained 22 triggers in total. The minimum two-test-adjusted FAP is 0.08. Inference covers only the evaluable held-out observations.

| ObsID | Trigger | Seed score | k/24 | Survey FAP | Exact 95% interval for exceedance | Two-test adjusted |
|---:|---|---:|---:|---:|---|---:|
| 30892 | real/0 | 24.2995 | 0/24 | 0.040 | 0–0.1425 | 0.080 |
| 23037 | real/0 | 9.14129 | 1/24 | 0.080 | 0.001054–0.2112 | 0.160 |
| 23473 | real/0 | 7.69243 | 9/24 | 0.400 | 0.188–0.5941 | 0.800 |

[Prospective protocol](../../HELD_OUT.md) · [Immutable sample](selection.json) · [Coverage and exclusions](coverage.json) · [Every evaluated raw classification](trigger_audit.csv.gz) · [Control maxima](controls.json) · [Positions, photometry, input receipts and scores](survivors.json)

[Original producer ledgers and artifact checksums](original_ledgers.json.gz) · [Sealed evaluation](original_result.json)
The source manifest identifies every original worker ledger. The replay command verifies its archived bytes, including observations retained in pre-correction history.
