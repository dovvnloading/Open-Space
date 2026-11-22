# Recovering a large Chandra evidence deposit

Use this procedure in the existing Chandra workflow if the ordinary deposit exhausts runner capacity or a Git push fails. The primary search, its thresholds and the existing shard run continue unchanged.

## Measure first

Start from a fresh checkout of the research branch. Record the Actions run ID, attempt, immutable source commit, artifact IDs and service-reported artifact SHA-256 digests. Keep the original artifacts available until independent replay and the final committed-tree verification succeed. Record `df -h`, `df -i`, `free -h`, `swapon --show`, the cgroup memory limit, swap limit and OOM counters, and `git count-objects -vH`. Measure a Git tree by summing uncompressed blob sizes **by path**, including repeated object IDs. Use `git ls-tree -r -l -z HEAD` and parse records as NUL-delimited bytes.

The repository verifier reads every committed blob through one `git cat-file --batch` buffer and copies the blobs into a dictionary. Its integrity pass hashes one level of ZIP members of up to 200,000,000 bytes each. It retains their digests and locations. A file ending in `.gz` is fully decompressed without that member-size cap. The trace pass decompresses eligible `.json.gz` and tabular evidence before its 60,000,000-byte size check. It retains decoded text and numerical indexes across result folders. Splitting one large Git blob limits individual file size but leaves the aggregate tree memory requirement. [Git tree reader](https://github.com/dovvnloading/Open-Space/blob/0655e5b54215a3d7bbac1f6b2561379bb50c81cf/open_space/verify/tree.py) · [Integrity pass](https://github.com/dovvnloading/Open-Space/blob/0655e5b54215a3d7bbac1f6b2561379bb50c81cf/open_space/verify/integrity.py) · [Trace evidence reader](https://github.com/dovvnloading/Open-Space/blob/0655e5b54215a3d7bbac1f6b2561379bb50c81cf/open_space/verify/evidence.py)

## Preserve all catalog bytes

From the repository root, use the frozen survey dependencies and select the exact evidence directory. The measurement command builds and verifies proposed packs in temporary storage and prints the actual native, retained-loose and compressed sizes.

```bash
survey=results/chandra_acis_fxt
python "$survey/catalog_pack.py" --evidence "$survey/evidence" measure
python "$survey/catalog_pack.py" --evidence "$survey/evidence" pack --prune
python "$survey/catalog_pack.py" --evidence "$survey/evidence" verify
```

Each producer under `evidence/catalogs/` has a separate, content-addressed set of ordinary ZIPs. Every original request memo remains readable at its original path, including error and retry receipts. The helper also preserves loose replies for every real instrumental pass and direct catalog-response links found in the survey Markdown. ObsID ledgers, diagnostics, candidate photometry, independent visits, prior history and the supplemental search are untouched. Missing or conflicting source receipts stop the operation.

The packer sorts exact logical paths, uses fixed ZIP metadata, and records the complete SHA-256 and byte length of every source and every part. Each ZIP holds at most 16 MiB of source chunks, at most 512 response members, and less than 24 MiB of ZIP bytes. Every member, including the inventory, is at most 8 MiB. Small replies retain their original paths inside the ZIP. Larger replies have ordered parts. Each ZIP contains a `MEMBERS.json`; each producer’s `index.json` records pack and logical-tree checksums, linked loose files, and a separate checksum over the untouched request-memo tree. ZIPs are one level deep. The packer rereads every byte and checks each reconstructed original before atomically publishing the new index and pruning verified loose duplicates. An interrupted run can repeat the same command. With the same Python/zlib runtime and inputs, repeated packs have identical bytes.

For an artifact-by-artifact recovery, securely unpack one verified producer artifact to temporary storage. Run the helper against that artifact’s evidence root. Transfer its ledgers, diagnostics, unchanged receipts, retained replies, ZIPs and index to the branch checkout. Compare SHA-256 when a previous artifact supplied the same path and stop on conflicting bytes. Keep the original working directory until the transferred index passes verification in the destination. Release that temporary disk space before acquiring the next artifact. Pass `--cache CACHE` to measure, pack or verify just one namespace. After all transfers, rerun packing once against the assembled store so that new written-review links are materialized at their original paths.

## Read or replay

Each ZIP can be downloaded directly from its producer directory on GitHub. The member inventory maps exact original paths to response bytes or ordered parts. The following commands locate one reply, restore a selected cache into a separate directory, and run the original observation replay against a temporary reconstructed catalog view.

```bash
python "$survey/catalog_pack.py" where catalogs/CACHE/RESPONSE.bin
python "$survey/catalog_pack.py" extract --cache CACHE --destination "$RUNNER_TEMP/one-cache"
python "$survey/catalog_pack.py" extract --obsid OBSID --destination "$RUNNER_TEMP/one-observation"
python "$survey/catalog_pack.py" replay OBSID
```

The selective NED repair reads packed catalogs directly. Its unchanged `followup.py prepare`, worker and `finish` commands verify each recorded original reply against its byte count and SHA-256, including HTTP attempts that retained a body. An affected event’s worker copies only its own responses into private staging at the paths expected by the frozen reader. The final pass checks carried observations against primary evidence and repaired observations against the primary and worker evidence together. The repair plan pins the catalog reader source. Use `extract --all --destination TEMP` only for a separate consumer that requires every response in loose form. Run the survey replay and both repository checks in the Chandra workflow.

## Push and verify the final commit

For disk pressure, process one artifact or namespace at a time, drop already-verified staging files and keep calibration installations out of the final verifier job. For memory pressure, measure the failed phase with `/usr/bin/time -v` and inspect cgroup OOM counters. Use a sufficiently large runner, or runner-local swap if permitted and disk capacity covers it. Catalog packing does not cap the numerical index built from the retained observation ledgers.

GitHub permits at most 2 GiB in one push. Stage moderate, disjoint groups of content-addressed ZIPs and evidence, then commit and fast-forward push each group on the same research branch. A sequence of local commits sent together still constitutes one push. Keep publication metadata and the root catalog change in the final commit after their evidence objects exist remotely. A practical starting point is several hundred MiB of newly introduced bytes per push; adapt after measuring. A local estimate before pushing uses `git pack-objects --revs --thin --delta-base-offset --stdout` with the candidate commit as the positive revision and the last confirmed remote commit as the negative revision, piping its binary output directly to `wc -c`. Re-fetch the branch and reconcile concurrent non-Chandra changes before each fast-forward. If a single local-only unpushed commit exceeds the server limit, construct the smaller commits from its preserved tree in a fresh worktree based on the remote branch. Keep published history intact. [GitHub push-limit documentation](https://docs.github.com/en/get-started/using-git/troubleshooting-the-2-gb-push-limit)

After the last commit, run `python -m open_space.verify` and `python -m unittest discover -s tests` against the committed tree. A default-HEAD check before the evidence commit verifies only its parent. Also run the survey fixtures, catalog-store verification and a representative archived observation replay. Confirm the final remote commit and the single open draft PR before releasing temporary source files.
