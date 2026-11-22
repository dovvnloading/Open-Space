# Archived catalog responses

Every HTTP memo remains at its original path under `evidence/catalogs/`. The loose binary files include replies needed by instrument-passing real triggers and directly linked review evidence. Each remaining reply retains its original logical pathname, byte length and SHA-256 in its producer’s ZIP member manifest. Each producer index records its sorted logical file inventory checksum and the compressed hash of each ZIP. Open a producer directory and download a ZIP from GitHub to read its `MEMBERS.json`. Large logical replies are split into ordered members. The acquisition memos include failed HTTP responses and retry histories.

From the repository root in the Chandra workflow, locate one reply, verify all bundles, or reconstruct a cache in a separate directory.

```bash
python results/chandra_acis_fxt/catalog_pack.py where catalogs/CACHE/REPLY.bin
python results/chandra_acis_fxt/catalog_pack.py verify
python results/chandra_acis_fxt/catalog_pack.py extract --cache CACHE --destination "$RUNNER_TEMP/catalog-view"
python results/chandra_acis_fxt/catalog_pack.py replay OBSID
```

The replay command restores the caches used by that observation in a temporary directory and invokes the original checksummed observation replay. The logical filenames and response bytes supplied to the frozen classifier are identical to the original working files. [Operational details](../../CATALOG_STORAGE.md)
