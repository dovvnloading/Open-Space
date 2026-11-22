# The public archive returned a CHEOPS light curve

The inventory contains 5,618 visit identifiers. At least one anonymous download yielded a FITS table with time and flux samples. The access run did not scan for transits.

## What the archive returned

| Request | HTTP status or recorded failure |
| --- | --- |
| DACE visit listing | 200 |
| ESA TAP | 200 |
| ESA inventory error | TAP returned XML or HTML: <p class="msg">Due to unexpected maintenance activities, <strong>ESASky is currently unavailable</strong>.</p>             <p class="msg">We are working to restore the service as soon as possib |
| DACE request lightcurves | 200 |
| DACE anonymous product lightcurves | 200 |
| DACE browse all product names | 200 |

The inventory records which service supplied each visit. The probes record exact queries, response status and SHA-256 digests of retrieved inputs. A successful metadata or portal request alone cannot establish photometry access. The recorded status was `ready`. Credential challenges, HTTP service failures and payload parsing failures appear separately in the receipts. The conclusion applies to these requests at the time of the run.

[Archive status](outputs/archive_status.json) · [Visit inventory](outputs/inventory.jsonl) · [Request and input receipts](outputs/archive_probes.jsonl) · [Checksums](outputs/SHA256SUMS)

## Reproduction

In the study folder, the commands used by the archive job were

```sh
python -m pip install -r requirements.txt
python archive_io.py inventory outputs
```

The archive is live and can change. Verify the deposited snapshot offline from the study folder with

```sh
(cd outputs && sha256sum -c SHA256SUMS)
```
