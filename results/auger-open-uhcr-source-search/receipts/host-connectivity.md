# Data-host connectivity receipt — 2026-11-19 03:20:51 UTC

Environment: local compute container (not GitHub Actions); Python `socket.getaddrinfo(host,443,type=SOCK_STREAM)` in concurrent threads. The probe is not a claim about the hosts' public uptime; this container has DNS disabled.

| Data host | DNS resolution from compute container | HTTPS data transfer |
| --- | --- | --- |
| opendata.auger.org | FAILED: `gaierror(-3, Temporary failure in name resolution)` | Blocked before connection |
| zenodo.org | FAILED: same error | Blocked before connection |
| zenodo-files.s3.amazonaws.com | FAILED: same error | Blocked before connection |
| raw.githubusercontent.com | FAILED: same error | Blocked before connection |
| api.github.com | FAILED: same error | Blocked before connection |
| cdsarc.cds.unistra.fr | FAILED: same error | Blocked before connection |
| fermi.gsfc.nasa.gov | FAILED: same error | Blocked before connection |

Probe invocation: Python `socket.getaddrinfo` on all hosts using port 443; outputs shown above. The first combined DNS/HTTPS probe hit a 60-second timeout, so the independent concurrent DNS-only rerun is the authoritative result. This prevents claiming raw event downloads or analyses occurred. Public metadata can separately be consulted via the web retrieval system, and repository writes via the GitHub connector; **neither is a direct download to the compute container**.
