# ACIS-S4 frame occupancy

**2 real instrument-passing S4 windows had at least three photons at one time tag.** The 24 scrambled streams had 0–0 such windows each, totaling 0.

The same count criterion was applied to all 25 streams across 4,344 fully scanned dated observations. Counts here are per event-file window, before cross-product deduplication. Multiple control realizations share the same images and event populations. Shuffling also disrupts real instrumental time correlations.

For the real windows, the photon tables additionally summarize how many photons shared both an identical time and a CHIPY row within a 256-column CHIPX node. 0 selected real windows lacked a deposited photon table. That row measurement is available only for the real events. An extended same-row pattern in the same exposure number can support an S4 readout explanation. A shared timestamp alone has several physical explanations. The comparisons were chosen after inspection of two events and do not change primary or held-out trial statistics.

[Selected windows and their exact inputs](flagged_windows.csv) · [Counts for all 25 streams](streams.csv) · [Contributing ledgers](ledgers.csv.gz) · [Measurement receipt](summary.json) · [CXC S4 readout description](https://cxc.cfa.harvard.edu/ciao/why/destreak.html)
