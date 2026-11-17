# JWST repeat-imaging transient search

No new extragalactic transient has been established at this checkpoint. The search has acquired actual single-epoch NEXUS science and error maps. The first full-footprint subtraction compares NIRCam F444W from 2026 May 28 and August 7. A coverage-selected central square also has F200W in both epochs. The newer September epoch has public MAST products and is being acquired separately.

[qdr-inputs.csv](qdr-inputs.csv) records the URLs and SHA-256 hashes of the first eight compressed archive inputs and the derived central cutouts. [screen.py](screen.py) fits a relative astrometric translation and empirical compact-source PSFs, cross-convolves the epochs, searches both signs of the difference and measures conditional point-source recovery by injection on the observed background. Calibration, measured area, candidate vetting and depths are pending.

Data credit: the NEXUS collaboration, [QDR archive and release documentation](https://nexus-jwst.github.io/updates/). [Survey overview (Shen et al.)](https://arxiv.org/abs/2408.12777) · [Release notes (Zhuang et al.)](https://arxiv.org/abs/2603.04586) · [Published year-one transients](https://arxiv.org/abs/2609.06985).
