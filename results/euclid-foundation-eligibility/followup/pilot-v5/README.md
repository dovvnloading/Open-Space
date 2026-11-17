# Five-arcminute development trial

The first five-arcminute trial produced 1,172 maxima, 269 selected centers, 54 positive paired-flux measurements and 43 tests with enough matched controls. Four raw crossings occur on visibly linear structures in the development images; their physical or instrumental origin is unassigned. The images are in `trigger-stamps.png`. The script conditioned on positive second-visit flux and was replaced by the signed-flux v6 method; these crossings do not provide a variable-source count or sensitivity.

The 677 fixed-grid control aperture centers and 269 source centers are deposited as CSV. Their corresponding 17-pixel science/RMS/flag stamps from both development visits are in `../pilot-v6/stamps.npz`. The frozen full-frame code is `../pilot_v5.py`, the revised signed-flux code is `../pilot_v6.py`; the latter reuses only first-visit positions. Full five-arcminute cutouts are retained on D: and have hashes in `cutout-receipts.json`.
