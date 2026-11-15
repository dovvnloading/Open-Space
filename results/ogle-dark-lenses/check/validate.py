"""Check pspl.py's annual-parallax magnification against MulensModel.

Usage: python -I check/validate.py
"""
import os
import sys

import numpy as np
import MulensModel as mm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pspl as lens  # noqa: E402

t = np.linspace(2460700.0, 2461500.0, 400)
dn, de = lens.parallax_offsets(t)
for u0, pn, pe in [(0.5, 0.05, 0.02), (-0.5, 0.05, 0.02), (0.5, -0.15, 0.1), (-0.3, 0.2, -0.05)]:
    mine = lens.magnification([2461213.0, u0, 275.0, pn, pe], t, dn, de)
    model = mm.Model({"t_0": 2461213.0, "u_0": u0, "t_E": 275.0, "pi_E_N": pn, "pi_E_E": pe,
                      "t_0_par": lens.TPAR}, coords=lens.RA_DEC)
    ref = model.get_magnification(t)
    flipped = lens.magnification([2461213.0, -u0, 275.0, pn, pe], t, dn, de)
    print("u0=%+.2f piEN=%+.2f piEE=%+.2f  max|dA| same-sign=%.2e  flipped-sign=%.2e"
          % (u0, pn, pe, np.abs(mine - ref).max(), np.abs(flipped - ref).max()))
