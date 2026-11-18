"""Containment in the stored 90% highest-posterior-density NUNIQ cells.

ICRS coordinates are in degrees; RA wrapping and polar geometry use native
nested HEALPix. The same membership test works for scalar or vector arrays.
"""
import numpy as np
from astropy_healpix import HEALPix
from astropy import units as u


def prepare(uniq):
    """Represent the sparse MOC by merged half-open ranges at its finest order."""
    uniq = np.asarray(uniq, dtype=np.int64)
    orders = np.floor(np.log2(uniq / 4) / 2).astype(int)
    max_order = int(np.max(orders))
    cells = uniq - 4 * (4 ** orders)
    factor = 4 ** (max_order - orders)
    starts, ends = [], []
    for left, right in sorted(zip(cells * factor, (cells + 1) * factor)):
        if ends and left <= ends[-1]:
            ends[-1] = max(ends[-1], right)
        else:
            starts.append(left)
            ends.append(right)
    return HEALPix(nside=2 ** max_order, order="nested"), np.array(starts), np.array(ends)


def contains(prepared, ra_deg, dec_deg):
    hp, starts, ends = prepared
    cells = hp.lonlat_to_healpix(np.asarray(ra_deg) * u.deg, np.asarray(dec_deg) * u.deg)
    i = np.searchsorted(starts, cells, side="right") - 1
    return (i >= 0) & (cells < ends[np.maximum(0, i)])
