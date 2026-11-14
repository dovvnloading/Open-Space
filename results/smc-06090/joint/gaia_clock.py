"""Gaia BJD(TCB) to HJD(UTC) for OGLE-SMC-ECL-06090. Requires pyerfa."""
import erfa
import numpy as np
RA, DEC = np.radians([21.06275, -73.32594444444444])
LINE_OF_SIGHT = np.array([np.cos(RA)*np.cos(DEC), np.sin(RA)*np.cos(DEC), np.sin(DEC)])

def to_hjd_utc(bjd_tcb):
    # ERFA epv00 gives position of Earth measured heliocentrically, barycentrically.
    # Their difference is the Sun's barycentric position. Subtract its projection.
    # Treat TDB as TT in the final clock conversion (sub-2-ms periodic error).
    a, b = erfa.tcbtdb(2451545., np.asarray(bjd_tcb)-2451545.)
    helio, bary = erfa.epv00(a, b)
    sun_projection = ((bary['p']-helio['p']) @ LINE_OF_SIGHT) * 499.004783836 / 86400
    tai1, tai2 = erfa.tttai(a,b)
    utc1, utc2 = erfa.taiutc(tai1,tai2)
    return utc1 + utc2 - sun_projection


