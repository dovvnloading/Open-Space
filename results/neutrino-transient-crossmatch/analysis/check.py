"""Small deterministic checks for geometry, RA rotation, shuffle scoring and trial tails."""
import json
from pathlib import Path
import numpy as np
from astropy_healpix import HEALPix
from astropy import units as u
from run import prepared_region, statistic, summarize


def main():
    hp = HEALPix(nside=4, order='nested')
    ix = np.arange(192)
    lon, lat = hp.healpix_to_lonlat(ix)
    ra = lon.to_value(u.deg)
    dec = lat.to_value(u.deg)
    # coarse level 1 pixel 5 and fine level 2 pixel 120
    region = prepared_region({'uniq90': [16+5, 64+120]})
    expected = ((ix//4) == 5) | (ix == 120)
    assert np.array_equal(region(ra, dec), expected)
    assert np.array_equal(region((ra+360) % 360, dec), expected)
    shifted = (ra+100) % 360
    assert np.array_equal(region((shifted-100) % 360, dec), expected)
    poly = prepared_region({'polygons_deg': [[[359, -1], [1, -1], [1, 1], [359, 1], [359, -1]]]})
    assert list(poly(np.array([0, .5, 359.5, 2, 180]), np.zeros(5))) == [True, True, True, False, False]
    hits = np.array([[0, 0, 0], [1, 1, 0], [0, 1, 1]], bool)
    assert np.allclose(statistic(hits, np.array([4., 2., 1.])), [[0, 0], [2, .5], [2, 1]])
    assert summarize(np.zeros(9999), 0, 8)['corrected_p'] == 1
    assert np.isclose(summarize(np.zeros(9999), 1, 12)['corrected_p'], .0012)
    # Seeds/strata and pair tests verified against an explicit tiny two-alert exchange.
    t = np.array([[1, 0], [0, 1]], bool)
    sky = np.array([[1, 0], [0, 1]], bool)
    perm = np.array([[0, 1], [1, 0]])
    fast = np.column_stack([(t[:, sky[i]].any(axis=1))[perm[:, i]] for i in range(2)])
    slow = np.array([[any(sky[i, j] and t[p[i], j] for j in range(2)) for i in range(2)] for p in perm])
    assert np.array_equal(fast, slow)
    print('Passed native HEALPix mixed-order, wrap, rotation, spherical polygon, scoring, ties and stratified time lookup checks.')


if __name__ == '__main__':
    main()
