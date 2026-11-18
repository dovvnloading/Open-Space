import unittest

import numpy as np

from corridor import altitude_arc, angular_residual
from geometry import A, B, intersect, surface


class CorridorTests(unittest.TestCase):
    def test_all_nested_shell_directions_fall_on_endpoint_arc(self):
        station_height = 35786.0234375
        east = surface(-75.2, 0, A + station_height, B + station_height)
        west = surface(-137.0, 0, A + station_height, B + station_height)
        for lon, lat in [(-100, -20), (-150, 30), (-75, 55), (-105, 0)]:
            arc = altitude_arc(east, west, lon, lat)
            mapped_point = surface(lon, lat, A + 14, B + 6)
            for height in np.linspace(40, 120, 51):
                point = intersect(east, mapped_point, A + height, B + height)
                self.assertLess(angular_residual(arc, point - west), 2e-14)

    def test_off_arc_uses_projection_or_endpoint(self):
        arc = ([1, 0, 0], [0, 1, 0])
        self.assertAlmostEqual(angular_residual(arc, [1, 1, 1]), np.arctan(1 / np.sqrt(2)))
        self.assertAlmostEqual(angular_residual(arc, [1, -1, 0]), np.pi / 4)
        self.assertAlmostEqual(angular_residual(arc, [-1, -1, 0]), 3 * np.pi / 4)
        self.assertAlmostEqual(angular_residual(arc, [0, 0, 1]), np.pi / 2)

    def test_point_arc_and_invalid_stations(self):
        self.assertAlmostEqual(angular_residual(([1, 0, 0], [1, 0, 0]), [0, 1, 0]), np.pi / 2)
        with self.assertRaises(ValueError):
            angular_residual(([1, 0, 0], [-1, 0, 0]), [0, 1, 0])
        with self.assertRaises(ValueError):
            altitude_arc([0, 0, 0], [0, 42000, 0], 0, 0)


if __name__ == '__main__':
    unittest.main()
