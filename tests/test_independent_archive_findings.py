"""Independent standard-library checks for the #36 archival result packets."""
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1] / "results"
COMET = ROOT / "comet-407p-asteroid-75547"
M33 = ROOT / "m33n-1996-12a"
AU_KM = 149597870.7


def load_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def closest_linear(a, b):
    """Project the relative heliocentric state to its local minimum."""
    delta = [x - y for x, y in zip(a["xyz_au"], b["xyz_au"])]
    velocity = [x - y for x, y in zip(a["vxyz_au_per_day"], b["vxyz_au_per_day"])]
    dt = -sum(x * v for x, v in zip(delta, velocity)) / sum(v * v for v in velocity)
    distance = math.dist([x + dt * v for x, v in zip(delta, velocity)], [0., 0., 0.]) * AU_KM
    speed = math.sqrt(sum(v * v for v in velocity)) * AU_KM / 86400
    return dt, distance, speed


def great_circle_arcsec(a, b):
    """Haversine separation without the discovery script or SciPy."""
    d2r = math.pi / 180
    ra = (float(a["ra_deg"]) - float(b["ra_deg"])) * d2r
    dec1, dec2 = float(a["dec_deg"]) * d2r, float(b["dec_deg"]) * d2r
    h = math.sin((dec1 - dec2) / 2)**2 + math.cos(dec1) * math.cos(dec2) * math.sin(ra / 2)**2
    return 2 * math.asin(math.sqrt(max(0., min(1., h)))) * 180 / math.pi * 3600


class IndependentArchiveFindings(unittest.TestCase):
    def test_two_orbit_services_recover_close_passage(self):
        evidence = json.loads((COMET / "ephemerides.json").read_text())
        left, right = evidence["jpl"]["407P"], evidence["jpl"]["75547"]
        self.assertEqual([r["jd_tdb"] for r in left], [r["jd_tdb"] for r in right])
        bracket = [math.dist(a["xyz_au"], b["xyz_au"]) * AU_KM for a, b in zip(left, right)]
        i = min(range(len(bracket)), key=bracket.__getitem__)
        self.assertTrue(0 < i < len(bracket)-1, "minimum not bracketed")
        dt, km, speed = closest_linear(left[i], right[i])
        self.assertAlmostEqual(left[i]["jd_tdb"] + dt, 2461094.330271346, delta=0.0001)
        self.assertAlmostEqual(km, 31373, delta=2)
        self.assertAlmostEqual(speed, 8.05951, delta=0.0001)
        dt2, km2, speed2 = closest_linear(evidence["miriade"]["407P"], evidence["miriade"]["75547"])
        self.assertAlmostEqual(evidence["miriade"]["407P"]["jd_tt"] + dt2, 2461094.329249124, delta=0.0001)
        self.assertAlmostEqual(km2, 31688, delta=2)
        self.assertAlmostEqual(speed2, 8.05947, delta=0.0001)
        self.assertLess(abs(km2-km), 500)

    def test_historical_m33_match_and_independent_i_variability(self):
        positions = {r["id"]: r for r in load_rows(M33 / "original_positions.csv")}
        stars = load_rows(M33 / "mira_fields.csv")
        for historical, expected_id, expected_dist in [
            ("1996-1", "01342067+3034391", 0.3939),
            ("1997-1", "01342428+3042433", 1.5440),
        ]:
            star = min(stars, key=lambda r: great_circle_arcsec(positions[historical], r))
            self.assertEqual(star["id"], expected_id)
            self.assertEqual(star["classifier_votes"], "6")
            self.assertAlmostEqual(great_circle_arcsec(positions[historical], star), expected_dist, delta=0.001)
        rows = [r for r in load_rows(M33 / "ztf_dr24.csv")
                if r["case"] == "96" and r["filtercode"] == "zi" and r["catflags"] == "0"]
        mags = [float(r["mag"]) for r in rows]
        self.assertEqual(len(mags), 24)
        self.assertAlmostEqual(max(mags)-min(mags), 2.31987, delta=0.001)

    def test_comet_non_detection_replays_from_released_maps(self):
        completed = subprocess.run([sys.executable, str(COMET / "recovery.py")],
                                   capture_output=True, text=True, check=True, timeout=40)
        result = json.loads(completed.stdout)
        expected = json.loads((COMET / "recovery.json").read_text())
        self.assertEqual(result, expected)
        self.assertEqual(result["all"]["frames"], 32)
        self.assertLess(result["all"]["search_peak_snr"], result["all"]["control_peak_snr"])
        for split in result["all"]["chronological_halves"]:
            self.assertLess(split["other_half_snr"], 1.0)


if __name__ == "__main__":
    unittest.main()
