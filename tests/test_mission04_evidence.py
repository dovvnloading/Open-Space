"""Offline, standard-library regression of the frozen Mission 04 O4b evidence.

Runs in the existing repository CI without downloading strain or requiring NumPy.
Replays ledger totals and statistical identities, not the signal-processing step.
"""
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import unittest
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1] / "results" / "mission04-o4b-bursts"
DATA = ROOT / "results"
SHIFTS = (0,) + tuple(s for magnitude in range(4, 257, 4)
                      for s in (magnitude, -magnitude))
RATE_THRESHOLD = 20.0


def read_json(path):
    return json.loads(path.read_text())


def binomial_upper_tail(n, k, p):
    """P(X >= k), accumulated from an exact log-gamma binomial mass."""
    if not 0 < p < 1 or not 0 < k <= n:
        raise ValueError("Invalid binomial parameters")
    log_mass = (math.lgamma(n + 1) - math.lgamma(k + 1) -
                math.lgamma(n - k + 1) + k * math.log(p) +
                (n - k) * math.log1p(-p))
    mass = math.exp(log_mass)
    result = mass
    for i in range(k, n):
        mass *= (n - i) * p / ((i + 1) * (1 - p))
        result += mass
    return result


class Mission04EvidenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.campaign = read_json(DATA / "holdout" / "campaign.json")
        cls.selection_path = DATA / "holdout_selection.json"
        cls.selection = read_json(cls.selection_path)
        cls.audit = read_json(DATA / "audit_results.json")
        cls.pilot = read_json(DATA / "sensitivity.json")

    def test_frozen_manifest_is_disjoint(self):
        campaign, selection = self.campaign, self.selection
        self.assertEqual(
            hashlib.sha256(self.selection_path.read_bytes()).hexdigest(),
            campaign["selection_sha256"])
        blocks = selection["blocks"]
        self.assertEqual(len(blocks), campaign["n_frames"])
        self.assertEqual(len(blocks), 125)
        self.assertEqual(campaign["scheduled_s"], 4096 * len(blocks))
        starts = sorted(b["gps"] for b in blocks)
        self.assertEqual(len(starts), len(set(starts)))
        self.assertTrue(all(b - a >= 4096 for a, b in zip(starts, starts[1:])))
        exclusions = set(selection["pilot_exclusions_gps"])
        exclusions.add(selection["control_exclusion_gps"])
        self.assertTrue(exclusions.isdisjoint(starts))
        self.assertTrue(all(
            len(b["sources"][ifo]["md5"]) == 32
            for b in blocks for ifo in ("H1", "L1")))

    def test_full_monthly_ledgers_reproduce_campaign(self):
        campaign = self.campaign
        totals = [[0] * 5 for _ in SHIFTS]
        maxima = [0.] * len(SHIFTS)
        tail_counts = [0] * len(SHIFTS)
        shift_index = {shift: i for i, shift in enumerate(SHIFTS)}
        found, trials = set(), {}
        foreground_count, peak = 0, 0.

        for receipt in campaign["monthly_ledgers"]:
            path = DATA / "holdout" / receipt["file"]
            payload = path.read_bytes()
            self.assertEqual(len(payload), receipt["bytes"])
            self.assertEqual(hashlib.sha256(payload).hexdigest(),
                             receipt["sha256"])
            with ZipFile(io.BytesIO(payload)) as z:
                self.assertIsNone(z.testzip())
                names = set(z.namelist())
                gps_set = {int(name.split("/")[0]) for name in names}
                self.assertEqual(len(names), 4 * len(gps_set))
                for gps in sorted(gps_set):
                    self.assertNotIn(gps, found)
                    found.add(gps)
                    base = str(gps) + "/"
                    self.assertTrue(all(base + name in names for name in (
                        "epoch.json", "foreground.csv", "tail.csv",
                        "injections.json")))
                    epoch = json.loads(z.read(base + "epoch.json"))
                    self.assertEqual(epoch["gps"], gps)
                    self.assertEqual(epoch["configuration_sha256"],
                                     campaign["configuration_sha256"])
                    rows = epoch["slides"]
                    self.assertEqual(len(rows), len(SHIFTS))
                    for i, row in enumerate(rows):
                        self.assertEqual(row[0], SHIFTS[i])
                        self.assertEqual(row[2], 88 * row[1])
                        for j in range(5):
                            totals[i][j] += int(row[1 + j])
                        maxima[i] = max(maxima[i], row[-1])

                    baseline = {}
                    with z.open(base + "foreground.csv") as raw:
                        for row in csv.DictReader(io.TextIOWrapper(
                                raw, encoding="utf-8")):
                            sec = int(row["offset_second"])
                            self.assertNotIn(sec, baseline)
                            value = float(row["score"])
                            self.assertTrue(0 <= sec < 4096 and
                                            math.isfinite(value) and value >= 0)
                            baseline[sec] = value
                    self.assertEqual(len(baseline), epoch["live_s"])
                    self.assertEqual(max(baseline.values(), default=0.),
                                     rows[0][-1])
                    foreground_count += sum(
                        value >= RATE_THRESHOLD for value in baseline.values())
                    peak = max(peak, max(baseline.values(), default=0.))

                    with z.open(base + "tail.csv") as raw:
                        for row in csv.DictReader(io.TextIOWrapper(
                                raw, encoding="utf-8")):
                            self.assertGreaterEqual(float(row["score"]), 3.)
                            if float(row["score"]) >= RATE_THRESHOLD:
                                tail_counts[shift_index[int(row["shift_s"])]] += 1
                    for draw in json.loads(z.read(base + "injections.json")):
                        number = draw["trial"]
                        self.assertNotIn(number, trials)
                        trials[number] = draw
                        self.assertEqual(draw["gps_start"], gps)
                        self.assertEqual(draw["baseline"],
                                         baseline.get(draw["second"]))
                        self.assertEqual(draw["live"],
                                         draw["second"] in baseline)
                        if draw["live"]:
                            self.assertEqual(draw["recovered"],
                                             baseline[draw["second"]] < RATE_THRESHOLD
                                             <= draw["injected_score"])
                        else:
                            self.assertFalse(draw["recovered"])
                            self.assertIsNone(draw["injected_score"])

        self.assertEqual(found, {b["gps"] for b in self.selection["blocks"]})
        self.assertEqual(sorted(trials), list(range(8192)))
        self.assertEqual(sum(x["recovered"] for x in trials.values()), 6950)
        self.assertEqual(foreground_count, 0)
        for i, (totals_i, maximum) in enumerate(zip(totals, maxima)):
            self.assertEqual(campaign["campaign_slides"][i][1:6], totals_i)
            self.assertEqual(campaign["campaign_slides"][i][-1], maximum)
            self.assertEqual(tail_counts[i],
                             campaign["campaign_slides"][i][-2])
        self.assertEqual(campaign["HL_live_s"], totals[0][0])
        self.assertEqual(campaign["background_s"],
                         sum(row[0] for row in totals[1:]))
        self.assertEqual(campaign["foreground_loudest"], peak)
        self.assertEqual(
            campaign["background_campaigns_ge_loudest"],
            sum(m >= peak for m in maxima[1:]))
        self.assertEqual(
            campaign["campaign_rank_fraction"],
            (1 + sum(m >= peak for m in maxima[1:])) / len(SHIFTS))
        self.assertEqual(self.audit["holdout"]["recovered"], 6950)

    def test_confidence_limits_from_binomial_and_poisson(self):
        for rate, n, recovered, seconds in (
                (self.campaign["rate"], 8192, 6950, 512000),
                ({"pointwise90": {
                    "one_sided_alpha_per_term": .05,
                    "efficiency_lower": self.pilot["cases"]["primary"][
                        "lower_efficiency_95"],
                    "count_upper": self.pilot["cases"]["primary"][
                        "poisson95_count_upper"],
                    "effective_time_s_lower": self.pilot["cases"]["primary"][
                        "effective_live_s_lower"],
                    "rate_per_day_upper": self.pilot["cases"]["primary"][
                        "rate_upper_per_day_pointwise90"]
                }}, 2048, 1792, 40315)):
            for result in rate.values():
                if not isinstance(result, dict) or "efficiency_lower" not in result:
                    continue
                alpha = result["one_sided_alpha_per_term"]
                eps = result["efficiency_lower"]
                self.assertAlmostEqual(
                    binomial_upper_tail(n, recovered, eps),
                    alpha, delta=2e-8)
                self.assertAlmostEqual(
                    result["count_upper"], -math.log(alpha), delta=1e-11)
                self.assertAlmostEqual(
                    result["effective_time_s_lower"],
                    seconds * eps, delta=1e-6)
                self.assertAlmostEqual(
                    result["rate_per_day_upper"],
                    -math.log(alpha) * 86400 / (seconds * eps),
                    delta=1e-8)


if __name__ == "__main__":
    unittest.main()
