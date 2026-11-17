"""Synthetic regression checks. No test fixture is a DESI observation."""
from dataclasses import replace
import unittest

import numpy as np

from laser_search import (Exposure, convolve_template, inject_line, isotropic_power_w,
                          scan_exposure, scan_target)
from recovery import measure_completeness, summarize_recovery


def fixture(expid=1, velocity=0, offset=0):
    w = 5000 + 0.8*np.arange(160) + offset
    x = np.arange(-5, 6)
    # Asymmetric line shape tests the sense of the DIA offset indexing.
    p = np.exp(-0.5*(x/1.18)**2) * (1 + 0.08*x)
    r = np.repeat((p/p.sum())[:, None], len(w), axis=1)
    return Exposure(expid, w, np.ones(w.size), np.ones(w.size)*4, np.zeros(w.size, int),
                    r, np.ones(w.size), velocity)


class LaserTests(unittest.TestCase):
    def test_asymmetric_column_orientation(self):
        e = fixture()
        e.resolution[3, 80] *= 2
        index = 80
        delta = np.zeros(len(e.wave))
        delta[index] = 1
        measured = convolve_template(e.resolution, delta)
        np.testing.assert_allclose(measured[index-5:index+6], e.resolution[:, index])
        injected = inject_line(e, e.wave[index], 45)
        np.testing.assert_allclose(measured / np.sum(measured * np.gradient(e.wave)),
                                   (injected.flux - e.flux) / 45, atol=1e-14)
        t = inject_line(e, e.wave[index]+0.31, 45)
        self.assertAlmostEqual((t.flux - e.flux) @ np.gradient(e.wave), 45)

    def test_persistent_offset_grids_rv_and_local_affine(self):
        e1, e2 = fixture(), fixture(2, 120, 0.3)
        truth = 5058.3436
        injected = [inject_line(e, truth, f) for e, f in [(e1, 35.), (e2, 25.)]]
        injected[0] = replace(injected[0], flux=injected[0].flux + 3 + 0.2*(e1.wave-truth))
        result = scan_target(injected, [1, 2])
        self.assertEqual(len(result.candidates), 1)
        c = result.candidates[0]
        self.assertAlmostEqual(c.rest_A, truth, places=6)
        self.assertAlmostEqual(c.hits[0].flux17, 35., places=6)
        self.assertAlmostEqual(c.hits[1].flux17, 25., places=6)
        self.assertEqual(len(scan_target([injected[0], e2], [1,2]).candidates), 0)
        drifting = inject_line(e2, truth+1.8, 25.)
        self.assertEqual(len(scan_target([injected[0], drifting], [1,2]).candidates), 0)

    def test_missing_duplicate_masked_and_edge(self):
        a, b = fixture(), fixture(3)
        with self.assertRaises(ValueError):
            scan_target([a], [1, 3])
        with self.assertRaises(ValueError):
            scan_target([a, a], [1, 3])
        with self.assertRaises(ValueError):
            scan_target([a], [1])
        b.mask[80] = 1
        truth = a.wave[80]
        result = scan_target([inject_line(e, truth, 40) for e in (a, b)], [1, 3])
        self.assertEqual(len(result.candidates), 0)
        self.assertFalse(any(lo <= truth <= hi for lo, hi in result.coverage))
        self.assertFalse(any(lo <= a.wave[0] <= hi for lo, hi in result.coverage))

    def test_shape_alternatives(self):
        e = fixture()
        cosmic = e.flux.copy()
        cosmic[77] += 60
        hits, _ = scan_exposure(replace(e, flux=cosmic))
        peak = max(hits, key=lambda h: h.snr)
        self.assertGreater(peak.spike_delta_chi2, 9)
        self.assertFalse(peak.shape_pass)
        broad = e.flux + 30*np.exp(-0.5*((e.wave-e.wave[76])/2.5)**2)
        peak = max(scan_exposure(replace(e, flux=broad))[0], key=lambda h: h.snr)
        self.assertGreater(peak.broad_delta_chi2, 9)
        self.assertFalse(peak.shape_pass)

    def test_recovery_veto_and_distance(self):
        exposures = [fixture(), fixture(3)]
        refits = []
        def refit(records):
            refits.append(True)
            return records
        def veto(records, candidate):
            return "site_control" if candidate.rest_A < 5050 else None
        records = measure_completeness(exposures, [1, 3], (4990, 5100), [0.1, 25], 20,
                                       17, refit=refit, veto=veto)
        self.assertEqual(len(records), 40)
        self.assertEqual(len(refits), sum(t["eligible"] for t in records))
        self.assertEqual({r["reason"] for r in records},
                         {"recovered", "missed_detection", "veto_site_control", "outside_joint_coverage"})
        for r in records:
            self.assertEqual(r["recovered"], r["flux17"] == 25 and r["rest_A"] >= 5050)
        summary, limit = summarize_recovery(records)
        self.assertIsNone(limit)
        self.assertEqual(summary[0]["trials"], 20)
        self.assertAlmostEqual(isotropic_power_w(1., 1.)/1.19649518e14, 1., places=6)
        self.assertAlmostEqual(isotropic_power_w(1., 2.)/isotropic_power_w(1., 1.), 4.)

    def test_no_threshold_without_evidence_at_all_brighter_levels(self):
        ledger = [dict(index=i, seed=0, eligible=True, flux17=f, recovered=(f>=2))
                  for i in range(100) for f in [1, 2, 3]]
        _, limit = summarize_recovery(ledger)
        self.assertEqual(limit, 2)
        for row in ledger[-18:]:
            if row["flux17"] == 3:
                row["recovered"] = False
        _, limit = summarize_recovery(ledger)
        self.assertIsNone(limit)


if __name__ == "__main__":
    unittest.main()
