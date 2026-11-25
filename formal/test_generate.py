"""Unit tests for the explicitly trusted CSV-to-Lean adapter."""
import tempfile
import unittest
from pathlib import Path

import generate as g


class AdapterTests(unittest.TestCase):
    def test_exact_scaling(self):
        self.assertEqual(g.scaled("0.21", 100), 21)
        self.assertEqual(g.scaled("0.200", 1000), 200)
        self.assertEqual(g.scaled("139.3", 1000), 139300)
        self.assertEqual(g.scaled("0.2100000000000000000000000000000000", 100), 21)
        for value in ("0.211", "0.2100000000000000000000000000000001",
                      "-1", "1e-3", "NaN", "", " 1"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                g.scaled(value, 100)

    def test_count_missing_is_not_zero(self):
        self.assertIsNone(g.nat_or_none(""))
        self.assertEqual(g.nat_or_none("0"), 0)
        for value in ("-1", "1.0", "nan", " 1"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                g.nat_or_none(value)
        self.assertEqual(
            g.block_summary([(0, 0, True), (1, 3, True), (4, 2, False),
                             (None, 5, True), (7, None, False)]),
            [5, 3, 1, 1, 1, 5, 5, 2, 1, 1, 0, 1, 3],
        )

    def test_utc_labels(self):
        day = "2026-12-24"
        self.assertEqual(g.milliseconds(day + "T22:55:02.324Z", day), 82502324)
        self.assertEqual(g.milliseconds(day + "T00:00:00.000Z", day), 0)
        for value in ("2026-12-25T00:00:00.000Z", day + "T24:00:00.000Z",
                      day + "T23:59:60.000Z", day + "T22:55:02.32Z"):
            with self.subTest(value=value), self.assertRaises(ValueError):
                g.milliseconds(value, day)

    def test_malformed_tables(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "input.csv"
            for content in ("a,a\n1,2\n", "a,b\n1\n", "a,b\n1,2,3\n",
                            "a,b\n1,\"unterminated\n"):
                path.write_text(content, encoding="utf-8")
                with self.subTest(content=content), self.assertRaises((ValueError, g.csv.Error)):
                    g.read_rows(path, ["a", "b"])

    def test_duplicate_and_invalid_bursts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "windows.csv"
            header = "obsid,number,source,narrow_literal,narrow_plus69\n"
            good = "0000000001,1,other,0,1\n"
            path.write_text(header + good, encoding="utf-8")
            self.assertEqual(g.nicer_rows(path), [(0, 1, True)])
            for row in (good, "0000000001,01,other,2,3\n",
                        "0000000002,0,other,0,1\n",
                        "0000000002,1,,0,1\n",
                        "0000000002,1,SGR 1935+2154 ,0,1\n"):
                path.write_text(header + good + row, encoding="utf-8")
                with self.subTest(row=row), self.assertRaises(ValueError):
                    g.nicer_rows(path)


if __name__ == "__main__":
    unittest.main()
