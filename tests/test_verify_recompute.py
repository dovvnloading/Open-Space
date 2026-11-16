"""Independent estimators, the safe expression language and the claim engine."""

import json
import math
import unittest

from open_space.verify import expr, recompute, runner, stats
from open_space.verify.model import FAIL, INFO, PASS, WARN
from open_space.verify.tree import RepoTree


def tree_of(files):
    return RepoTree({path: (data if isinstance(data, bytes) else data.encode("utf-8")) for path, data in files.items()})


class StatisticsTest(unittest.TestCase):
    """Each estimator against a value that can be checked by hand or in any reference table."""

    def test_correlation_and_line(self):
        x, y = [1, 2, 3, 4, 5], [2.0, 4.1, 5.9, 8.2, 9.8]
        # By hand: Sxy = 19.7, Sxx = 10, Syy = 38.9, so r = 19.7 / sqrt(389).
        self.assertAlmostEqual(stats.pearson(x, y), 19.7 / math.sqrt(389.0), places=12)
        slope, intercept = stats.linear_fit(x, y)
        self.assertAlmostEqual(slope, 1.97, places=9)
        self.assertAlmostEqual(intercept, 0.09, places=9)
        self.assertAlmostEqual(stats.pearson(x, [-value for value in y]), -19.7 / math.sqrt(389.0), places=12)
        self.assertEqual(stats.spearman([1, 2, 3, 4], [10, 30, 20, 40]), 0.8)
        self.assertEqual(stats.ranks([10, 20, 20, 30]), [1.0, 2.5, 2.5, 4.0])
        with self.assertRaises(ValueError):
            stats.pearson([1, 1, 1], [1, 2, 3])

    def test_descriptive_statistics_skip_missing_values(self):
        data = [2.0, 4.0, math.nan, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0]
        self.assertEqual(stats.mean(data), 5.0)
        self.assertEqual(stats.median(data), 4.5)
        self.assertEqual(stats.std(data, 0), 2.0)
        self.assertAlmostEqual(stats.std(data), math.sqrt(32 / 7))
        self.assertEqual(stats.quantile([1, 2, 3, 4, 5], 0.25), 2.0)
        self.assertEqual(stats.mad([1, 1, 2, 2, 4, 6, 9]), 1.0)

    def test_tail_probabilities(self):
        self.assertAlmostEqual(stats.norm_sf(1.959963985), 0.025, places=9)
        self.assertAlmostEqual(stats.norm_sf(5.0), 2.866515719e-7, delta=1e-15)
        self.assertAlmostEqual(stats.norm_isf(0.025), 1.959963985, places=8)
        self.assertAlmostEqual(stats.chi2_sf(3.841458821, 1), 0.05, places=9)
        self.assertAlmostEqual(stats.chi2_sf(18.30703805, 10), 0.05, places=9)
        self.assertAlmostEqual(stats.chi2_sf(6.6, 1), 0.010195, places=5)
        self.assertAlmostEqual(stats.binom_sf(8, 10, 0.5), 56 / 1024)
        self.assertAlmostEqual(stats.binom_cdf(0, 300, 0.01), 0.99 ** 300)
        self.assertAlmostEqual(stats.poisson_sf(3, 1.0), 1 - math.exp(-1) * 2.5)
        self.assertAlmostEqual(stats.poisson_cdf(2, 1.0), math.exp(-1) * 2.5)
        self.assertAlmostEqual(sum(stats.poisson_pmf(k, 4.2) for k in range(60)), 1.0)

    def test_periodogram_recovers_an_injected_period(self):
        period = 3.83 / 24.0                                        # days
        times = [0.0137 * index + 0.004 * math.sin(index) for index in range(400)]
        values = [10.0 + 0.3 * math.sin(2 * math.pi * time / period + 0.4) + 0.02 * math.cos(17.0 * index)
                  for index, time in enumerate(times)]
        found = stats.lomb_scargle_peak(times, values, 0.12, 0.2, samples=600)
        self.assertAlmostEqual(found, period, delta=2e-4)
        self.assertGreater(stats.lomb_scargle_power(times, values, period),
                           20 * stats.lomb_scargle_power(times, values, 0.137))

    def test_sky_and_time_scales(self):
        self.assertAlmostEqual(stats.angular_separation(10, 0, 10, 1), 1.0)
        self.assertAlmostEqual(stats.angular_separation(0, 90, 123, 0), 90.0)
        self.assertAlmostEqual(stats.angular_separation(0, 0, 180, 0), 180.0)
        by_cosine_rule = math.degrees(math.acos(0.75 + 0.25 * math.cos(math.radians(1.0))))
        self.assertAlmostEqual(stats.angular_separation(10, 60, 11, 60), by_cosine_rule, places=9)
        self.assertEqual(stats.tai_minus_utc("2016-12-31"), 36)
        self.assertEqual(stats.tai_minus_utc("2017-01-01"), 37)
        self.assertEqual(stats.tt_minus_utc("2018-11-22T11:58:14"), 69.184)
        self.assertEqual(stats.tt_minus_utc("1999-06-01"), 64.184)
        self.assertEqual(stats.julian_date("2000-01-01T12:00:00"), 2451545.0)
        self.assertEqual(stats.modified_julian_date("2021-12-11"), 59559.0)


class ExpressionTest(unittest.TestCase):
    def test_arithmetic_constants_and_identities(self):
        self.assertAlmostEqual(expr.evaluate("(2*pi)**(5/3)"), 21.39438370764715)
        self.assertEqual(expr.evaluate("comb(18, 3)"), 816)
        self.assertAlmostEqual(expr.evaluate("1 - 0.05**(1/300)"), 0.00993608194445772)
        self.assertEqual(expr.evaluate("tt_minus_utc('2020-05-01')"), 69.184)
        self.assertEqual(expr.evaluate("3 if 2 > 1 else 4"), 3)
        self.assertTrue(expr.evaluate("1 < 2 <= 2 and not 3 == 4"))

    def test_columns_masks_and_comprehensions(self):
        self.assertEqual(expr.evaluate("sum((vec(1, 2, 3, 4) * 2)[vec(1, 2, 3, 4) > 2])"), 14.0)
        self.assertEqual(expr.evaluate("count((vec(1, 2, 3, 4) > 1) & ~(vec(1, 2, 3, 4) == 3))"), 2)
        self.assertEqual(expr.evaluate("[x * x for x in range(5) if x != 2]"), [0, 1, 9, 16])
        self.assertEqual(expr.evaluate("mean(vec(1, nan, 3))"), 2.0)
        self.assertEqual(expr.evaluate("dot([1, 2, 3], unit(0, 90))"), 3.0)

    def test_tables_and_json_come_only_from_the_claim_folder(self):
        tree = tree_of({
            "results/r/t.csv": "band,flux,flag\ng,1.5,0\nr,2.5,0\ng,,1\nr,4.0,0\n",
            "results/r/m.json": json.dumps({"a": {"b": [10, 20, 30]}}),
            "results/other/secret.csv": "x,y\n1,2\n",
        })
        evaluator = expr.Evaluator(tree.read, "results/r")
        evaluator.define("t", "table('t.csv')")
        self.assertEqual(evaluator.evaluate("len(t)"), 4)
        self.assertEqual(evaluator.evaluate("count(t['band'] == 'r')"), 2)
        self.assertEqual(evaluator.evaluate("sum(t[t['flag'] == 0]['flux'])"), 8.0)
        self.assertEqual(evaluator.evaluate("count(~isnan(t['flux']))"), 3)
        self.assertEqual(evaluator.evaluate("json('m.json', '/a/b/1')"), 20)
        self.assertEqual(evaluator.evaluate("json('m.json')['a']['b'][2]"), 30)
        for escape in ("table('../other/secret.csv')", "table('/etc/passwd')", "table('C:/x.csv')", "table('nope.csv')"):
            with self.assertRaises(expr.ExprError):
                evaluator.evaluate(escape)

    def test_long_identifiers_stay_exact_and_text_can_be_tested(self):
        tree = tree_of({"results/r/ids.csv": "system,published_id,matched_id\n"
                                             "WDWDMS a,5041268200697001600,5041268200697001600\n"
                                             "WDWDWD b,5041268200697001600,5041268200697001601\n"
                                             "WDWDWD c,4190500054845023488,4190500054845023488\n"})
        evaluator = expr.Evaluator(tree.read, "results/r")
        evaluator.define("t", "table('ids.csv')")
        # These two identifiers are the same number once rounded to a float; they must still differ here.
        self.assertEqual(float(5041268200697001600), float(5041268200697001601))
        self.assertEqual(evaluator.evaluate("count(t['published_id'] != t['matched_id'])"), 1)
        self.assertEqual(evaluator.evaluate("nunique(t['matched_id'])"), 3)
        self.assertEqual(evaluator.evaluate("count(startswith(t['system'], 'WDWDWD'))"), 2)
        self.assertEqual(evaluator.evaluate("count(contains(t['system'], ' a'))"), 1)
        self.assertEqual(evaluator.evaluate("count(endswith(t['system'], 'c'))"), 1)

    def test_nothing_outside_the_language_can_run(self):
        for source in ("__import__('os')", "(1).__class__", "open('x')", "'a' * 10**9", "9**9**9",
                       "[x for x in range(10**12)]", "lambda: 1", "exec('1')", "f'{1}'", "(x := 3)",
                       "factorial(10**6)", "[1, 2] * 5", "vec(1, 2).data"):
            with self.assertRaises(expr.ExprError, msg=source):
                expr.evaluate(source)

    def test_runaway_expressions_stop(self):
        evaluator = expr.Evaluator()
        with self.assertRaises(expr.ExprError):
            evaluator.evaluate("sum([1 for a in range(3000) for b in range(3000) for c in range(3000)])")


class ClaimEngineTest(unittest.TestCase):
    README = "# Result\n\nPearson r = 0.9988 in 5 samples. The window count is **843 /** of 944.\n"
    TABLE = "x,y,win\n1,2.0,1\n2,4.1,1\n3,5.9,0\n4,8.2,1\n5,9.8,1\n"

    def run_claims(self, claims, readme=None, extra=None):
        document = {"schema": recompute.SCHEMA, "let": {"t": "table('t.csv')"}, "claims": claims}
        if extra:
            document.update(extra)
        tree = tree_of({"results/r/README.md": readme or self.README, "results/r/t.csv": self.TABLE,
                        "results/r/claims.json": json.dumps(document)})
        findings, summary = recompute.run(tree)
        return findings, summary, tree

    def test_a_claim_holds_when_the_recomputed_value_rounds_to_the_printed_one(self):
        findings, summary, _ = self.run_claims([
            {"id": "r", "printed": "0.9988", "compute": "pearson(t['x'], t['y'])"},
            {"id": "n", "printed": "5 samples", "compute": "len(t)"},
            {"id": "wins", "expect": 4, "compute": "count(t['win'] == 1)"},
            {"id": "bound", "expect": True, "compute": "slope(t['x'], t['y']) > 1.9"},
        ])
        self.assertEqual([finding.severity for finding in findings], [PASS] * 4)
        self.assertEqual(summary, {"claim_files": 1, "claims": 4, "reproduced": 4, "failed": 0})
        self.assertAlmostEqual(findings[0].data["computed"], 19.7 / math.sqrt(389.0), places=12)
        self.assertEqual(findings[0].data["line"], 3)

    def test_a_number_that_disagrees_fails(self):
        findings, _, _ = self.run_claims([{"id": "r", "printed": "0.9988", "compute": "pearson(t['x'], t['win'])"}])
        self.assertEqual((findings[0].severity, findings[0].rule), (FAIL, "recompute/disagrees"))
        self.assertIn("independent recomputation gives", findings[0].message)

    def test_editing_the_readme_makes_the_claim_stale(self):
        findings, _, _ = self.run_claims([{"id": "r", "printed": "0.9988", "compute": "pearson(t['x'], t['y'])"}],
                                         readme="# Result\n\nPearson r = 0.9987 in 5 samples.\n")
        self.assertEqual(findings[0].rule, "recompute/claim-not-in-document")

    def test_tolerances_and_exact_expectations(self):
        findings, _, _ = self.run_claims([
            {"id": "loose", "printed": "0.9988", "compute": "0.9991", "tolerance": {"ulps": 3}},
            {"id": "tight", "printed": "0.9988", "compute": "0.9991", "tolerance": {"abs": 1e-5}},
            {"id": "float", "expect": 0.1, "compute": "0.1 + 1e-15"},
            {"id": "true-is-not-one", "expect": 1, "compute": "2 > 1"},
        ])
        self.assertEqual([finding.severity for finding in findings], [PASS, FAIL, PASS, FAIL])

    def test_broken_specifications_are_failures_not_crashes(self):
        findings, _, _ = self.run_claims([
            {"id": "a", "printed": "0.9988"},
            {"id": "b", "printed": "0.9988", "expect": 1, "compute": "1"},
            {"id": "c", "printed": "843 /** of 944", "compute": "1"},
            {"id": "d", "printed": "0.9988", "compute": "table('missing.csv')"},
            {"id": "d", "printed": "0.9988", "compute": "1"},
            {"id": "e", "printed": "0.9988", "compute": "__import__('os').system('x')"},
        ])
        self.assertEqual([finding.rule for finding in findings],
                         ["recompute/spec-invalid", "recompute/spec-invalid", "recompute/spec-invalid",
                          "recompute/error", "recompute/spec-invalid", "recompute/error"])
        tree = tree_of({"results/r/claims.json": "{not json"})
        (finding,) = recompute.check_file(tree, "results/r/claims.json")
        self.assertEqual(finding.rule, "recompute/spec-invalid")

    def test_acknowledged_warnings_are_downgraded_with_their_reason(self):
        readme = "# R\n\nThe catalogue period is 0.3604030 d.\n"
        document = {"schema": recompute.SCHEMA, "claims": [],
                    "acknowledged": [{"rule": "trace/mismatch", "key": "0.3604030",
                                      "reason": "catalogue value, quoted for comparison"}]}
        tree = tree_of({"results/r/README.md": readme, "results/r/m.json": json.dumps({"period": 0.3604025}),
                        "results/r/claims.json": json.dumps(document)})
        ledger = runner.run_checks(tree, ("trace", "recompute"))
        (finding,) = [item for item in ledger["findings"] if item.rule == "trace/mismatch"]
        self.assertEqual(finding.severity, INFO)
        self.assertEqual(finding.data["acknowledged"], "catalogue value, quoted for comparison")
        self.assertEqual(ledger["summary"]["findings"][WARN], 0)
        del document["acknowledged"]
        tree = tree_of({"results/r/README.md": readme, "results/r/m.json": json.dumps({"period": 0.3604025}),
                        "results/r/claims.json": json.dumps(document)})
        ledger = runner.run_checks(tree, ("trace", "recompute"))
        self.assertEqual(ledger["summary"]["findings"][WARN], 1)


if __name__ == "__main__":
    unittest.main()
