"""Catalog, checksum and number-tracing checks on small synthetic repositories."""

import gzip
import hashlib
import io
import json
import shutil
import subprocess
import tempfile
import unittest
import zipfile
from pathlib import Path

from open_space.verify import catalog, integrity, markdown, trace
from open_space.verify.evidence import Evidence
from open_space.verify.model import FAIL, INFO, PASS, WARN
from open_space.verify.tree import RepoTree


def sha(data):
    return hashlib.sha256(data).hexdigest()


def tree_of(files):
    return RepoTree({path: (data if isinstance(data, bytes) else data.encode("utf-8")) for path, data in files.items()})


def rules(findings, severity=None):
    return sorted(finding.rule for finding in findings if severity is None or finding.severity == severity)


class TreeTest(unittest.TestCase):
    def test_directory_listing_and_digest_are_reproducible(self):
        tree = tree_of({"results/b/x.csv": "1\n", "results/a/README.md": "# a\n", "results/ab/y": "2", "README.md": "r"})
        self.assertEqual(tree.result_dirs(), ["a", "ab", "b"])
        self.assertEqual(tree.under("results/a"), ["results/a/README.md"])
        self.assertTrue(tree.has_dir("results/ab"))
        self.assertFalse(tree.has_dir("results/c"))
        listing = "".join(f"{tree.digest(path)}  {path}\n" for path in tree.paths)
        self.assertEqual(tree.tree_digest(), sha(listing.encode()))

    @unittest.skipUnless(shutil.which("git"), "git is not installed")
    def test_git_view_returns_committed_bytes_not_checkout_bytes(self):
        with tempfile.TemporaryDirectory() as scratch:
            run = lambda *args: subprocess.run(["git", "-C", scratch, *args], check=True, capture_output=True)
            run("init", "-q")
            run("config", "user.email", "t@example.org")
            run("config", "user.name", "t")
            run("config", "core.autocrlf", "false")
            run("config", "commit.gpgsign", "false")
            Path(scratch, "data.csv").write_bytes(b"a,b\n1,2\n")
            run("add", "-A")
            run("commit", "-q", "-m", "first")
            Path(scratch, "data.csv").write_bytes(b"a,b\r\n1,2\r\n")      # what a converting checkout would hold
            Path(scratch, "new.txt").write_bytes(b"draft\n")
            committed = RepoTree.open(scratch)
            self.assertEqual(committed.read("data.csv"), b"a,b\n1,2\n")
            self.assertNotIn("new.txt", committed)
            self.assertEqual(committed.origin, "commit")
            working = RepoTree.open(scratch, worktree=True)
            self.assertIn("new.txt", working)
            self.assertEqual(working.read("data.csv"), b"a,b\r\n1,2\r\n")
            Path(scratch, "data.csv").write_bytes(b"a,b\n9,9\n")
            run("commit", "-q", "-am", "second")
            history = RepoTree.open(scratch).history("data.csv")
            self.assertEqual([digest for _, digest in history], [sha(b"a,b\n9,9\n"), sha(b"a,b\n1,2\n")])


class MarkdownTest(unittest.TestCase):
    def test_code_is_blanked_without_moving_anything(self):
        text = "see `python x.py 12`\n\n```sh\npython y.py 34\n```\nafter 56\n"
        prose, blocks = markdown.split_code(text)
        self.assertEqual(len(prose), len(text))
        self.assertNotIn("12", prose)
        self.assertNotIn("34", prose)
        self.assertIn("56", prose)
        self.assertEqual([(block.fenced, block.language, block.content.strip()) for block in blocks],
                         [(False, "", "python x.py 12"), (True, "sh", "python y.py 34")])

    def test_anchor_slugs_follow_github(self):
        self.assertEqual(markdown.slug("How to read the results"), "how-to-read-the-results")
        self.assertEqual(markdown.slug("Upcoming events you can watch for"), "upcoming-events-you-can-watch-for")
        offered = markdown.anchors("# Same\n\n## Same\n\n<a id=\"no-vanishing-supergiant\"></a>\n")
        self.assertEqual(offered, {"same", "same-1", "no-vanishing-supergiant"})


class CatalogTest(unittest.TestCase):
    README = ("# Project\n\nThe catalog below has 3 entries: 1 new measurements, 1 candidate discoveries, "
              "0 forecasts of specific events, 0 corrections to published catalogs and papers, and 1 carefully "
              "recorded searches that found nothing.\n\n| Label | Meaning |\n| --- | --- |\n| **Measured** | new |\n\n"
              "- **One.** Text. *Measured* · [Details](results/one/README.md)\n"
              "- **Two.** Text. *Candidate* · [Details](results/two/README.md#limits)\n"
              "- **Three.** Text. *Null result* · [Details](results/three/)\n")

    def files(self, **extra):
        files = {"README.md": self.README, "results/one/README.md": "# One\n", "results/two/README.md": "# Two\n\n## Limits\n",
                 "results/three/README.md": "# Three\n"}
        files.update(extra)
        return files

    def test_consistent_catalog_passes(self):
        findings = catalog.run(tree_of(self.files()))
        self.assertEqual(rules(findings, FAIL), [])
        self.assertEqual(rules(findings, WARN), [])
        self.assertIn("catalog/counts-agree", rules(findings, PASS))

    def test_legend_table_is_not_counted_as_an_entry(self):
        entries = catalog.catalog_entries(tree_of(self.files()))
        self.assertEqual([entry["label"] for entry in entries], ["Measured", "Candidate", "Null result"])

    def test_stated_counts_must_match_the_entries(self):
        wrong = self.README.replace("has 3 entries: 1 new", "has 4 entries: 2 new")
        findings = catalog.run(tree_of(self.files(**{"README.md": wrong})))
        (failure,) = [finding for finding in findings if finding.severity == FAIL]
        self.assertEqual(failure.rule, "catalog/count-mismatch")
        self.assertIn("states 4 entries, catalog has 3", failure.message)
        self.assertIn("states 2 'Measured', catalog has 1", failure.message)

    def test_broken_link_unlisted_folder_and_missing_anchor(self):
        files = self.files(**{"results/four/README.md": "# Four\n[gone](missing.csv) [up](../../README.md#nope)\n"})
        findings = catalog.run(tree_of(files))
        self.assertEqual(rules(findings, FAIL), ["catalog/link-missing"])
        self.assertEqual(rules(findings, WARN), ["catalog/anchor-missing", "catalog/result-unlisted"])
        self.assertEqual(catalog.labels_by_result(tree_of(files))["two"], ["Candidate"])


class IntegrityTest(unittest.TestCase):
    def check(self, files):
        findings, extra = integrity.run(tree_of(files))
        return findings, extra["census"]

    def test_manifest_lines_are_verified_and_mismatches_fail(self):
        good, bad = b"alpha\n", b"beta\n"
        manifest = f"{sha(good)}  data/good.csv\n{sha(b'other')}  data/bad.csv\n{sha(b'far away')}  inputs/big.fits\n"
        findings, census = self.check({"results/r/SHA256SUMS": manifest, "results/r/data/good.csv": good,
                                       "results/r/data/bad.csv": bad})
        self.assertEqual(census["path_bound_verified"], 1)
        self.assertEqual(census["listed_but_absent"], 1)
        self.assertEqual(rules(findings, FAIL), ["integrity/digest-mismatch"])
        self.assertEqual(rules(findings, INFO), ["integrity/listed-file-not-deposited"])

    def test_json_digest_key_decides_which_bytes_it_describes(self):
        data = b"x,y\n1,2\n"
        document = {
            "inputs": [{"file": "inputs/a.csv", "sha256": sha(data), "raw_sha256": sha(b"upstream response")}],
            "deposited": {"inputs/a.csv": sha(data)},
            "table": {"csv_file": "inputs/a.csv", "csv_sha256": sha(data), "source_sha256": sha(b"the fits file")},
            "note": {"upstream": "download.csv", "sha256": sha(b"a download"), "used_for": "inputs/a.csv"},
        }
        findings, census = self.check({"results/r/sources.json": json.dumps(document), "results/r/inputs/a.csv": data})
        self.assertEqual(rules(findings, FAIL), [])
        self.assertEqual(rules(findings, WARN), [])
        self.assertEqual(census["path_bound_verified"], 3)
        self.assertEqual(census["external"], 3)

    def test_nested_map_and_table_forms(self):
        data = b"payload"
        document = {"data/x.bin": {"url": "https://example.org/x", "sha256": sha(data), "raw_response_sha256": sha(b"r")}}
        table = f"path,sha256,source_sha256,url\ndata/x.bin,{sha(data)},{sha(b'up')},https://example.org/x\n"
        findings, census = self.check({"results/r/sources.json": json.dumps(document), "results/r/files.csv": table,
                                       "results/r/data/x.bin": data})
        self.assertEqual(rules(findings, FAIL), [])
        self.assertEqual(census["path_bound_verified"], 2)
        self.assertEqual(census["external"], 2)                  # the two upstream digests, each with its URL
        self.assertEqual(census["external_with_url"], 2)

    def test_mismatch_is_explained_when_it_is_a_line_ending_conversion(self):
        committed = b"a,b\n1,2\n"
        recorded = sha(committed.replace(b"\n", b"\r\n"))
        findings, _ = self.check({"results/r/SHA256SUMS": f"{recorded}  t.csv\n", "results/r/t.csv": committed})
        (warning,) = [finding for finding in findings if finding.severity == WARN]
        self.assertEqual(warning.rule, "integrity/digest-of-converted-bytes")
        self.assertIn("CRLF", warning.message)
        self.assertEqual(rules(findings, FAIL), [])

    def test_same_name_as_a_larger_upstream_file_is_a_warning_not_a_failure(self):
        document = [{"file": "gaia.csv", "url": "https://example.org/q", "bytes": 53558888, "sha256": sha(b"big")}]
        findings, _ = self.check({"results/r/sources.json": json.dumps(document), "results/r/gaia.csv": b"small extract"})
        self.assertEqual(rules(findings, WARN), ["integrity/name-shadows-upstream-file"])
        self.assertEqual(rules(findings, FAIL), [])

    def test_upstream_namesake_is_a_receipt_once_the_committed_file_has_its_own_digest(self):
        small = b"small extract"
        document = [{"file": "gaia.csv", "url": "https://example.org/q", "bytes": 53558888, "sha256": sha(b"big"),
                     "committed_table": {"file": "gaia.csv", "bytes": len(small), "sha256": sha(small)}}]
        findings, census = self.check({"results/r/sources.json": json.dumps(document), "results/r/gaia.csv": small})
        self.assertEqual(rules(findings, WARN) + rules(findings, FAIL), [])
        self.assertEqual(census["path_bound_verified"], 1)
        self.assertEqual(census["external"], 1)
        frozen = {"upstream": {"hosts.csv": {"url": "https://example.org/h", "bytes": 394503, "sha256": sha(b"big")}},
                  "frozen_files": {"hosts.csv": sha(small)}}
        findings, _ = self.check({"results/r/provenance.json": json.dumps(frozen), "results/r/hosts.csv": small})
        self.assertEqual(rules(findings, WARN) + rules(findings, FAIL), [])

    def test_digests_of_archive_members_and_gzip_payloads_count_as_committed(self):
        member, payload = b"light curve", b"decompressed table"
        bundle = io.BytesIO()
        with zipfile.ZipFile(bundle, "w") as archive:
            archive.writestr("lightcurves/29489.lc", member)
        manifest = f"{sha(member)}  lightcurves/29489.lc\n"
        document = {"raw_sha256": sha(payload)}
        findings, census = self.check({"results/r/sources.sha256": manifest, "results/r/inputs.zip": bundle.getvalue(),
                                       "results/r/t.csv.gz": gzip.compress(payload), "results/r/p.json": json.dumps(document)})
        self.assertEqual(census["content_verified"], 2)
        self.assertEqual(rules(findings, WARN), [])
        self.assertEqual(rules(findings, FAIL), [])

    @unittest.skipUnless(shutil.which("git"), "git is not installed")
    def test_stale_digest_is_traced_to_the_commit_it_matches(self):
        with tempfile.TemporaryDirectory() as scratch:
            run = lambda *args: subprocess.run(["git", "-C", scratch, *args], check=True, capture_output=True)
            run("init", "-q")
            run("config", "user.email", "t@example.org")
            run("config", "user.name", "t")
            run("config", "core.autocrlf", "false")
            run("config", "commit.gpgsign", "false")
            folder = Path(scratch, "results", "r")
            folder.mkdir(parents=True)
            (folder / "out.csv").write_bytes(b"v\n1\n")
            (folder / "SHA256SUMS").write_bytes(f"{sha(b'v' + bytes([10]) + b'1' + bytes([10]))}  out.csv\n".encode())
            run("add", "-A")
            run("commit", "-q", "-m", "first")
            (folder / "out.csv").write_bytes(b"v\n2\n")
            run("commit", "-q", "-am", "regenerate without updating the manifest")
            findings, _ = integrity.run(RepoTree.open(scratch))
            (failure,) = [finding for finding in findings if finding.severity == FAIL]
            self.assertEqual(failure.rule, "integrity/stale-digest")
            self.assertIn("matches_commit", failure.data)


class EvidenceTest(unittest.TestCase):
    def test_pools_separate_stored_values_from_computed_facts(self):
        rows = "".join(f"{index},{'g' if index % 3 else 'r'},{index * 0.5}\n" for index in range(1, 121))
        tree = tree_of({
            "results/r/README.md": "ignored 999\n",
            "results/r/summary.json": json.dumps({"period_d": 3.166, "events": [1, 2, 3], "note": "kept 17 of 20"}),
            "results/r/small.csv": "name,value\na,1.250\nb,2.5\n",
            "results/r/big.csv": "id,band,flux\n" + rows,
            "results/r/run.py": "THRESHOLD = 0.0025\nLIMIT = -7\nYEAR = 365.25 * 86400\n",
        })
        evidence = Evidence.build(tree, "results/r")
        self.assertIn(3.166, evidence.direct.values)
        self.assertIn(1.25, evidence.direct.values)
        self.assertNotIn(999.0, evidence.direct.values + evidence.derived.values + evidence.bulk.values)
        self.assertIn(120.0, evidence.derived.values)            # rows of big.csv
        self.assertIn(80.0, evidence.derived.values)             # rows with band == 'g'
        self.assertIn(17.0, evidence.derived.values)             # number inside a JSON string
        self.assertIn(60.0, evidence.bulk.values)                # a raw cell of the large table
        self.assertNotIn(60.0, evidence.direct.values)
        self.assertEqual(sorted(evidence.code.values), [-7.0, 0.0025, 365.25, 86400.0])
        index = evidence.direct.values.index(1.25)
        self.assertEqual(evidence.direct.places(index), 3)
        self.assertEqual(evidence.describe(evidence.direct, index), "small.csv [value] row 1")
        self.assertEqual(evidence.find_literal("0.0025"), "run.py")


class CodeConstantsTest(unittest.TestCase):
    def test_literals_are_read_without_parsing_the_script(self):
        from open_space.verify.evidence import code_numbers
        source = (
            "# threshold 99 in a comment\n"
            "THRESHOLD = 0.0025   # trailing 77\n"
            "LIMIT = -7\n"
            "x = a - 3\n"
            "y = (-.5, 1e-3, 0x1F, 1_000, 2j, f(-4), b[-1])\n"
            "s = 'text 123' + \"more 456\" + f'{value:.3f} 789'\n"
            'doc = """multi\nline 555 string"""\n'
            "name2 = v1 + abc3d\n"
            "result = 12 if z else -13\n"
            "match value:\n    case 5 -6: pass\n"          # newer syntax must not stop the reader
        )
        self.assertEqual(list(code_numbers(source)),
                         [(0.0025, 2), (-7.0, 3), (3.0, 4), (-0.5, 5), (0.001, 5), (31.0, 5), (1000.0, 5),
                          (-4.0, 5), (-1.0, 5), (12.0, 10), (-13.0, 10), (5.0, 12), (6.0, 12)])


class TraceTest(unittest.TestCase):
    def setUp(self):
        self.files = {
            "results/r/measurements.json": json.dumps(
                {"delay_s": 30.6996, "r": 0.968565, "rate": -8.19e-12, "fraction": 0.4251, "coordinate": 264.526816,
                 "count": 868, "depth": 7.3468}),
            "results/r/windows.csv": "id,flag\n" + "".join(f"{index},{'y' if index <= 843 else 'n'}\n"
                                                             for index in range(1, 945)),
        }

    def traces(self, text):
        tree = tree_of(dict(self.files, **{"results/r/README.md": text}))
        evidence = Evidence.build(tree, "results/r")
        _, traces, identifiers = trace.trace_document(text, evidence)
        return {item.number.text: item for item in traces}, identifiers

    def test_values_are_confirmed_at_their_printed_precision(self):
        found, _ = self.traces("The delay is 30.6996 s with r = 0.96857; the rate is 8.19×10⁻¹² and 42.51% remain.")
        self.assertEqual(found["30.6996"].status, "confirmed")
        self.assertEqual(found["30.6996"].where, "measurements.json #/delay_s")
        self.assertEqual(found["0.96857"].status, "confirmed")
        self.assertEqual(found["8.19×10⁻¹²"].how, "sign dropped in prose")
        self.assertEqual(found["42.51"].how, "percent as a fraction")

    def test_table_facts_such_as_counts_are_evidence(self):
        found, _ = self.traces("A stronger signal in 843 of 944 windows.")
        self.assertEqual(found["843"].where, "windows.csv [flag] (rows equal to 'y')")
        self.assertEqual(found["944"].where, "windows.csv (number of rows)")

    def test_a_wrong_number_is_unmatched(self):
        found, _ = self.traces("The delay is 12.3456 s and r = 0.95857.")
        self.assertEqual(found["12.3456"].status, "unmatched")
        self.assertEqual(found["0.95857"].status, "unmatched")

    def test_a_single_wrong_digit_of_a_stored_value_is_flagged(self):
        found, _ = self.traces("The delay is 30.7996 s.")           # stored: 30.6996
        self.assertEqual(found["30.7996"].status, "mismatch")
        self.assertEqual(found["30.7996"].note, "one digit differs")
        self.assertLessEqual(found["30.7996"].chance, trace.NEAR_ALPHA)

    def test_last_digit_rounding_is_reported_as_such(self):
        found, _ = self.traces("a depth of 7.346 here")            # 7.3468 truncated instead of rounded to 7.347
        self.assertEqual(found["7.346"].status, "confirmed")
        self.assertEqual(found["7.346"].note, trace.LOOSE_NOTE)

    def test_near_coincidence_with_a_stored_value_raises_a_mismatch(self):
        found, _ = self.traces("The expected position is 264.526819 degrees.")
        self.assertEqual(found["264.526819"].status, "mismatch")
        self.assertIn("3.0 units in the last printed place", found["264.526819"].note)

    def test_a_count_off_by_one_is_not_explained_away_as_rounding(self):
        found, _ = self.traces("We kept 867 sources.")
        self.assertNotEqual(found["867"].note, trace.LOOSE_NOTE)

    def test_crowded_evidence_weakens_a_match(self):
        self.files["results/r/dense.csv"] = "x,y\n" + "".join(f"{value / 100:.2f},1\n" for value in range(100, 140))
        found, _ = self.traces("One value is 1.23 and another is 30.6996.")
        self.assertEqual(found["1.23"].status, "found")
        self.assertGreater(found["1.23"].chance, trace.ALPHA)
        self.assertEqual(found["30.6996"].status, "confirmed")

    def test_occupancy_measures_how_crowded_the_neighbourhood_is(self):
        self.assertEqual(trace.occupancy([5.0], 5.0, 0.5), 0.0)
        everywhere = [float(value) for value in range(-100, 200)]
        self.assertEqual(trace.occupancy(everywhere, 50.0, 0.5), 1.0)
        half = [float(value) for value in range(0, 200, 2)]
        self.assertAlmostEqual(trace.occupancy(half, 50.0, 0.5), 0.5)

    def test_identifiers_must_appear_literally(self):
        self.files["results/r/targets.csv"] = "tic,x\n384244835,1\n"
        _, identifiers = self.traces("The source is TIC 384244835, not TIC 384244836.")
        self.assertEqual([(number.text, where) for number, where in identifiers],
                         [("384244835", "targets.csv"), ("384244836", None)])

    def test_run_reports_per_result_totals_and_checks_the_catalog(self):
        readme = "# R\n\nThe delay is 30.6996 s and the count is 12.345.\n"
        root = ("# Project\n\n- **R.** A delay of 30.6996 seconds. *Measured* · [Details](results/r/README.md)\n"
                "- **R again.** A delay of 31.6996 seconds. *Candidate* · [Details](results/r/README.md)\n")
        tree = tree_of(dict(self.files, **{"results/r/README.md": readme, "README.md": root}))
        findings, extra = trace.run(tree)
        totals = extra["per_result"]["r"]
        self.assertEqual((totals["precise"], totals["confirmed"], totals["unmatched"]), (2, 1, 1))
        self.assertEqual(extra["catalog"], {"entries": 2, "precise": 2, "confirmed": 1, "found": 0, "unmatched": 1})
        self.assertIn("trace/catalog-number-unmatched", rules(findings, WARN))
        self.assertIn("trace/unmatched", rules(findings, INFO))


if __name__ == "__main__":
    unittest.main()
