"""Replay of documented commands in an offline sandbox, the ledger and the command line."""

import contextlib
import io
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from open_space.verify import __main__ as cli
from open_space.verify import replay, report, runner
from open_space.verify.model import FAIL, INFO, PASS, WARN
from open_space.verify.tree import RepoTree

README = """# Synthetic result

Run `python make.py` to regenerate the table and `python check.py --check` to verify it.

```sh
python drift.py
python fetch.py --cache /tmp/cache
python needs_package.py
python gone.py
python broken.py
python slow.py
python tool.py <observation-id>
python summary.py > summary.json
python quiet.py
```
"""

FILES = {
    "results/demo/README.md": README,
    "results/demo/table.csv": "name,value\na,1.25\nb,2.5\n",
    "results/demo/drifting.csv": "name,value\na,1.25\n",
    "results/demo/summary.json": '{"n": 2, "mean": 1.875}\n',
    "results/demo/make.py": "from pathlib import Path\nPath('table.csv').write_text('name,value\\na,1.25\\nb,2.5\\n')\n",
    "results/demo/check.py": "import sys\nfrom pathlib import Path\nassert '--check' in sys.argv\n"
                             "assert Path('table.csv').read_text().splitlines()[1] == 'a,1.25'\n",
    "results/demo/drift.py": "from pathlib import Path\nPath('drifting.csv').write_text('name,value\\na,1.26\\n')\n",
    "results/demo/fetch.py": "import urllib.request\nurllib.request.urlopen('http://example.org/data.csv', timeout=5)\n",
    "results/demo/needs_package.py": "import package_that_does_not_exist_anywhere\n",
    "results/demo/broken.py": "raise RuntimeError('the fit did not converge')\n",
    "results/demo/slow.py": "import time\ntime.sleep(60)\n",
    "results/demo/tool.py": "print('never run')\n",
    "results/demo/summary.py": "import json\nprint(json.dumps({'n': 2, 'mean': 1.875}))\n",
    "results/demo/quiet.py": "x = 1\n",
}


def tree_of(files):
    return RepoTree({path: data.encode("utf-8") for path, data in files.items()})


class DiscoveryTest(unittest.TestCase):
    def test_commands_are_read_from_code_spans_and_shell_blocks(self):
        commands = replay.discover(tree_of(FILES))
        self.assertEqual([" ".join(command.argv) for command in commands],
                         ["make.py", "check.py --check", "drift.py", "fetch.py --cache /tmp/cache", "needs_package.py",
                          "gone.py", "broken.py", "slow.py", "tool.py <observation-id>", "summary.py", "quiet.py"])
        by_script = {command.argv[0]: command for command in commands}
        self.assertEqual(by_script["make.py"].script, "results/demo/make.py")
        self.assertEqual(by_script["make.py"].cwd, "results/demo")
        self.assertIsNone(by_script["gone.py"].script)
        self.assertEqual(by_script["summary.py"].redirect, "summary.json")

    def test_cd_inside_a_shell_block_moves_the_following_commands(self):
        text = "~~~sh\ncd screen\npython fetch.py\ncd ..\npython top.py\n~~~\n\n`cd tools && python build.py`\n"
        files = {"results/demo/README.md": text, "results/demo/screen/fetch.py": "", "results/demo/top.py": "",
                 "results/demo/tools/build.py": ""}
        found = {command.argv[0]: (command.script, command.cwd) for command in replay.discover(tree_of(files))}
        self.assertEqual(found["fetch.py"], ("results/demo/screen/fetch.py", "results/demo/screen"))
        self.assertEqual(found["top.py"], ("results/demo/top.py", "results/demo"))
        self.assertEqual(found["build.py"], ("results/demo/tools/build.py", "results/demo/tools"))
        from_root = "```sh\ncd results/demo\npython top.py --check\n```\n"
        (command,) = replay.discover(tree_of({"results/demo/README.md": from_root, "results/demo/top.py": ""}))
        self.assertEqual((command.script, command.cwd), ("results/demo/top.py", "results/demo"))

    def test_placeholders_are_recognised(self):
        for token in ("<observation-id>", "PATH_TO_FITS_DIRECTORY", "/path/to/cache", "JOB_ID", "..."):
            self.assertTrue(replay.PLACEHOLDER.search(token), token)
        for token in ("data/ews", "2024-BLG-0062", "--steps", "8000", "DE440", "README.md", "A_b"):
            self.assertFalse(replay.PLACEHOLDER.search(token), token)

    def test_python_blocks_and_module_invocations_are_not_commands(self):
        text = "```python\npython = 3\n```\n\n`python -m pip install x` and `python -c 'print(1)'`\n"
        self.assertEqual(replay.discover(tree_of({"results/demo/README.md": text})), [])

    def test_imports_are_split_into_own_modules_and_packages(self):
        files = {"results/demo/a.py": "import json, helper, numpy\nfrom scipy import stats\nfrom lib.deep import thing\n",
                 "results/demo/helper.py": "", "results/demo/lib/deep.py": ""}
        self.assertEqual(replay.third_party_imports(tree_of(files), "results/demo/a.py"), {"numpy", "scipy"})

    def test_numeric_difference_separates_spelling_from_value(self):
        self.assertEqual(replay.numeric_difference(b"a,1.250\nb,2\n", b"a,1.25\r\nb,2.0\r\n"), (True, 0.0))
        same_shape, deviation = replay.numeric_difference(b"x 1.0 y 100.0", b"x 1.0 y 101.0")
        self.assertTrue(same_shape)
        self.assertAlmostEqual(deviation, 1 / 101)
        self.assertEqual(replay.numeric_difference(b"a 1 b 2", b"a 1 c 2")[0], False)
        self.assertEqual(replay.changed_numbers(b"1.0 2.0 3.0", b"1.0 2.0000000001 3.3"), (1, 3))

    def test_gzip_members_are_compared_by_content(self):
        import gzip
        first = gzip.compress(b"a,b\n1,2\n", mtime=1)
        second = gzip.compress(b"a,b\n1,2\n", mtime=999)
        self.assertNotEqual(first, second)
        self.assertEqual(replay._payload("t.csv.gz", first), replay._payload("t.csv.gz", second))


class SandboxTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tree = tree_of(FILES)
        cls.findings, cls.summary = replay.run(cls.tree, timeout=6.0, jobs=1)
        cls.by_script = {finding.data["command"].split()[1]: finding for finding in cls.findings}

    def status(self, script):
        return self.by_script[script].data["status"]

    def test_every_outcome_is_classified(self):
        self.assertEqual(self.status("make.py"), "reproduced")
        self.assertEqual(self.status("check.py"), "reproduced")
        self.assertEqual(self.status("drift.py"), "drift")
        self.assertEqual(self.status("fetch.py"), "needs-network")
        self.assertEqual(self.status("needs_package.py"), "needs-dependency")
        self.assertEqual(self.status("gone.py"), "script-missing")
        self.assertEqual(self.status("broken.py"), "failed")
        self.assertEqual(self.status("slow.py"), "timeout")
        self.assertEqual(self.status("tool.py"), "template")
        self.assertEqual(self.status("summary.py"), "reproduced")
        self.assertEqual(self.status("quiet.py"), "ran")
        self.assertEqual(self.summary["commands"], 11)

    def test_reproduction_names_the_regenerated_files(self):
        self.assertEqual(self.by_script["make.py"].data["regenerated_identically"], ["results/demo/table.csv"])
        self.assertEqual(self.by_script["summary.py"].data["stdout_compared_with"], "results/demo/summary.json")
        self.assertEqual(self.by_script["make.py"].severity, PASS)

    def test_drift_reports_the_file_and_the_size_of_the_change(self):
        finding = self.by_script["drift.py"]
        self.assertEqual(finding.severity, WARN)
        (item,) = finding.data["drifted"]
        self.assertEqual(item["file"], "results/demo/drifting.csv")
        self.assertAlmostEqual(item["largest_relative_deviation"], 0.01 / 1.26)
        self.assertTrue(item["same_layout"])
        self.assertEqual((item["numbers_changed"], item["numbers_compared"]), (1, 1))
        self.assertIn("1 of 1 numbers moved, by up to 7.9e-03", finding.message)

    def test_the_network_is_refused_not_merely_unreachable(self):
        self.assertIn("network access is disabled", self.by_script["fetch.py"].message)
        self.assertEqual(self.by_script["fetch.py"].severity, INFO)

    def test_failures_keep_the_last_line_of_the_error(self):
        self.assertIn("the fit did not converge", self.by_script["broken.py"].message)

    def test_the_sandbox_is_restored_between_commands(self):
        with tempfile.TemporaryDirectory() as scratch:
            sandbox = replay.Sandbox(self.tree, Path(scratch), runner.interpreter())
            commands = {command.argv[0]: command for command in replay.discover(self.tree)}
            self.assertEqual(sandbox.run(commands["drift.py"], 30).status, "drift")
            self.assertEqual((sandbox.root / "results/demo/drifting.csv").read_bytes(),
                             self.tree.read("results/demo/drifting.csv"))
            self.assertEqual(sandbox.run(commands["check.py"], 30).status, "reproduced")


class WhyACommandStoppedTest(unittest.TestCase):
    """A crash is explained from what the sandbox saw the script do, never from where it sits in a document."""

    FILES = {
        "results/demo/README.md": "```sh\npython swallow.py\npython stamps.py\npython listing.py\n"
                                  "python thesis.py thesis.pdf out.json\npython prepare.py /tmp/wide\n"
                                  "python scratch.py /tmp/wide\npython bug.py\npython carryon.py\n"
                                  "python honest.py\npython native.py\npython wrapper.py\npython optional.py\npython netdrift.py\npython required.py\n```\n",
        "results/demo/table.csv": "n\n3\n",
        "results/demo/other.csv": "n\n3\n",
        "results/demo/carryon.py": "from pathlib import Path\nrows = sorted(Path('cache').glob('*.json'))\n"
                                   "Path('table.csv').write_text('n\\n%d\\n' % len(rows))\n",
        "results/demo/honest.py": "from pathlib import Path\nPath('other.csv').write_text('n\\n4\\n')\n",
        "results/demo/native.py": "import ctypes, ctypes.util\n"
                                  "LIB = ctypes.CDLL(ctypes.util.find_library('no_such_library_xyz'))\n",
        "results/demo/wrapper.py": "import native\n",
        "results/demo/required.py": "import ctypes\nctypes.CDLL('libmissing_example_xyz_123.so')\n",
        "results/demo/optional.py": "import ctypes.util\nif False: ctypes.util.find_library('no_such_library_xyz')\nprint('ok')\n",
        "results/demo/netdrift.py": "import urllib.request\nfrom pathlib import Path\ntry: urllib.request.urlopen('http://example.org')\nexcept Exception: pass\nPath('other.csv').write_text('n\\n4\\n')\n",
        "results/demo/inputs/sources.json": "{}\n",
        "results/demo/swallow.py": "import urllib.request\ntry:\n    urllib.request.urlopen('http://example.org/x', timeout=5)\n"
                                   "except Exception:\n    rows = []\nprint(rows[0])\n",
        "results/demo/stamps.py": "from pathlib import Path\nprint(sorted(Path('stamps').glob('*.npz'))[0])\n",
        "results/demo/listing.py": "from pathlib import Path\nprint(sorted(Path('inputs').glob('P*FTZ'))[0])\n",
        "results/demo/thesis.py": "import subprocess, sys\n"
                                  "subprocess.run([sys.executable, '-c', 'raise SystemExit(1)', sys.argv[1]], check=True)\n",
        "results/demo/prepare.py": "import os, sys\nos.makedirs(sys.argv[1], exist_ok=True)\n",
        "results/demo/scratch.py": "import os, sys\nprint(os.listdir(sys.argv[1])[0])\n",
        "results/demo/bug.py": "from pathlib import Path\nrows = Path('inputs/sources.json').read_text().split()\nprint(rows[5])\n",
    }

    @classmethod
    def setUpClass(cls):
        findings, _ = replay.run(tree_of(cls.FILES), timeout=30.0, jobs=1)
        cls.by_script = {finding.data["command"].split()[1]: finding for finding in findings}

    def outcome(self, script):
        finding = self.by_script[script]
        return finding.data["status"], finding.message

    def test_a_swallowed_network_refusal_is_still_a_network_need(self):
        status, message = self.outcome("swallow.py")
        self.assertEqual(status, "needs-network")
        self.assertIn("tried to reach the network", message)
        self.assertIn("IndexError", message)
        self.assertGreaterEqual(self.by_script["swallow.py"].data["network_attempts"], 1)

    def test_searching_for_files_that_are_not_deposited_is_a_missing_input(self):
        status, message = self.outcome("stamps.py")
        self.assertEqual(status, "needs-input")
        self.assertIn("results/demo/stamps", message)
        status, message = self.outcome("listing.py")
        self.assertEqual(status, "needs-input")
        self.assertIn("results/demo/inputs/P*FTZ", message)

    def test_an_argument_that_names_an_absent_file_is_a_missing_input(self):
        status, message = self.outcome("thesis.py")
        self.assertEqual(status, "needs-input")
        self.assertIn("thesis.pdf", message)
        self.assertEqual(self.by_script["thesis.py"].data["undeposited"], ["thesis.pdf"])

    def test_an_empty_scratch_folder_is_a_missing_input(self):
        self.assertEqual(self.outcome("prepare.py")[0], "ran")
        status, message = self.outcome("scratch.py")
        self.assertEqual(status, "needs-input")
        self.assertIn("/tmp/wide", message)

    def test_missing_probe_does_not_hide_drift(self):
        status, message = self.outcome("carryon.py")
        self.assertEqual(status, "drift")
        self.assertEqual(self.by_script["carryon.py"].severity, WARN)
        self.assertIn("results/demo/cache", message)
        self.assertIn("effect on output is unknown", message)

    def test_network_probe_does_not_hide_drift(self):
        status, message = self.outcome("netdrift.py")
        self.assertEqual(status, "drift")
        self.assertEqual(self.by_script["netdrift.py"].severity, WARN)
        self.assertIn("network attempt", message)

    def test_a_difference_with_every_input_present_stays_a_difference(self):
        self.assertEqual(self.outcome("honest.py")[0], "drift")
        self.assertEqual(self.by_script["honest.py"].severity, WARN)

    def test_unexecuted_library_lookup_does_not_block_replay(self):
        self.assertEqual(self.outcome("optional.py")[0], "ran")

    def test_loader_error_identifies_required_library(self):
        status, message = self.outcome("required.py")
        self.assertEqual(status, "needs-dependency")
        self.assertIn("libmissing_example_xyz_123.so", message)

    def test_a_crash_on_deposited_data_stays_a_failure(self):
        status, message = self.outcome("bug.py")
        self.assertEqual(status, "failed")
        self.assertEqual(self.by_script["bug.py"].severity, WARN)


class Mission03CacheTest(unittest.TestCase):
    def test_empty_cache_does_not_overwrite_results(self):
        script = Path(__file__).resolve().parents[1] / "results/mission03-dr15-thermal-ins/analysis/ledger.py"
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            results, cache = root / "results", root / "cache"
            results.mkdir()
            cache.mkdir()
            target = results / "counterparts.csv"
            target.write_bytes(b"preserve this evidence\n")
            done = subprocess.run([sys.executable, str(script), str(results), str(root / "overlap"), str(cache)],
                                  capture_output=True, text=True)
            self.assertNotEqual(done.returncode, 0)
            self.assertIn("required counterpart cache is missing or empty", done.stderr)
            self.assertEqual(target.read_bytes(), b"preserve this evidence\n")


class LedgerAndCommandLineTest(unittest.TestCase):
    FILES = {
        "README.md": "# P\n\nThe catalog below has 1 entries: 1 new measurements.\n\n"
                     "- **Demo.** A delay of 30.6996 seconds. *Measured* · [Details](results/demo/README.md)\n",
        "results/demo/README.md": "# Demo\n\nThe delay is 30.6996 s over 944 windows.\n",
        "results/demo/measurements.json": json.dumps({"delay_s": 30.6996, "windows": 944}),
        "results/demo/claims.json": json.dumps({
            "schema": "open-space/claims/1",
            "claims": [{"id": "windows", "printed": "944 windows", "compute": "json('measurements.json')['windows']"}]}),
    }

    def write(self, folder, files):
        for path, data in files.items():
            target = Path(folder, path)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data.encode("utf-8"))

    def test_ledger_summarises_every_check(self):
        ledger = runner.run_checks(tree_of(self.FILES))
        summary = ledger["summary"]
        self.assertEqual(summary["findings"][FAIL], 0)
        self.assertEqual(summary["numbers"]["precise"], 2)
        self.assertEqual(summary["numbers"]["found"], 2)
        self.assertEqual(summary["numbers"]["catalog"]["precise"], 1)
        self.assertEqual(summary["claims"]["reproduced"], 1)
        (card,) = ledger["results"]
        self.assertEqual(card["labels"], ["Measured"])
        self.assertEqual(card["numbers"]["found_percent"], 100.0)
        self.assertEqual(card["claims"], {"declared": 1, "reproduced": 1, "failed": 0})
        self.assertEqual(ledger["tree"]["files"], 4)
        self.assertEqual(len(ledger["tool"]["source_sha256"]), 64)
        text = report.markdown(ledger)
        self.assertIn("| Claims recomputed independently | 1 of 1 declared claims", text)
        self.assertIn("Nothing. No check failed and no warning is open.", text)
        self.assertEqual(text, report.markdown(runner.run_checks(tree_of(self.FILES))))
        self.assertEqual(runner.exit_code(ledger), 0)
        json.dumps(runner.to_json(ledger, full=True))

    def test_command_line_writes_report_and_ledger_and_sets_the_exit_code(self):
        with tempfile.TemporaryDirectory() as scratch:
            self.write(scratch, self.FILES)
            ledger_path, report_path = Path(scratch, "out", "ledger.json"), Path(scratch, "out", "report.md")
            with contextlib.redirect_stdout(io.StringIO()) as printed, contextlib.redirect_stderr(io.StringIO()):
                code = cli.main([scratch, "--json", str(ledger_path), "--report", str(report_path)])
            self.assertEqual(code, 0)
            self.assertIn("0 failed", printed.getvalue())
            ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
            self.assertEqual(ledger["schema"], runner.LEDGER_SCHEMA)
            self.assertEqual(ledger["tree"]["origin"], "directory")
            self.assertIn("# Verification report", report_path.read_text(encoding="utf-8"))

            stale = dict(self.FILES)
            stale["results/demo/README.md"] = "# Demo\n\nThe delay is 30.6996 s over 945 windows.\n"
            self.write(scratch, stale)
            with contextlib.redirect_stdout(io.StringIO()) as printed, contextlib.redirect_stderr(io.StringIO()):
                code = cli.main([scratch, "--quiet"])
            self.assertEqual(code, 1)
            self.assertEqual(printed.getvalue(), "")
            with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as errors:
                self.assertEqual(cli.main([scratch, "--result", "nope"]), 2)
            self.assertIn("no such result folder", errors.getvalue())

    def test_warnings_fail_only_in_strict_mode(self):
        files = dict(self.FILES)
        files["results/extra/README.md"] = "# Not in the catalog\n"
        ledger = runner.run_checks(tree_of(files), ("catalog",))
        self.assertEqual(ledger["summary"]["findings"][WARN], 1)
        self.assertEqual(runner.exit_code(ledger), 0)
        self.assertEqual(runner.exit_code(ledger, strict=True), 1)
        with self.assertRaises(ValueError):
            runner.run_checks(tree_of(files), ("catalog", "astrology"))


if __name__ == "__main__":
    unittest.main()
