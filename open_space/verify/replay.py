"""Replay: run the commands each README documents and compare with what is committed.

Every ``python ...`` command found in a result folder's Markdown is run in a
sandbox that holds an exact copy of the committed tree. The network is switched
off inside the sandbox, so the outcome does not depend on any archive being
reachable, and a command that needs a download says so instead of hanging.

After a command exits, every file it changed is compared with the committed
bytes. A check that passes, or outputs that regenerate identically, count as a
reproduction. Outputs that regenerate with different numbers are drift. A
command that cannot run offline is recorded with the reason: it needs the
network, a package that is not installed, or an input that was never deposited.
"""

from __future__ import annotations

import ast
from concurrent.futures import ThreadPoolExecutor
import gzip
from dataclasses import dataclass, field
import hashlib
import math
import os
from pathlib import Path
import posixpath
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import threading
import time

from . import markdown
from .model import FAIL, INFO, PASS, WARN, Finding
from .tree import RepoTree

CHECK = "replay"
MARKER = "OPEN_SPACE_VERIFY: network access is disabled during replay"
SHELL_LANGUAGES = {"", "sh", "bash", "shell", "console", "zsh", "text", "powershell", "ps1", "bat", "cmd"}
IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg", ".svg", ".pdf", ".gif")
# Containers whose bytes carry timestamps or compression choices: rewriting one proves nothing either way.
OPAQUE_SUFFIXES = (".npz", ".npy", ".zip", ".fits", ".fit", ".h5", ".hdf5", ".pkl", ".pickle", ".parquet", ".sqlite")
CHECK_WORDS = ("check", "verify", "validate", "replay", "reproduce", "test", "audit", "holdout")
NUMBER = re.compile(r"[+\-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+\-]?\d+)?")
TOLERANCE = 1e-6        # largest relative change in any number that still counts as the same output
PLACEHOLDER = re.compile(r"<[^<>]+>|\.\.\.|/path/to|^[A-Z][A-Z0-9]*(?:_[A-Z0-9]+)+$")

BLOCKER = f'''"""Installed by open_space.verify for the duration of a replay: no network."""
import os as _os
import socket as _socket

_LOG = _os.environ.get("OPEN_SPACE_VERIFY_NETLOG")


def _refuse(*arguments, **keywords):
    if _LOG:
        try:
            with open(_LOG, "a", encoding="utf-8") as _handle:
                _handle.write("refused\\n")
        except OSError:
            pass
    raise OSError({MARKER!r})


_socket.socket.connect = _refuse
_socket.socket.connect_ex = _refuse
_socket.create_connection = _refuse
_socket.getaddrinfo = _refuse

# What the script looked for inside the tree and did not find: a file it tried to read, a folder it
# tried to list, a pattern that matched nothing. Recorded so that a later crash can be explained.
import sys as _sys

_LOOKS = _os.environ.get("OPEN_SPACE_VERIFY_LOOKLOG")
_TREE = _os.path.normcase(_os.path.abspath(_os.environ.get("OPEN_SPACE_VERIFY_TREE", "")))
_WRITING = _os.O_WRONLY | _os.O_RDWR | _os.O_CREAT
_busy = False


def _inside(path):
    try:
        full = _os.path.abspath(_os.fsdecode(path))
    except (TypeError, ValueError):
        return None
    return full if _os.path.normcase(full).startswith(_TREE + _os.sep) else None


def _watch(event, arguments):
    global _busy
    if _busy or event not in ("open", "os.scandir", "os.listdir", "glob.glob", "pathlib.Path.glob",
                              "pathlib.Path.rglob"):
        return
    _busy = True
    try:
        found = None
        if event == "open":
            path, mode, flags = (tuple(arguments) + (None, None))[:3]
            reading = not any(letter in mode for letter in "wax+") if isinstance(mode, str) \\
                else not (isinstance(flags, int) and flags & _WRITING)
            full = None if isinstance(path, int) or not reading else _inside(path)
            if full and not _os.path.exists(full):
                found = "file\\t" + full
        elif event in ("os.scandir", "os.listdir"):
            full = _inside(arguments[0] if arguments and arguments[0] is not None else ".")
            if full and (not _os.path.isdir(full) or not _os.listdir(full)):
                found = "folder\\t" + full
        elif event == "glob.glob":
            import glob as _glob
            full = _inside(arguments[0])
            if full and not _glob.glob(full, recursive=bool(arguments[1]) if len(arguments) > 1 else False):
                found = "pattern\\t" + full
        else:
            base, pattern = arguments[0], str(arguments[1])
            full = _inside(base)
            walk = base.rglob if event.endswith("rglob") else base.glob
            if full and next(iter(walk(pattern)), None) is None:
                found = "pattern\\t" + _os.path.join(full, pattern)
        if found:
            with open(_LOOKS, "a", encoding="utf-8") as _handle:
                _handle.write(found + "\\n")
    except Exception:
        pass
    finally:
        _busy = False


if _LOOKS and _TREE:
    _sys.addaudithook(_watch)
'''

NETWORK_SIGNS = (MARKER, "URLError", "ConnectionError", "ConnectionRefusedError", "gaierror", "ProxyError",
                 "Max retries exceeded", "Temporary failure in name resolution", "getaddrinfo failed",
                 "RemoteDisconnected", "ConnectTimeout", "HTTPError", "Network is unreachable",
                 "Command '['curl'", "Command '['wget'", "Could not resolve host", "Failed to connect to")
INPUT_SIGNS = ("FileNotFoundError", "No such file or directory", "the following arguments are required",
               "NotADirectoryError", "IsADirectoryError", "does not exist", "not found in cache", "usage: ")
FILE_NAME = re.compile(r"[^/\\.]\.[A-Za-z][A-Za-z0-9]{0,4}$")     # an argument that names a file
PACKAGE_SIGN = re.compile(r"(?:ModuleNotFoundError|ImportError): No module named '([^']+)'"
                          r"|[`']?[Ii]mport ([A-Za-z_][A-Za-z0-9_]*)[`']? failed"
                          r"|Missing optional dependency '([^']+)'")


@dataclass
class Command:
    """One documented command, resolved to a script inside the tree."""

    argv: list[str]                 # as written, starting with the script path
    text: str                       # the command as it appears in the document
    document: str
    line: int
    result: str
    script: str | None = None       # tree path of the script, or None when it does not exist
    cwd: str = ""                   # tree folder the command runs in
    redirect: str | None = None     # file the document redirects stdout to, if any

    @property
    def key(self) -> str:
        return f"{self.cwd}$ {' '.join(self.argv)}"


@dataclass
class Outcome:
    command: Command
    status: str                     # reproduced, ran, drift, failed, needs-network, needs-dependency, needs-input, timeout
    detail: str = ""
    seconds: float = 0.0
    changed: list[str] = field(default_factory=list)
    data: dict = field(default_factory=dict)


# ---------------------------------------------------------------- finding commands


def _python_lines(block: markdown.CodeBlock):
    """Yield (line number, command text, folder) for python invocations in a code block or code span.

    ``folder`` is where earlier ``cd`` lines of the same block have moved to, relative to the document.
    """
    if block.fenced and block.language not in SHELL_LANGUAGES:
        return
    pending = ""
    folder = ""
    for offset, raw in enumerate(block.content.splitlines() if block.fenced else [block.content]):
        line = raw.strip()
        if pending:
            line = pending + " " + line
            pending = ""
        if line.endswith("\\"):
            pending = line[:-1].strip()
            continue
        line = re.sub(r"^(?:\$|>|PS>)\s+", "", line)
        moved = re.match(r"^cd\s+([^\s;&|]+)\s*(?:&&\s*(.*))?$", line)
        if moved:
            target = moved.group(1).strip("\"'").replace("\\", "/")
            if not target.startswith(("/", "~", "$")) and ":" not in target:
                folder = posixpath.normpath(posixpath.join(folder, target))
                folder = "" if folder == "." else folder
            line = (moved.group(2) or "").strip()
        if re.match(r"^(?:python3?|py)(?:\.exe)?\s+\S", line):
            yield block.line + offset, line, folder


def _split(command: str) -> tuple[list[str], str | None]:
    """Tokenise one shell line, keep only the python invocation, note a stdout redirect."""
    head = re.split(r"\s(?:&&|\|\||;|\|)\s", command, maxsplit=1)[0]
    try:
        tokens = shlex.split(head, posix=True)
    except ValueError:
        return [], None
    redirect = None
    kept: list[str] = []
    skip = False
    for index, token in enumerate(tokens):
        if skip:
            skip = False
            continue
        if token in (">", "1>", ">>"):
            redirect = tokens[index + 1] if index + 1 < len(tokens) else None
            skip = True
            continue
        if token.startswith(">") and len(token) > 1:
            redirect = token.lstrip(">")
            continue
        if token in ("2>", "2>&1", "<") or token.startswith("2>"):
            skip = token in ("2>", "<")
            continue
        kept.append(token)
    return kept[1:], redirect


def discover(tree: RepoTree, only: set[str] | None = None) -> list[Command]:
    """Every distinct documented python command in the result folders."""
    commands: dict[str, Command] = {}
    for folder in tree.result_dirs():
        if only is not None and folder not in only:
            continue
        root = f"results/{folder}"
        for path in tree.under(root):
            if not path.lower().endswith(".md"):
                continue
            _, blocks = markdown.split_code(tree.text(path))
            for block in blocks:
                for line, text, moved in _python_lines(block):
                    argv, redirect = _split(text)
                    if not argv or argv[0].startswith("-"):
                        continue                      # python -m ..., python -c ...: not a script in this folder
                    command = Command(argv, text, path, line, folder, redirect=redirect)
                    written = argv[0].replace("\\", "/")
                    starts = (posixpath.dirname(path), root, "")
                    bases: list[str] = []
                    if moved:                         # a `cd` may be meant from the document or from the repository root
                        for start in starts:
                            target = posixpath.normpath(posixpath.join(start, moved))
                            target = "" if target == "." else target
                            if not target.startswith("..") and (target == "" or tree.has_dir(target)) \
                                    and target not in bases:
                                bases.append(target)
                    bases += [start for start in starts if start not in bases]
                    for base in bases:
                        candidate = posixpath.normpath(posixpath.join(base, written))
                        if candidate in tree and not candidate.startswith(".."):
                            command.script, command.cwd = candidate, base
                            break
                    if command.script is None:
                        command.cwd = bases[0]
                    commands.setdefault(command.key, command)
    return sorted(commands.values(), key=lambda command: (command.result, command.document, command.line))


# ---------------------------------------------------------------- static profile


def third_party_imports(tree: RepoTree, script: str) -> set[str]:
    """Top-level modules a script imports that are neither standard library nor files beside it."""
    try:
        module = ast.parse(tree.text(script))
    except (SyntaxError, ValueError, RecursionError):
        return set()
    names: set[str] = set()
    for node in ast.walk(module):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names.add(node.module.split(".")[0])
    standard = set(sys.stdlib_module_names)
    parts = script.split("/")
    home = "/".join(parts[:2]) if len(parts) > 2 and parts[0] == "results" else posixpath.dirname(script)
    siblings = tree.under(home) if home else tree.paths
    local = set()
    for name in names - standard:
        # A module that lives anywhere in the same result folder, or at the repository root, is the project's own.
        if f"{name}.py" in tree or tree.has_dir(name) or any(
                path.endswith(f"/{name}.py") or f"/{name}/" in path for path in siblings):
            local.add(name)
    return names - standard - local


def missing_modules(python: str, names: set[str]) -> set[str]:
    """Which of ``names`` the target interpreter cannot import."""
    if not names:
        return set()
    probe = ("import importlib.util, sys\n"
             "print('\\n'.join(n for n in sys.argv[1:] if importlib.util.find_spec(n) is None))")
    try:
        done = subprocess.run([python, "-c", probe, *sorted(names)], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return set(names)
    return {line.strip() for line in done.stdout.splitlines() if line.strip()}


# ---------------------------------------------------------------- comparing outputs


def numeric_difference(before: bytes, after: bytes) -> tuple[bool, float]:
    """Compare two text files as token streams. Returns (same shape, largest relative deviation).

    The text between the numbers must be identical; the numbers may differ. A deviation
    of 0.0 with the same shape means the files differ only in how a number is spelled.
    """
    try:
        old, new = before.decode("utf-8"), after.decode("utf-8")
    except UnicodeDecodeError:
        return False, math.inf
    old, new = old.replace("\r\n", "\n"), new.replace("\r\n", "\n")
    if NUMBER.sub("#", old).split() != NUMBER.sub("#", new).split():
        return False, math.inf
    worst = 0.0
    for a, b in zip(NUMBER.findall(old), NUMBER.findall(new)):
        x, y = float(a), float(b)
        if x == y:
            continue
        scale = max(abs(x), abs(y))
        worst = max(worst, abs(x - y) / scale if scale else 0.0)
    return True, worst


def changed_numbers(before: bytes, after: bytes, tolerance: float = TOLERANCE) -> tuple[int, int]:
    """For two files of the same layout: (numbers that moved by more than the tolerance, numbers compared)."""
    old = NUMBER.findall(before.decode("utf-8", "replace"))
    new = NUMBER.findall(after.decode("utf-8", "replace"))
    moved = 0
    for a, b in zip(old, new):
        x, y = float(a), float(b)
        scale = max(abs(x), abs(y))
        if scale and abs(x - y) / scale > tolerance:
            moved += 1
    return moved, min(len(old), len(new))


def _payload(path: str, data: bytes) -> bytes:
    """What a file holds: gzip members are compared by content, since their header carries a timestamp."""
    if path.lower().endswith(".gz") and data[:2] == b"\x1f\x8b":
        try:
            return gzip.decompress(data)
        except (OSError, EOFError, ValueError):
            return data
    return data


# ---------------------------------------------------------------- sandbox


class Sandbox:
    """An exact copy of the tree in a scratch directory, with the network switched off."""

    def __init__(self, tree: RepoTree, base: Path, python: str) -> None:
        self.tree = tree
        self.root = base / "tree"
        self.scratch = base / "tmp"
        self.blocker = base / "blocker"
        self.netlog = base / "network-attempts.log"
        self.looklog = base / "looked-for.log"
        self.python = python
        for folder in (self.root, self.scratch, self.blocker):
            folder.mkdir(parents=True, exist_ok=True)
        (self.blocker / "sitecustomize.py").write_text(BLOCKER, encoding="utf-8")
        for path in tree.paths:
            target = self.root / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(tree.read(path))
        self._stamp = self._snapshot()

    def _snapshot(self) -> dict[str, tuple[int, int]]:
        state: dict[str, tuple[int, int]] = {}
        for folder, subdirs, names in os.walk(self.root):
            subdirs[:] = [name for name in subdirs if name != "__pycache__"]
            for name in names:
                full = os.path.join(folder, name)
                try:
                    info = os.stat(full)
                except OSError:
                    continue
                state[os.path.relpath(full, self.root).replace(os.sep, "/")] = (info.st_size, info.st_mtime_ns)
        return state

    def environment(self) -> dict[str, str]:
        env = {key: value for key, value in os.environ.items()
               if key.upper() not in ("PYTHONPATH", "PYTHONSTARTUP", "PYTHONHOME", "NO_PROXY")}
        dead = "http://127.0.0.1:9"
        env.update({
            "PYTHONPATH": os.pathsep.join([str(self.blocker), str(self.root)]),
            "PYTHONDONTWRITEBYTECODE": "1", "PYTHONHASHSEED": "0", "PYTHONIOENCODING": "utf-8", "PYTHONUTF8": "1",
            "MPLBACKEND": "Agg", "TZ": "UTC", "OPEN_SPACE_VERIFY_REPLAY": "1",
            "OPEN_SPACE_VERIFY_NETLOG": str(self.netlog), "OPEN_SPACE_VERIFY_LOOKLOG": str(self.looklog),
            "OPEN_SPACE_VERIFY_TREE": str(self.root),
            "HOME": str(self.scratch), "USERPROFILE": str(self.scratch),
            "TMPDIR": str(self.scratch), "TEMP": str(self.scratch), "TMP": str(self.scratch),
            "http_proxy": dead, "https_proxy": dead, "ftp_proxy": dead, "all_proxy": dead,
            "HTTP_PROXY": dead, "HTTPS_PROXY": dead, "FTP_PROXY": dead, "ALL_PROXY": dead,
        })
        return env

    def run(self, command: Command, timeout: float) -> Outcome:
        cwd = self.root / command.cwd if command.cwd else self.root
        argv = [self.python]
        unfilled = []                   # scratch paths the command names that nothing has filled
        undeposited = []                # files the command names that are not in the tree
        for position, token in enumerate(command.argv):
            if token.startswith("/tmp/") or token == "/tmp":
                mapped = self.scratch / token[5:] if len(token) > 5 else self.scratch
                if not mapped.exists() or (mapped.is_dir() and not any(mapped.iterdir())):
                    unfilled.append(token)
                token = str(mapped)
            elif position and not token.startswith("-") and FILE_NAME.search(token) and not (cwd / token).exists():
                undeposited.append(token)
            argv.append(token)
        for log in (self.netlog, self.looklog):
            try:
                log.unlink()
            except OSError:
                pass
        started = time.monotonic()
        kwargs = dict(cwd=str(cwd), env=self.environment(), stdin=subprocess.DEVNULL,
                      stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if os.name == "nt":
            kwargs["creationflags"] = getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)
        else:
            kwargs["start_new_session"] = True
        try:
            process = subprocess.Popen(argv, **kwargs)
        except OSError as error:
            return Outcome(command, "failed", f"could not start: {error}")
        try:
            out, err = process.communicate(timeout=timeout)
            timed_out = False
        except subprocess.TimeoutExpired:
            _kill(process)
            out, err = process.communicate()
            timed_out = True
        seconds = time.monotonic() - started
        changed = self._changes()
        outcome = _classify(self.tree, command, process.returncode, out, err, timed_out, changed, self.root)
        try:
            attempts = len(self.netlog.read_text(encoding="utf-8").splitlines())
        except OSError:
            attempts = 0
        if outcome.status == "timeout" and attempts:
            outcome = Outcome(command, "needs-network", f"kept trying to reach the network ({attempts} "
                              f"{'attempt' if attempts == 1 else 'attempts'}) until the time limit",
                              data=dict(outcome.data, network_attempts=attempts))
        elif outcome.status == "failed":
            said = (err + b"\n" + out).decode("utf-8", "replace")
            named = [token for token in undeposited if posixpath.basename(token.replace("\\", "/")) in said]
            if attempts:
                # The script swallowed the refusal and stopped later, on the data it never received.
                outcome = Outcome(command, "needs-network", f"tried to reach the network ({attempts} "
                                  f"{'attempt' if attempts == 1 else 'attempts'}) and then stopped: " + outcome.detail,
                                  data=dict(outcome.data, network_attempts=attempts))
            elif unfilled:
                # It reads a scratch folder that an earlier download step would have filled.
                outcome = Outcome(command, "needs-input", "reads " + ", ".join(unfilled) + ", which only a download "
                                  "step fills: " + outcome.detail, data=outcome.data)
            elif named:
                outcome = Outcome(command, "needs-input", "names " + ", ".join(named) + ", which is not deposited: "
                                  + outcome.detail, data=dict(outcome.data, undeposited=named))
            else:
                looked = self._looked_for()
                if looked:
                    # It went looking for inputs that are not in the tree and crashed on their absence.
                    outcome = Outcome(command, "needs-input", "looked for " + ", ".join(looked[:3])
                                      + (f" and {len(looked) - 3} more" if len(looked) > 3 else "")
                                      + ", not deposited: " + outcome.detail,
                                      data=dict(outcome.data, looked_for=looked[:20]))
        outcome.seconds = round(seconds, 2)
        self._restore(changed)
        return outcome

    def _looked_for(self) -> list[str]:
        """Tree paths and patterns the last command searched for without finding anything."""
        try:
            lines = self.looklog.read_text(encoding="utf-8").splitlines()
        except OSError:
            return []
        looked: list[str] = []
        for line in lines:
            kind, _, full = line.partition("\t")
            try:
                relative = os.path.relpath(full, self.root).replace(os.sep, "/")
            except ValueError:
                continue
            if relative.startswith("..") or relative in self.tree or "__pycache__" in relative                     or posixpath.basename(relative).startswith("<"):
                continue                               # "<string>", "<unknown>": the traceback printer's probes
            if kind == "folder":
                relative += "/"
            if relative not in looked:
                looked.append(relative)
        return looked

    def _changes(self) -> list[str]:
        now = self._snapshot()
        changed = [path for path, stamp in now.items() if self._stamp.get(path) != stamp]
        changed += [path for path in self._stamp if path not in now]
        return sorted(changed)

    def _restore(self, changed: list[str]) -> None:
        for path in changed:
            target = self.root / path
            if path in self.tree:
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(self.tree.read(path))
            elif target.exists():
                try:
                    target.unlink()
                except OSError:
                    pass
        self._stamp = self._snapshot()


def _kill(process: subprocess.Popen) -> None:
    try:
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(process.pid)], capture_output=True)
        else:
            os.killpg(os.getpgid(process.pid), signal.SIGKILL)
    except (OSError, ProcessLookupError):
        process.kill()


def _last_line(text: str) -> str:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    return lines[-1][:300] if lines else ""


def _classify(tree: RepoTree, command: Command, code: int | None, out: bytes, err: bytes, timed_out: bool,
              changed: list[str], root: Path) -> Outcome:
    stderr = err.decode("utf-8", "replace")
    stdout = out.decode("utf-8", "replace")
    if timed_out:
        return Outcome(command, "timeout", "did not finish within the time limit", changed=changed)
    if code != 0:
        both = stderr + "\n" + stdout
        last = _last_line(stderr) or _last_line(stdout)
        missing = PACKAGE_SIGN.search(both)
        if missing:
            package = next(group for group in missing.groups() if group).split(".")[0]
            return Outcome(command, "needs-dependency", f"needs the package {package}", data={"modules": [package]})
        if any(sign in both for sign in NETWORK_SIGNS):
            return Outcome(command, "needs-network", "tried to reach the network: " + last)
        if any(sign in both for sign in INPUT_SIGNS):
            return Outcome(command, "needs-input", "needs an input or argument that is not deposited: " + last)
        return Outcome(command, "failed", f"exit code {code}: {last}", changed=changed, data={"exit_code": code})
    identical, close, drifted, images, opaque, created = [], [], [], [], [], []
    worst = 0.0
    for path in changed:
        if path not in tree:
            created.append(path)
            continue
        target = root / path
        after = _payload(path, target.read_bytes() if target.exists() else b"")
        before = _payload(path, tree.read(path))
        if after == before or after.replace(b"\r\n", b"\n") == before.replace(b"\r\n", b"\n"):
            identical.append(path)                     # line endings follow the platform that wrote the file
            continue
        if path.lower().endswith(IMAGE_SUFFIXES):
            images.append(path)
            continue
        if path.lower().endswith(OPAQUE_SUFFIXES):
            opaque.append(path)
            continue
        same_shape, deviation = numeric_difference(before, after)
        if same_shape and deviation <= TOLERANCE:
            close.append(path)
            worst = max(worst, deviation)
        elif same_shape:
            moved, compared = changed_numbers(before, after)
            drifted.append({"file": path, "same_layout": True, "largest_relative_deviation": deviation,
                            "numbers_changed": moved, "numbers_compared": compared})
        else:
            drifted.append({"file": path, "same_layout": False, "largest_relative_deviation": None})
    data = {"regenerated_identically": identical, "regenerated_within_tolerance": close,
            "images_rewritten": images, "binary_files_rewritten": opaque, "new_files": len(created)}
    if drifted:
        data["drifted"] = drifted
        parts = []
        for item in drifted[:4]:
            if item["same_layout"]:
                parts.append(f"{posixpath.basename(item['file'])} ({item['numbers_changed']} of "
                             f"{item['numbers_compared']} numbers moved, by up to "
                             f"{item['largest_relative_deviation']:.1e} of their value)")
            else:
                parts.append(posixpath.basename(item["file"]) + " (different text or a different number of values)")
        return Outcome(command, "drift", "regenerated outputs differ from the committed ones: " + ", ".join(parts),
                       changed=changed, data=data)
    if command.redirect and stdout.strip():
        written = command.redirect.replace("\\", "/")
        for base in (command.cwd, ""):
            candidate = posixpath.normpath(posixpath.join(base, written))
            if candidate in tree:
                same_shape, deviation = numeric_difference(tree.read(candidate), out)
                data["stdout_compared_with"] = candidate
                if not (same_shape and deviation <= TOLERANCE):
                    data["drifted"] = [{"file": candidate, "same_layout": same_shape, "largest_relative_deviation":
                                        None if math.isinf(deviation) else deviation}]
                    return Outcome(command, "drift", f"standard output differs from committed {candidate}",
                                   changed=changed, data=data)
                identical.append(candidate)
                break
    stdout_digest = hashlib.sha256(out).hexdigest()
    twins = tree.by_digest().get(stdout_digest)
    if twins and out.strip():
        data["stdout_equals"] = twins[0]
        identical.append(twins[0])
    words = " ".join(command.argv).lower()
    is_check = any(word in words for word in CHECK_WORDS)
    if identical or close:
        note = f"{len(identical) + len(close)} committed outputs regenerated"
        note += " identically" if not close or worst == 0 else f"; numbers agree to {worst:.1e} of their value"
        data["largest_relative_deviation"] = worst
        return Outcome(command, "reproduced", note, changed=changed, data=data)
    if is_check:
        return Outcome(command, "reproduced", "the documented check passed offline", changed=changed, data=data)
    return Outcome(command, "ran", "exit code 0; it changed no committed file, so there was nothing to compare",
                   changed=changed, data=data)


# ---------------------------------------------------------------- driver

SEVERITY = {
    "reproduced": PASS, "ran": INFO, "drift": WARN, "failed": WARN, "needs-network": INFO,
    "needs-dependency": INFO, "needs-input": INFO, "timeout": INFO, "script-missing": WARN, "template": INFO,
}


def _finding(command: Command, status: str, detail: str, data: dict | None = None, seconds: float | None = None) -> Finding:
    payload = {"command": command.text, "status": status}
    payload.update(data or {})
    if seconds is not None:
        payload["seconds"] = seconds
    return Finding(CHECK, f"replay/{status}", SEVERITY[status], command.document, f"`{command.text}`: {detail}",
                   line=command.line, key=command.key, result=command.result, data=payload)


def run(tree: RepoTree, only: set[str] | None = None, *, python: str | None = None, timeout: float = 120.0,
        jobs: int = 2, limit: int | None = None, progress=None) -> tuple[list[Finding], dict]:
    """Replay every documented command. Returns findings and a summary."""
    python = python or sys.executable
    located = shutil.which(python) or (str(Path(python).resolve()) if Path(python).exists() else None)
    if located is None:
        raise FileNotFoundError(f"interpreter {python!r} was not found")
    python = str(Path(located).resolve()) if os.sep in located or "/" in located else located
    commands = discover(tree, only)
    findings: list[Finding] = []
    runnable: list[Command] = []
    needed: dict[str, set[str]] = {}
    for command in commands:
        if command.script is None:
            findings.append(_finding(command, "script-missing",
                                     f"the script {command.argv[0]} does not exist in the repository"))
            continue
        if any(PLACEHOLDER.search(token) for token in command.argv[1:]):
            findings.append(_finding(command, "template", "contains a placeholder, so it was not run"))
            continue
        needed[command.key] = third_party_imports(tree, command.script)
        runnable.append(command)
    absent = missing_modules(python, set().union(*needed.values()) if needed else set())
    queue: list[Command] = []
    for command in runnable:
        lacking = sorted(needed[command.key] & absent)
        if lacking:
            findings.append(_finding(command, "needs-dependency",
                                     "needs packages that are not installed here: " + ", ".join(lacking),
                                     {"modules": lacking}))
        else:
            queue.append(command)
    if limit is not None:
        queue = queue[:limit]
    outcomes: list[Outcome] = []
    if queue:
        base = Path(tempfile.mkdtemp(prefix="osv-"))
        try:
            # One result folder is replayed by one worker, in document order, so its commands never
            # race each other. Workers take the next waiting folder as they become free.
            groups: dict[str, list[Command]] = {}
            for command in queue:
                groups.setdefault(command.result, []).append(command)
            waiting = sorted(groups.values(), key=len, reverse=True)
            jobs = max(1, min(jobs, len(waiting)))
            lock = threading.Lock()

            def work(index: int) -> list[Outcome]:
                sandbox = None
                done = []
                while True:
                    with lock:
                        if not waiting:
                            return done
                        group = waiting.pop(0)
                    if sandbox is None:
                        sandbox = Sandbox(tree, base / f"w{index}", python)
                    for command in group:
                        outcome = sandbox.run(command, timeout)
                        done.append(outcome)
                        if progress:
                            progress(outcome)

            with ThreadPoolExecutor(max_workers=jobs) as pool:
                for batch in pool.map(work, range(jobs)):
                    outcomes.extend(batch)
            order = {command.key: position for position, command in enumerate(queue)}
            outcomes.sort(key=lambda outcome: order[outcome.command.key])
        finally:
            shutil.rmtree(base, ignore_errors=True)
    for outcome in outcomes:
        findings.append(_finding(outcome.command, outcome.status, outcome.detail, outcome.data, outcome.seconds))
    counts: dict[str, int] = {}
    per_result: dict[str, dict[str, int]] = {}
    for finding in findings:
        status = finding.data["status"]
        counts[status] = counts.get(status, 0) + 1
        bucket = per_result.setdefault(finding.result, {})
        bucket[status] = bucket.get(status, 0) + 1
    summary = {"commands": len(commands), "by_status": dict(sorted(counts.items())), "per_result": per_result,
               "python": python, "timeout_seconds": timeout,
               "seconds": round(sum(outcome.seconds for outcome in outcomes), 1)}
    return findings, summary
