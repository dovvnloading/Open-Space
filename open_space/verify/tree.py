"""Exact-byte views of the files under verification.

Checksums only mean something for the bytes that were committed. A checkout can
differ from them (Git rewrites line endings on some platforms), so the default
view reads blobs straight from the Git object store. Every check works on a
``RepoTree`` and never touches the filesystem itself.
"""

from __future__ import annotations

import bisect
import hashlib
import os
from pathlib import Path
import subprocess
import tempfile

SKIP_DIRS = {".git", ".venv", "venv", "runs", "__pycache__", "build", "dist", "node_modules",
             ".mypy_cache", ".pytest_cache", ".ruff_cache"}


class TreeError(RuntimeError):
    """The requested revision or directory could not be read."""


def _git(repo: Path, *args: str, stdin: bytes | None = None, env: dict | None = None) -> bytes:
    try:
        done = subprocess.run(["git", "-C", str(repo), *args], input=stdin, capture_output=True, env=env)
    except FileNotFoundError as error:
        raise TreeError("git is not installed") from error
    if done.returncode != 0:
        detail = done.stderr.decode("utf-8", "replace").strip().splitlines()
        raise TreeError(f"git {' '.join(args[:2])} failed: {detail[-1] if detail else done.returncode}")
    return done.stdout


def is_git_repo(root: Path) -> bool:
    try:
        return _git(Path(root), "rev-parse", "--is-inside-work-tree").strip() == b"true"
    except TreeError:
        return False


class RepoTree:
    """An immutable mapping from POSIX path to the exact bytes under verification."""

    def __init__(self, files: dict[str, bytes], *, root: Path | None = None, origin: str = "memory",
                 commit: str | None = None, git_tree: str | None = None) -> None:
        self._files = files
        self.paths = sorted(files)
        self.root = Path(root) if root is not None else None
        self.origin = origin            # "commit", "worktree", "directory" or "memory"
        self.commit = commit
        self.git_tree = git_tree
        self._digests: dict[str, str] = {}
        self._by_digest: dict[str, list[str]] | None = None

    # ------------------------------------------------------------------ construction

    @classmethod
    def from_git(cls, repo: Path | str, rev: str = "HEAD", *, origin: str = "commit") -> "RepoTree":
        """Read every blob of ``rev`` from the object store. No checkout conversion applies."""
        repo = Path(repo)
        try:
            commit = _git(repo, "rev-parse", "--verify", "--quiet", rev + "^{commit}").decode().strip()
        except TreeError:
            commit = None                # ``rev`` may name a bare tree (the worktree snapshot does)
        git_tree = _git(repo, "rev-parse", "--verify", rev + "^{tree}").decode().strip()
        entries: list[tuple[bytes, str]] = []
        for item in _git(repo, "ls-tree", "-r", "-z", "--full-tree", git_tree).split(b"\0"):
            if not item:
                continue
            meta, _, raw_path = item.partition(b"\t")
            mode, kind, sha = meta.split()
            if kind != b"blob" or mode == b"120000":      # submodules and symbolic links carry no evidence
                continue
            entries.append((sha, raw_path.decode("utf-8", "surrogateescape")))
        packed = _git(repo, "cat-file", "--batch", stdin=b"".join(sha + b"\n" for sha, _ in entries))
        files: dict[str, bytes] = {}
        cursor = 0
        for sha, path in entries:
            end = packed.index(b"\n", cursor)
            header = packed[cursor:end].split()
            if len(header) != 3 or header[1] != b"blob":
                raise TreeError(f"unexpected object header for {path}: {packed[cursor:end]!r}")
            size = int(header[2])
            files[path] = packed[end + 1:end + 1 + size]
            cursor = end + 1 + size + 1
        return cls(files, root=repo, origin=origin, commit=commit, git_tree=git_tree)

    @classmethod
    def from_worktree(cls, repo: Path | str) -> "RepoTree":
        """Snapshot what ``git add -A && git commit`` would record, without touching the real index.

        The snapshot is staged into a throwaway index file. Git stores the new blobs as
        unreferenced loose objects, which its garbage collection removes later.
        """
        repo = Path(repo)
        with tempfile.TemporaryDirectory(prefix="open-space-verify-") as scratch:
            env = dict(os.environ, GIT_INDEX_FILE=os.path.join(scratch, "index"))
            try:
                _git(repo, "read-tree", "HEAD", env=env)
            except TreeError:
                pass                     # a repository without commits still has a working tree
            _git(repo, "add", "-A", env=env)
            snapshot = _git(repo, "write-tree", env=env).decode().strip()
        tree = cls.from_git(repo, snapshot, origin="worktree")
        try:
            tree.commit = _git(repo, "rev-parse", "--verify", "--quiet", "HEAD^{commit}").decode().strip()
        except TreeError:
            tree.commit = None
        return tree

    @classmethod
    def from_directory(cls, root: Path | str) -> "RepoTree":
        """Read a plain directory. Bytes are taken as found, including any converted line endings."""
        root = Path(root)
        if not root.is_dir():
            raise TreeError(f"{root} is not a directory")
        files: dict[str, bytes] = {}
        for folder, subdirs, names in os.walk(root):
            subdirs[:] = sorted(d for d in subdirs if d not in SKIP_DIRS and not d.endswith(".egg-info"))
            for name in sorted(names):
                full = Path(folder) / name
                if full.is_symlink() or not full.is_file():
                    continue
                files[full.relative_to(root).as_posix()] = full.read_bytes()
        return cls(files, root=root, origin="directory")

    @classmethod
    def open(cls, root: Path | str, *, rev: str | None = None, worktree: bool = False) -> "RepoTree":
        """Pick the most exact view available: a named revision, the worktree snapshot or HEAD."""
        root = Path(root)
        if is_git_repo(root):
            if worktree:
                return cls.from_worktree(root)
            return cls.from_git(root, rev or "HEAD")
        if rev:
            raise TreeError(f"{root} is not a Git repository, so revision {rev!r} cannot be read")
        return cls.from_directory(root)

    # ------------------------------------------------------------------ access

    def __contains__(self, path: str) -> bool:
        return path in self._files

    def __len__(self) -> int:
        return len(self.paths)

    def read(self, path: str) -> bytes:
        return self._files[path]

    def text(self, path: str) -> str:
        """Decode for parsing. A byte-order mark is dropped; undecodable bytes are replaced."""
        return self._files[path].decode("utf-8-sig", "replace")

    def size(self, path: str) -> int:
        return len(self._files[path])

    def digest(self, path: str) -> str:
        found = self._digests.get(path)
        if found is None:
            found = self._digests[path] = hashlib.sha256(self._files[path]).hexdigest()
        return found

    def by_digest(self) -> dict[str, list[str]]:
        """Map each SHA-256 to the paths holding those bytes (content addressing)."""
        if self._by_digest is None:
            index: dict[str, list[str]] = {}
            for path in self.paths:
                index.setdefault(self.digest(path), []).append(path)
            self._by_digest = index
        return self._by_digest

    def under(self, prefix: str) -> list[str]:
        """Paths inside a directory, found by bisection on the sorted listing."""
        prefix = prefix.rstrip("/") + "/"
        start = bisect.bisect_left(self.paths, prefix)
        stop = bisect.bisect_left(self.paths, prefix[:-1] + "0")     # "0" sorts directly after "/"
        return self.paths[start:stop]

    def has_dir(self, prefix: str) -> bool:
        prefix = prefix.rstrip("/") + "/"
        at = bisect.bisect_left(self.paths, prefix)
        return at < len(self.paths) and self.paths[at].startswith(prefix)

    def result_dirs(self) -> list[str]:
        """Names of the folders directly under results/."""
        seen: list[str] = []
        for path in self.under("results"):
            parts = path.split("/")
            if len(parts) >= 3 and (not seen or seen[-1] != parts[1]):
                seen.append(parts[1])
        return seen

    # ------------------------------------------------------------------ identity and history

    def tree_digest(self) -> str:
        """SHA-256 of the tree's own ``SHA256SUMS`` listing, so anyone can recompute it."""
        listing = hashlib.sha256()
        for path in self.paths:
            listing.update(f"{self.digest(path)}  {path}\n".encode("utf-8", "surrogateescape"))
        return listing.hexdigest()

    def identity(self) -> dict:
        return {
            "origin": self.origin,
            "commit": self.commit,
            "git_tree": self.git_tree,
            "tree_sha256": self.tree_digest(),
            "files": len(self.paths),
            "bytes": sum(len(data) for data in self._files.values()),
        }

    def history(self, path: str, limit: int = 30) -> list[tuple[str, str]]:
        """Earlier committed versions of ``path`` as (commit, sha256), newest first.

        Used to explain a stale checksum. Returns nothing outside Git or in a shallow clone.
        """
        if self.root is None or self.origin not in ("commit", "worktree") or not self.commit:
            return []
        try:
            commits = _git(self.root, "log", f"-n{limit}", "--format=%H", self.commit, "--", path).decode().split()
        except TreeError:
            return []
        versions: list[tuple[str, str]] = []
        for commit in commits:
            try:
                blob = _git(self.root, "cat-file", "blob", f"{commit}:{path}")
            except TreeError:
                continue
            versions.append((commit, hashlib.sha256(blob).hexdigest()))
        return versions
