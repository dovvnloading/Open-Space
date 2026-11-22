"""Bring main into the Chandra feature branch and reconcile its single catalog entry."""
import os
from pathlib import Path
import re
import subprocess

import report


def git(*args):
    return subprocess.check_output(["git", *args], text=True)


def sync():
    if (os.environ.get("GITHUB_ACTIONS") != "true" or
            os.environ.get("GITHUB_REF") != "refs/heads/research/chandra-acis-fast-transients"):
        raise RuntimeError("Use the Chandra feature-branch workflow")
    if git("status","--porcelain","--untracked-files=no").strip():
        raise RuntimeError("Catalog synchronization requires a clean tracked worktree")
    subprocess.run(["git","fetch","--no-tags","origin","+refs/heads/main:refs/remotes/origin/main"],check=True)
    if subprocess.run(["git","merge-base","--is-ancestor","origin/main","HEAD"]).returncode == 0:
        print("The feature branch already contains the current main catalog.")
        return
    feature=git("show","HEAD:README.md")
    published=git("show","origin/main:README.md")
    def entries(body):
        return [line for line in body.splitlines() if line.startswith("- ") and report.TAG in line]
    prior=entries(feature)
    if len(prior)!=1 or entries(published):
        raise ValueError("Unexpected number of Chandra entries before synchronization")
    pattern=(r"^- \*\*Chandra ACIS search for minute-to-hour X-ray flares\.\*\* (.*?) "
             r"\*(Candidate|Measured|Null result|Prediction|Correction)\* ·")
    m=re.match(pattern,prior[0])
    section="**Black holes, outbursts and supernovae**\n\n"
    if m is None or published.count(section)!=1:
        raise ValueError("Current catalog layout requires inspection")
    result=subprocess.run(["git","merge","--no-ff","--no-commit","origin/main"])
    conflicts=git("diff","--name-only","--diff-filter=U").splitlines()
    if (result.returncode and conflicts!=["README.md"]) or any(p!="README.md" for p in conflicts):
        raise RuntimeError(f"Unexpected upstream conflicts: {conflicts}")
    Path("README.md").write_text(published.replace(section,section+prior[0]+"\n\n",1))
    report.update_catalog(m[1],m[2])
    subprocess.run(["git","add","--","README.md"],check=True)
    subprocess.run(["git","commit","-qm","Synchronize the current catalog into the Chandra research branch"],check=True)


if __name__=="__main__":
    sync()
