"""Replay core analyses in a new disposable copy and check frozen hashes."""
import argparse
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def verify(root):
    for line in (root / 'SHA256SUMS').read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        path = root / name
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f'Frozen byte mismatch: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('destination', type=Path, help='new disposable directory; must not exist')
    dest = parser.parse_args().destination.resolve()
    if dest.exists() or dest == ROOT or ROOT in dest.parents or dest in ROOT.parents:
        parser.error('destination must be a new directory outside the evidence tree')
    verify(ROOT)
    shutil.copytree(ROOT, dest, ignore=shutil.ignore_patterns('__pycache__', '*.pyc'))
    for args in [('analysis/check.py',), ('analysis/run.py',), ('analysis/run.py', '--archival'),
                 ('analysis/exposure.py',), ('analysis/supplement.py',)]:
        subprocess.run([sys.executable, *args], cwd=dest, check=True)
    verify(dest)
    verify(ROOT)
    print('Core replay completed; copied and original frozen files match SHA256SUMS.')


if __name__ == '__main__':
    main()