"""Verify the exact historical implementation named in the original campaign receipt.

python results/chandra_acis_fxt/archive_attestation.py
python results/chandra_acis_fxt/archive_attestation.py --rebuild

Rebuilding uses the original Git objects from a checkout containing the pinned
commit. --source-dir can instead supply the nine files directly; every input
must match its independently recorded hash. No archive or network access.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import subprocess
import zipfile

from migrate import BEFORE, COMMIT, RUN

HERE = Path(__file__).resolve().parent
PREFIX = "results/chandra_acis_fxt"
SNAPSHOT = HERE / "source_history" / f"{COMMIT}.zip"


def checked_sources(read_file):
    sources = {}
    for name, expected in sorted(BEFORE.items()):
        data = read_file(name)
        actual = hashlib.sha256(data).hexdigest()
        if actual != expected:
            raise ValueError(f"{name}: historical source digest {actual} differs from {expected}")
        sources[f"{PREFIX}/{name}"] = data
    return sources


def manifest():
    return (json.dumps({"repository":"dovvnloading/Open-Space","git_commit":COMMIT,
                        "workflow_run":RUN,"files":[{"path":f"{PREFIX}/{name}","sha256":digest}
                        for name,digest in sorted(BEFORE.items())]},indent=2)+"\n").encode()


def verify(path=SNAPSHOT):
    with zipfile.ZipFile(path) as archive:
        expected = {"manifest.json", *(f"{PREFIX}/{name}" for name in BEFORE)}
        if len(archive.namelist()) != len(expected) or set(archive.namelist()) != expected:
            raise ValueError("Historical source ZIP has missing, extra or duplicate entries")
        if archive.read("manifest.json") != manifest():
            raise ValueError("Historical source ZIP provenance differs")
        checked_sources(lambda name: archive.read(f"{PREFIX}/{name}"))


def rebuild(source_dir=None):
    def original(name):
        if source_dir is not None:
            return (source_dir/name).read_bytes()
        return subprocess.check_output(["git","show",f"{COMMIT}:{PREFIX}/{name}"],cwd=HERE)
    members = checked_sources(original)
    members["manifest.json"] = manifest()
    container = io.BytesIO()
    with zipfile.ZipFile(container,"w",compression=zipfile.ZIP_STORED) as archive:
        for name,data in sorted(members.items()):
            info = zipfile.ZipInfo(name,date_time=(1980,1,1,0,0,0))
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            archive.writestr(info,data)
    container.seek(0)
    verify(container)
    SNAPSHOT.parent.mkdir(parents=True,exist_ok=True)
    staged = SNAPSHOT.with_suffix(".tmp")
    staged.write_bytes(container.getvalue())
    staged.replace(SNAPSHOT)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild",action="store_true")
    parser.add_argument("--source-dir",type=Path)
    args=parser.parse_args()
    if args.source_dir and not args.rebuild:
        parser.error("--source-dir is an alternative input for --rebuild")
    if args.rebuild:
        rebuild(args.source_dir)
    verify()
    print(f"Verified {len(BEFORE)} frozen source files from {COMMIT}.")
