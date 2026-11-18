"""Check observatory plate bytes and source metadata. Network replay is opt-in."""

import argparse
import csv
import hashlib
import io
import json
import math
from datetime import datetime
from pathlib import Path
import subprocess
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "observatory/assets"
PLATES = {"mars", "sun", "quasar", "dwarf"}


def validate(directory, *, revision="HEAD", remote=False):
    plates = json.loads((directory / "plates.json").read_text(encoding="utf-8"))
    if set(plates) != PLATES:
        raise ValueError("Expected one plate for each featured dossier.")
    for key, plate in plates.items():
        for field in ("name", "tag", "field", "color", "file", "kind", "credit", "source", "description", "sha256"):
            if not isinstance(plate.get(field), str) or not plate[field].strip():
                raise ValueError(f"{key}: missing {field}")
        name = plate["file"]
        if Path(name).name != name or "/" in name or "\\" in name:
            raise ValueError(f"{key}: expected a single image filename")
        if urlparse(plate["source"]).scheme != "https":
            raise ValueError(f"{key}: missing source URL")
        raw = (directory / name).read_bytes()
        if hashlib.sha256(raw).hexdigest() != plate["sha256"]:
            raise ValueError(f"{key}: {name} does not match the manifest digest")
        if key == "dwarf":
            if raw[:4] != b"RIFF" or raw[8:12] != b"WEBP":
                raise ValueError("Expected the credited WebP illustration.")
            if plate.get("archiveUrl") or plate["kind"] != "Artist’s impression":
                raise ValueError("The generated dwarf illustration requires its artwork credit.")
            continue
        url = plate["archiveUrl"]
        origin = urlparse(url)
        hosts = {"mars": "assets.science.nasa.gov", "sun": "stereo-ssc.nascom.nasa.gov",
                 "quasar": "ps1images.stsci.edu"}
        if origin.scheme != "https" or origin.hostname != hosts[key] or raw[:3] != b"\xff\xd8\xff":
            raise ValueError(f"{key}: expected an archival JPEG and HTTPS download URL")
        if remote:
            request = Request(url, headers={"User-Agent": "Open-Space-Observatory/1.0 (science image replay)"})
            with urlopen(request, timeout=45) as response:
                fresh = response.read(len(raw) + 1)
            if fresh != raw:
                raise ValueError(f"{key}: {name} differs from bytes currently returned by the archive")

    # Check the actual request encoded in the download receipt against the
    # deposited A/B coordinates even in an offline build of the static site.
    target_rows = subprocess.check_output(
        ["git", "-C", str(ROOT), "show", f"{revision}:results/j0826-lens-delay/targets.csv"], encoding="utf-8"
    )
    pair = [r for r in csv.DictReader(io.StringIO(target_rows)) if r["system"] == "J0826+7002"]
    if len(pair) != 2 or {r["image"] for r in pair} != {"A", "B"}:
        raise ValueError("Expected the A and B lens positions for J0826+7002")
    query = parse_qs(urlparse(plates["quasar"]["archiveUrl"]).query, strict_parsing=True)
    for axis in ("ra", "dec"):
        midpoint = sum(float(r[axis + "_deg"]) for r in pair) / 2
        if len(query[axis]) != 1 or not math.isclose(float(query[axis][0]), midpoint, rel_tol=0, abs_tol=1e-8):
            raise ValueError(f"PS1 request {axis} differs from the J0826+7002 midpoint")
    for channel, band in (("red", "i"), ("green", "r"), ("blue", "g")):
        if not query[channel][0].endswith(f".stk.{band}.unconv.fits"):
            raise ValueError(f"Unexpected PS1 band in the {channel} channel")
    if query.get("size") != ["160"] or query.get("output_size") != ["512"]:
        raise ValueError("PS1 field dimensions differ from the image credit")

    cor2 = plates["sun"]
    frame = datetime.fromisoformat(cor2["observedAt"].replace("Z", "+00:00"))
    expected = frame.strftime("/browse/%Y/%m/%d/ahead/cor2/1024/%Y%m%d_%H%M%S_n4c2A.jpg")
    if urlparse(cor2["archiveUrl"]).path != expected:
        raise ValueError("COR2 spacecraft, instrument or image time differs from the archive receipt")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=ASSETS, help="Source or built assets directory")
    parser.add_argument("--revision", default="HEAD", help="Git revision containing the research target table")
    parser.add_argument("--remote", action="store_true", help="Also compare the three downloaded archive files")
    args = parser.parse_args()
    validate(args.directory, revision=args.revision, remote=args.remote)
    print("4 plate digests and credits, COR2 image timestamp and instrument, and PS1 pointing and bands passed.")
    if args.remote:
        print("3 archive JPEGs match fresh downloads byte-for-byte.")
