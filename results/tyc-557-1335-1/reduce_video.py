"""Extract aperture flux from Sposetti's display video, omitting its slow-motion replay."""
import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path

import numpy as np

BASE = Path(__file__).resolve().parent
STARS = {"target": (317, 249), "ref1": (348, 70), "ref2": (57, 128), "ref3": (263, 202)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("video", type=Path)
    parser.add_argument("--output", type=Path, default=BASE / "inputs/video_photometry.csv")
    args = parser.parse_args()
    expected = json.loads((BASE / "inputs/sources.json").read_text())["video"]["sha256"]
    if hashlib.sha256(args.video.read_bytes()).hexdigest() != expected:
        raise SystemExit("Video differs from the pinned source.")
    # Frame 0 is stamped 20:13:49.3489022 UTC and frame 120 is 20:13:51.7675897.
    # Only the first uninterrupted, normal-speed segment contains independent exposures.
    decoded = subprocess.check_output([
        "ffmpeg", "-v", "error", "-i", str(args.video), "-frames:v", "120",
        "-pix_fmt", "gray", "-f", "rawvideo", "-",
    ])
    frames = np.frombuffer(decoded, dtype="u1").reshape(120, 300, 480)
    dy, dx = np.mgrid[-8:9, -8:9]
    distance = np.hypot(dy, dx)
    data = {"frame": np.arange(len(frames))}
    for name, (x, y) in STARS.items():
        patch = frames[:, y-8:y+9, x-8:x+9]
        sky = np.median(patch[:, (distance >= 5) & (distance <= 8)], axis=1)
        for radius in (2, 3, 4):
            aperture = distance <= radius
            data[f"{name}_r{radius}"] = patch[:, aperture].sum(axis=1) - sky * aperture.sum()
    with args.output.open("w") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(data)
        writer.writerows(zip(*data.values()))


if __name__ == "__main__":
    main()
