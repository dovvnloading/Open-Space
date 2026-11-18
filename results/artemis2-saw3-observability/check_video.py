"""Read native SAW3 MP4 presentation times and decoded luma statistics."""

from argparse import ArgumentParser
import csv
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import file_digest, sha256
import json
from pathlib import Path
from xml.etree import ElementTree

import av
import numpy as np

from range_source import RangeSource


FIRST_CUE = "2026-04-07T00:41:01Z"
SECOND_CUE = "2026-04-07T00:51:15Z"
ONBOARD_START = "2026-04-07T01:19:13Z"
LABEL_NS = {"pds": "http://pds.nasa.gov/pds4/pds/v1"}


def utc(value):
    return datetime.fromisoformat(value.replace("Z", "+00:00")).astimezone(timezone.utc)


def native_input(source_dir, item, use_ranges):
    label = Path(__file__).with_name(item["label"]).read_bytes()
    if sha256(label).hexdigest() != item["label_sha256"]:
        raise ValueError(f"PDS label checksum mismatch: {item['label']}")
    xml = ElementTree.fromstring(label)
    start = utc(xml.findtext(".//pds:Time_Coordinates/pds:start_date_time", namespaces=LABEL_NS))
    stop = utc(xml.findtext(".//pds:Time_Coordinates/pds:stop_date_time", namespaces=LABEL_NS))
    file = source_dir / item["file"]
    if use_ranges and (not file.is_file() or file.stat().st_size != item["bytes"]):
        return RangeSource(item["url"], item["bytes"], source_dir / file.stem), start, stop
    with file.open("rb") as stream:
        md5 = file_digest(stream, "md5").hexdigest()
    if file.stat().st_size != item["bytes"] or md5 != item["md5"]:
        raise ValueError(f"PDS video checksum or size mismatch: {file.name}")
    return str(file), start, stop


def luma_pixels(frame):
    if frame.format.name not in ("yuv420p", "yuvj420p"):
        raise ValueError(f"Unexpected video pixel format: {frame.format.name}")
    plane = frame.planes[0]
    return np.frombuffer(plane, dtype=np.uint8).reshape(frame.height, plane.line_size)[:, : frame.width]


def luma_summary(image, previous):
    height, width = image.shape
    inner = np.s_[height // 10 : height - height // 10, width // 10 : width - width // 10]
    interior = image[inner]
    cdf = np.bincount(interior.ravel(), minlength=256).cumsum()
    levels = np.searchsorted(cdf, interior.size * np.array([0.01, 0.50, 0.99])).tolist()
    fraction_dark = float(cdf[16] / interior.size)
    fraction_changed = None
    exact_repeat = None
    if previous is not None:
        exact_repeat = bool(np.array_equal(image, previous))
        change = np.abs(np.subtract(interior, previous[inner], dtype=np.int16))
        fraction_changed = float(np.count_nonzero(change > 2) / interior.size)
    return {"luma_p01": levels[0], "luma_p50": levels[1], "luma_p99": levels[2],
            "fraction_luma_at_most_16": fraction_dark, "exact_repeat": exact_repeat,
            "fraction_interior_changed_gt2": fraction_changed}


def stream_windows(source, windows, output_dir):
    result = {}
    with av.open(source, mode="r") as video:
        stream = video.streams.video[0]
        origin = stream.start_time or 0
        details = {
            "frame_rate": str(stream.average_rate), "time_base": str(stream.time_base),
            "start_pts": origin, "duration": float(stream.duration * stream.time_base),
            "width": stream.width, "height": stream.height, "codec": stream.codec_context.name,
            "container_duration": video.duration / av.time_base,
        }
        for name, left, right in windows:
            video.seek(origin + int(Fraction(left) / stream.time_base),
                       stream=stream, backward=True, any_frame=False)
            previous = None
            rows = []
            for frame in video.decode(stream):
                if frame.pts is None:
                    raise ValueError("Video frame has no presentation timestamp")
                time = (frame.pts - origin) * frame.time_base
                if time < left:
                    continue
                if time >= right:
                    break
                pixels = luma_pixels(frame)
                rows.append(dict(pts=frame.pts, elapsed_seconds=round(float(time), 9),
                                 **luma_summary(pixels, previous)))
                if len(rows) == 1:
                    # Preserve one exact luma frame for each frozen interval.
                    (output_dir / f"{name}.pgm").write_bytes(
                        f"P5\n{frame.width} {frame.height}\n255\n".encode("ascii") + pixels.tobytes())
                previous = pixels.copy()
            if not rows:
                raise ValueError(f"No frames decoded for {name}")
            times = np.array([row["elapsed_seconds"] for row in rows])
            if np.any(np.diff(times) <= 0):
                raise ValueError(f"Presentation times are not increasing in {name}")
            changes = np.array([row["fraction_interior_changed_gt2"] is not None
                                and row["fraction_interior_changed_gt2"] > 0.001 for row in rows])
            update_intervals = np.diff(times[changes])
            result[name] = {
                "requested_elapsed": [left, right], "first_elapsed": float(times[0]),
                "last_elapsed": float(times[-1]), "frames": len(rows),
                "repeats": sum(row["exact_repeat"] is True for row in rows),
                "low_change_pairs": len(rows) - 1 - int(changes.sum()),
                "interior_updates": int(changes.sum()),
                "median_update_interval_seconds": float(np.median(update_intervals)) if len(update_intervals) else None,
                "median_interior_luma": float(np.median([r["luma_p50"] for r in rows])),
                "max_interior_p99": max(r["luma_p99"] for r in rows),
                "median_fraction_luma_at_most_16": float(np.median([r["fraction_luma_at_most_16"] for r in rows])),
                "largest_pts_gap_seconds": float(np.max(np.diff(times))),
            }
            with (output_dir / f"{name}.csv").open("x", encoding="utf-8", newline="") as out:
                writer = csv.DictWriter(out, fieldnames=list(rows[0]), lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
    return details, result


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--ranges", action="store_true")
    args = parser.parse_args()
    source = json.loads(args.manifest.read_bytes())
    args.output_dir.mkdir()
    result = {"decoder": {"pyav": av.__version__, "libraries": av.library_versions}}
    receipts = {}
    for name, item in source.items():
        file, start, stop = native_input(args.source_dir, item, args.ranges)
        span = (stop - start).total_seconds()
        cues = [(utc(cue) - start).total_seconds() for cue in (FIRST_CUE, SECOND_CUE, ONBOARD_START)]
        if name == "downlink":
            windows = [("first-cue", cues[0]-15, cues[0]+15), ("second-cue", cues[1]-15, cues[1]+15),
                       ("same-time-downlink", cues[2], cues[2]+30)]
        else:
            windows = [("onboard", 0, span)]
        details, stats = stream_windows(file, windows, args.output_dir)
        result[name] = {"pds_start_utc": start.isoformat(), "pds_stop_utc": stop.isoformat(),
                        "pds_span_seconds": span, "decoder": details, "windows": stats,
                        "cue_to_start_seconds": cues[:2]}
        receipts[name] = file.receipt() if isinstance(file, RangeSource) else {"verified_md5": item["md5"]}
    for name, payload in (("video-results.json", result), ("acquisition.json", receipts)):
        with (args.output_dir / name).open("x", encoding="utf-8", newline="\n") as out:
            out.write(json.dumps(payload, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
