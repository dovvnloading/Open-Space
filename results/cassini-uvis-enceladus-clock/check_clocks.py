"""Reconstruct the Cassini UVIS time labels from PDS and the two NAIF kernels."""

import argparse
import json
from pathlib import Path
import re
from xml.etree import ElementTree

import spiceypy as spice


HERE = Path(__file__).resolve().parent
ARCHIVE = HERE / "archive"
ZERO = "2017-03-27T14:30:00"
KERNELS = ("cas00171.tsc", "cas00172.tsc")
EVENT = ("FUV2017_086_14_30.LBL", "fuv2017_086_14_30.lblx", "HSP2017_086_14_30.LBL")
CONTROLS = [
    (f"FUV2016_071_{cut}_PDS3.LBL", f"FUV2016_071_{cut}_PDS4.lblx", f"HSP2016_071_{cut}_PDS3.LBL")
    for cut in ("11_48", "11_52")
]


def odl(path, key):
    text = (ARCHIVE / path).read_text(encoding="ascii", errors="replace")
    value = re.search(rf"^{key}\s*=\s*([^\r\n]+)", text, re.MULTILINE)
    if value is None:
        raise ValueError(f"Missing {key} in {path}")
    return value.group(1).strip().strip('"')


def xml(path, key):
    tree = ElementTree.parse(ARCHIVE / path)
    return next(item.text for item in tree.iter() if item.tag.rsplit("}", 1)[-1] == key)


def instant(t):
    return spice.str2et(t.removesuffix("Z"))


def relative(a, b):
    return round(float(a - b), 9)


def normalized_count(raw):
    prefix, ticks = raw.split(".")
    fine = int(ticks)
    if not 0 <= fine <= 255:
        raise ValueError("Invalid Cassini fine clock component")
    return f"{prefix}.{fine:03}"


def calculate():
    event_clock = {}
    instances = []
    for index, (fuv, migrated, hsp) in enumerate((EVENT, *CONTROLS)):
        pair = {"fuv_pds3_label": odl(fuv, "START_TIME"), "fuv_pds4_label": xml(migrated, "start_date_time"),
                "hsp_pds3_label": odl(hsp, "START_TIME"), "fuv_sclk": normalized_count(odl(fuv, "SPACECRAFT_CLOCK_START_COUNT")),
                "hsp_sclk": normalized_count(odl(hsp, "SPACECRAFT_CLOCK_START_COUNT"))}
        migrated_clock = normalized_count("1/" + xml(migrated, "spacecraft_clock_start_count"))
        if migrated_clock != pair["fuv_sclk"]:
            raise ValueError("FUV migrated count differs from PDS3")
        for kernel in KERNELS:
            spice.kclear()
            for file in ("naif0012.tls", kernel):
                spice.furnsh(str(ARCHIVE / file))
            target = {key: instant(pair[key]) for key in ("fuv_pds3_label", "fuv_pds4_label", "hsp_pds3_label")}
            zero = instant(ZERO)
            a = spice.scs2e(-82, pair["fuv_sclk"])
            b = spice.scs2e(-82, pair["hsp_sclk"])
            pair[kernel] = {
                "fuv_recomputed_utc": spice.et2utc(a, "ISOC", 9),
                "hsp_recomputed_utc": spice.et2utc(b, "ISOC", 9),
                "legacy_fuv_minus_recomputed_s": relative(target["fuv_pds3_label"], a),
                "migrated_fuv_minus_recomputed_s": relative(target["fuv_pds4_label"], a),
                "legacy_hsp_minus_recomputed_s": relative(target["hsp_pds3_label"], b),
                "recomputed_fuv_minus_hsp_s": relative(a, b),
            }
            if index == 0 and kernel == KERNELS[-1]:
                event_clock = {
                    "pds3_pair": [relative(target["fuv_pds3_label"], zero), relative(target["hsp_pds3_label"], zero)],
                    "mixed_archive": [relative(target["fuv_pds4_label"], zero), relative(target["hsp_pds3_label"], zero)],
                    "sclk_pair": [relative(a, zero), relative(b, zero)],
                }
        pair["legacy_fuv_minus_hsp_labels_s"] = relative(target["fuv_pds3_label"], target["hsp_pds3_label"])
        pair["mixed_fuv_minus_hsp_labels_s"] = relative(target["fuv_pds4_label"], target["hsp_pds3_label"])
        pair["legacy_minus_migrated_fuv_labels_s"] = relative(target["fuv_pds3_label"], target["fuv_pds4_label"])
        instances.append(pair)
    spice.kclear()
    naive = instances[0]["mixed_fuv_minus_hsp_labels_s"]
    consistent = instances[0][KERNELS[-1]]["recomputed_fuv_minus_hsp_s"]
    return {
        "event": instances[0], "earlier_controls": instances[1:],
        "first_end_seconds_since_utc_143000": event_clock,
        "mixed_archive_channel_shift_magnitude_s": round(abs(naive - consistent), 6),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = calculate()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(result["mixed_archive_channel_shift_magnitude_s"])


if __name__ == "__main__":
    main()
