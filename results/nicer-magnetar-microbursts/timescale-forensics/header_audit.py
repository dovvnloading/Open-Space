"""Dump time-reference headers from the actual two original NICER event files.

Run with astropy installed and a verified PR #23 cache:
    python header_audit.py /path/to/nicer-input-cache

Input files are reused from ../inputs.json and ../clock.py. This reads only
original CALIBRATED local EVENTS extensions: not a newly barycentered product.
It also prints GTI time systems and offsets; no values are invented from
mission-wide documentation in place of reading individual headers.
"""
import argparse
from pathlib import Path

KEYS = ("TIMESYS", "TIMEREF", "TASSIGN", "CLOCKAPP", "MJDREFI",
        "MJDREFF", "TIMEZERO", "LEAPINIT", "TSTART", "TSTOP",
        "RA_OBJ", "DEC_OBJ")


def main(cache):
    from astropy.io import fits
    for obs in ("1020240141", "5505010502"):
        filename = cache / f"ni{obs}_0mpu7_cl.evt.gz"
        with fits.open(filename) as hdus:
            print(f"Observation {obs}: {filename.name}")
            for name in ("EVENTS", "GTI"):
                header = hdus[name].header
                print(f"  [{name}]")
                for key in KEYS:
                    if key in header:
                        print(f"    {key}: {header[key]!r}")
                if name == "EVENTS":
                    assert header.get("TIMESYS", "").strip() == "TT"
                    assert header.get("CLOCKAPP") is True
                    assert "TIME" in hdus[name].columns.names
                else:
                    assert "START" in hdus[name].columns.names
                    assert "STOP" in hdus[name].columns.names
        print()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cache", type=Path)
    main(p.parse_args().cache)
