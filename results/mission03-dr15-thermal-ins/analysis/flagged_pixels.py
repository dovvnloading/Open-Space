from astropy.io import fits
from collections import Counter
import json, numpy as np
from pathlib import Path
import argparse


from event_apertures import event_wcs


def audit(ra, dec, directory):
    results = []
    for p in sorted(directory.glob("*[MP]IEVLI*.FTZ")):
        with fits.open(p) as f:
            h = f["EVENTS"].header
            e = f["EVENTS"].data
            w = event_wcs(h)
            x, y = w.all_world2pix([[ra, dec]], 1)[0]
            r = np.hypot(
                (e["X"] - x) * h["REFXCDLT"] * 3600, (e["Y"] - y) * h["REFYCDLT"] * 3600
            )
            pn = h["INSTRUME"] == "EPN"
            m = (
                (r < 15)
                & (e["PI"] >= 200)
                & (e["PI"] < 1000)
                & (e["PATTERN"] <= (4 if pn else 12))
            )
            flags = Counter(e["FLAG"][m].tolist())
            pix = Counter(zip(e["CCDNR"][m], e["RAWX"][m], e["RAWY"][m]))
            top = pix.most_common(3)
            events = e[m]
            keep = [
                ((events["CCDNR"] == a) & (events["RAWX"] == b) & (events["RAWY"] == c))
                for (a, b, c), n in top
            ]
            q = {
                "instrument": h["INSTRUME"],
                "all_flag_count": int(m.sum()),
                "flags_hex": {hex(k): v for k, v in flags.most_common(8)},
                "all_flag_top_pixels": [[int(i) for i in k] + [n] for k, n in top],
                "peak_pixel_energy_eV_percentiles": np.percentile(
                    events["PI"][keep[0]], [0, 25, 50, 75, 100]
                ).tolist()
                if len(top)
                else None,
                "peak_pixel_time_span_s": float(np.ptp(events["TIME"][keep[0]]))
                if len(top)
                else None,
            }
            results.append(q)
    return results


if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Inspect accepted-pattern soft photons before flag filters"
    )
    p.add_argument("ra", type=float)
    p.add_argument("dec", type=float)
    p.add_argument("directory", type=Path)
    a = p.parse_args()
    print(json.dumps(audit(a.ra, a.dec, a.directory), indent=2))
