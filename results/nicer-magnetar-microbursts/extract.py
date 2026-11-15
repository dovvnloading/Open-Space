"""Extract the NICER event and instrument comparisons from hash-pinned archives."""
import argparse
import csv
import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
MICRO = ["1020240141", "5505010502"]
STORM = "3020560101"
# Raw spacecraft event TIME in seconds for identifying the three longer
# excesses. A +/-5 s bracket is wider than the two source-position-dependent
# barycentric conversions and the catalog durations. The more conservative
# full-GTI pointing test below does not use these times.
LONG_TIMES = {"5505011602": 299075147, "5505011703": 302222214, "5505011704": 302300261}


def write_rows(directory, name, rows):
    with (directory / name).open("w", newline="") as out:
        writer = csv.DictWriter(out, list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def inputs(cache, download):
    cache.mkdir(parents=True, exist_ok=True)
    for name, info in json.loads((HERE / "inputs.json").read_text()).items():
        path = cache / name
        if not path.exists() and download:
            temporary = cache / (name + ".partial")
            with urllib.request.urlopen(info["url"], timeout=180) as source, temporary.open("wb") as out:
                shutil.copyfileobj(source, out)
            with temporary.open("rb") as fp:
                if hashlib.file_digest(fp, "sha256").hexdigest() != info["sha256"]:
                    raise ValueError(f"Archive hash changed: {name}")
            temporary.replace(path)
        with path.open("rb") as fp:
            if hashlib.file_digest(fp, "sha256").hexdigest() != info["sha256"]:
                raise ValueError(f"Archive hash changed: {name}")


def flags(table):
    # FITS 8X storage, high-order bit first. Reject missing, forced and reset
    # flags separately from whether the slow and fast measurement chains fired.
    return np.recarray.field(table, "EVENT_FLAGS").reshape(-1)


def band(pi):
    return (pi >= 50) & (pi <= 800)


def clusters(times):
    """Merge overlapping five-count, 10-microsecond windows."""
    indices = np.flatnonzero(times[4:] - times[:-4] <= 10e-6)
    if not len(indices):
        return []
    groups = np.split(indices, np.flatnonzero(np.diff(indices) > 4) + 1)
    return [(group[0], group[-1] + 5) for group in groups]


def position(text, hours=False):
    major, minute, second = map(float, text.split())
    value = abs(major) + minute / 60 + second / 3600
    return np.copysign(value * (15 if hours else 1), major)


def separation(ra, dec, source_ra, source_dec):
    ra, dec, r0, d0 = [np.deg2rad(np.asarray(v, dtype=float)) for v in (ra, dec, source_ra, source_dec)]
    hav = np.sin((dec - d0) / 2)**2 + np.cos(d0) * np.cos(dec) * np.sin((ra - r0) / 2)**2
    return np.rad2deg(2 * np.arcsin(np.sqrt(hav))) * 60


def main(cache, output):
    output.mkdir(parents=True, exist_ok=True)
    event_rows, matching, scan_rows, pointing, hard = [], [], [], [], []
    summary = {}
    footprint = []
    with (cache / "mcgill.csv").open() as source:
        sources = {r["Name"]: r for r in csv.DictReader(source)}
    catalog = {}
    for line in (cache / "pubtable.txt").read_text().splitlines():
        obs = line[26:36]
        if obs in [*MICRO, *LONG_TIMES]:
            catalog[obs] = (line[:21].strip(), line[37:60], float(line[61:69]))
    if len(catalog) != 5:
        raise ValueError("Expected all four SMC entries and the Galactic microburst")

    for obs, (name, label, t90) in catalog.items():
        cl = cache / f"ni{obs}_0mpu7_cl.evt.gz"
        with fits.open(cl) as hdus:
            d = hdus["EVENTS"].data
            h = hdus["EVENTS"].header
            gti = hdus["GTI"]
            t, pi = d["TIME"], d["PI"]
            det, f = d["DET_ID"], flags(d)
            tzero = h.get("TIMEZERO", 0)
            good = band(pi)
            if obs in MICRO:
                times = t[good]
                clumps = clusters(times)
                if len(clumps) != 1:
                    raise ValueError(f"{obs}: expected one cleaned cluster, found {len(clumps)}")
                i, j = clumps[0]
                center = times[i]
                sig = np.abs(t-center) < 10e-6
                comparison = (np.abs(t-center) > 1) & good
                scan_rows.append(dict(obsid=obs, product="clean", eligible=len(times),
                                      events=j-i, first_time=repr(float(center))))
                if obs=="1020240141":
                    for det_id in np.unique(det[comparison]):
                        ii = np.flatnonzero(comparison & (det==det_id))
                        footprint.append(dict(det_id=int(det_id), rawx=int(d["RAWX"][ii[0]]),
                                              rawy=int(d["RAWY"][ii[0]]), count=len(ii)))
                for row in d[sig & good]:
                    det_id, channel = int(row["DET_ID"]), int(row["PI"])
                    pool = np.flatnonzero(comparison & (det == det_id))
                    delta = np.abs(pi[pool].astype(int) - channel)
                    use = pool[np.argsort(delta, kind="stable")[:20]]
                    matching.append(dict(obsid=obs, det_id=det_id, pi=channel,
                                         controls=len(use), fast=int(np.sum((f[use] & 8) != 0)),
                                         min_pi=int(pi[use].min()), max_pi=int(pi[use].max()),
                                         within25pct=int(np.sum(np.abs(pi[use].astype(int)-channel)<=channel*0.25))))
                for scope, take in [("flash", sig), ("rest", np.abs(t-center)>1)]:
                    for lo in (50,100,150,200,300):
                        m = take & (pi >= lo) & (pi <= 800)
                        hard.append(dict(obsid=obs, sample=scope, min_pi=lo, max_pi=800,
                                         events=int(m.sum()), fast=int(np.sum((f[m]&8)!=0))))
                with fits.open(cache / f"ni{obs}_0mpu7_ufa.evt.gz") as raw:
                    r = raw["EVENTS"].data
                    ft = flags(r)
                    rt = r["TIME"] + (raw["EVENTS"].header.get("TIMEZERO", 0)-tzero)
                    mask = (ft & 16 != 0) & (ft & 7 == 0) & band(r["PI"])
                    tt = rt[mask]
                    for i, j in clusters(tt):
                        scan_rows.append(dict(obsid=obs, product="UFA", eligible=len(tt),
                                              events=j-i, first_time=repr(float(tt[i]))))
                    keep = np.flatnonzero(np.abs(rt-center) <= 500e-6)
                    clean_keys = set(zip(det[np.abs(t-center)<=500e-6], t[np.abs(t-center)<=500e-6]))
                    for k in keep:
                        row = r[k]
                        event_rows.append(dict(obsid=obs, time=repr(float(rt[k])),
                            delta_us=round(float(rt[k]-center)*1e6,5), det_id=int(row["DET_ID"]),
                            rawx=int(row["RAWX"]), rawy=int(row["RAWY"]), flags=int(ft[k]),
                            pi=int(row["PI"]), pi_fast=int(row["PI_FAST"]),
                            pi_ratio=round(float(row["PI_RATIO"]),5) if np.isfinite(row["PI_RATIO"]) else "",
                            deadtime_us=round(float(row["DEADTIME"])*1e6,3),
                            cleaned=int((row["DET_ID"], rt[k]) in clean_keys)))
            else:
                center = LONG_TIMES[obs]

            with fits.open(cache / f"ni{obs}.mkf.gz") as mkf:
                m = mkf["PREFILTER"].data
                mt = m["TIME"] + mkf["PREFILTER"].header.get("TIMEZERO", 0)
                gtiz = gti.header.get("TIMEZERO", 0)
                use = np.zeros(len(mt), dtype=bool)
                for start, stop in gti.data:
                    use |= (mt >= start + gtiz) & (mt < stop + gtiz)
                src = sources[name]
                ra, dec = position(src["RA"], True), position(src["Decl"])
                sep = separation(m["RA"], m["DEC"], ra, dec)
                k = int(np.argmin(np.abs(mt-(center+tzero))))
                near = np.abs(mt - (center+tzero)) <= 5
                pointing.append(dict(obsid=obs, catalog_source=name, catalog_start=label,
                    catalog_t90_ms=t90, fits_target=h["OBJECT"], source_ra=round(ra,8),
                    source_dec=round(dec,8), spacecraft_event_time=repr(float(center)),
                    file_timezero=tzero, ra=round(float(m["RA"][k]),6), dec=round(float(m["DEC"][k]),6),
                    offset_arcmin=round(float(sep[k]),4), bracket_min=round(float(sep[near].min()),4),
                    bracket_max=round(float(sep[near].max()),4), gti_min=round(float(sep[use].min()),4),
                    gti_max=round(float(sep[use].max()),4), gti_samples=int(use.sum()),
                    cutoff_rigidity=round(float(m["COR_SAX"][k]),4),
                    sunlight=int(m["SUNSHINE"][k]), enabled=int(m["NUM_FPM_ON"][k])))

    for filename, rows in [("events.csv",event_rows),("matched_controls.csv",matching),
                           ("pointing.csv",pointing),("cluster_scan.csv",scan_rows),("chain_counts.csv",hard),("detector_counts_2018.csv", footprint)]:
        write_rows(output, filename, rows)

    with fits.open(cache / f"ni{STORM}_0mpu7_cl.evt.gz") as hdus:
        d = hdus["EVENTS"].data
        mask = band(d["PI"])
        bins, n = np.unique(np.floor(d["TIME"][mask]*10).astype(int), return_counts=True)
        dense = bins[n>=100]
        use = np.isin(np.floor(d["TIME"]*10).astype(int), dense) & mask
        fast = (flags(d) & 8) != 0
        hard = use & (d["PI"]>=200)
        summary["burst_storm"] = dict(obsid=STORM, bins=int(len(dense)),
            counts=int(use.sum()),fast=int(np.sum(use&fast)),hard=int(hard.sum()),hard_fast=int(np.sum(hard&fast)))
    # CALDB maps tabulate a finite angular grid at one energy. Report only
    # sampled grid cells and use the ratio to the same shell's on-axis value.
    with fits.open(cache / "nixtivignette20170601v001.fits") as hdus:
        d, h = hdus[1].data, hdus[1].header
        v = d["VIGNET"]
        axes = [(np.arange(1,n+1)-h[f"{i}CRPX6"])*h[f"{i}CDLT6"]+h[f"{i}CRVL6"]
                for i,n in [(1,v.shape[2]),(2,v.shape[1])]]
        radius = np.hypot(*np.meshgrid(*axes))
        m = radius>=720
        v0 = v[:,int(h["2CRPX6"]-1),int(h["1CRPX6"]-1)]
        summary["vignetting"] = dict(energy_kev=float(d["ENERG_LO"][0]), profiles=len(v),
             grid_points_ge_12_arcmin=int(m.sum()), largest_ratio=float((v[:,m]/v0[:,None]).max()))
    (output/"comparison.json").write_text(json.dumps(summary,indent=2)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path,help="Directory for the original input files")
    parser.add_argument("output",type=Path,help="Directory for the derived tables")
    parser.add_argument("--download",action="store_true",help="Retrieve missing inputs from their archive URLs")
    args=parser.parse_args()
    inputs(args.cache,args.download)
    main(args.cache,args.output)
