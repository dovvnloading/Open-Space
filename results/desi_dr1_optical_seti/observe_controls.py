"""Force each persistent candidate on the target, adjacent slit traces and nearest sky.

Run after scan_shards.py. Reads public cframes in memory; saves hashes and
measurements for each hit. The row index is resolved from the physical FIBER ID.
"""
import argparse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess

from astropy.io import fits
import numpy as np

from vetting import forced_line, select_controls

BASE = "https://desidata.s3.amazonaws.com/dr1/spectro/redux/iron/exposures/"


def download(url):
    return subprocess.run(["curl", "-fsSL", "--retry", "3", "--max-time", "450",
                           url], check=True, stdout=subprocess.PIPE).stdout


def extract(hdu, data, entry):
    targetid, exposure, hit, arm = entry
    fiber = exposure["fiber"]
    ids, fib = data["FIBER"], data
    camera = hdu[0].header["CAMERA"]
    expid = int(hdu[0].header["EXPID"])
    if expid != exposure["expid"] or camera != arm + str(exposure["petal_loc"]):
        raise ValueError("wrong camera/exposure")
    rowmap = {int(f): i for i, f in enumerate(ids)}
    if int(fib["TARGETID"][rowmap[fiber]]) != targetid:
        raise ValueError("target ID does not match fibre")
    controls = select_controls(ids, np.char.strip(fib["OBJTYPE"].astype(str)) == "SKY",
                               fib["FIBER_RA"], fib["FIBER_DEC"], fiber, nsky=len(ids))
    wave = hdu["WAVELENGTH"].data.astype(float)
    center = hit["observed_A"]
    left = np.searchsorted(wave, center, side="right") - 1
    t = (center - wave[left]) / (wave[left+1] - wave[left])
    resolution = hdu["RESOLUTION"].data
    half = resolution.shape[1] // 2
    offsets = np.arange(-half, half + 1)
    bins = np.arange(left - half - 4, left + half + 6)
    widths = np.gradient(wave)
    def measurement(f):
        if f not in rowmap:
            return {"fiber": f, "unavailable": "missing_fibermap_row"}
        i = rowmap[f]
        profile = np.zeros(len(bins))
        for column, weight in [(left, 1-t), (left+1, t)]:
            response = resolution[i, :, column]
            norm = np.dot(response, widths[column+offsets])
            if not np.isfinite(norm) or norm <= 0:
                return {"fiber": f, "unavailable": "invalid_resolution"}
            profile[column + offsets - bins[0]] += weight * response / norm
        mask = hdu["MASK"].data[i, bins]
        flags = int(fib["FIBERSTATUS"][i])
        flux = hdu["FLUX"].data[i, bins]
        result = dict(fiber=int(f), targetid=int(fib["TARGETID"][i]),
                      objtype=str(fib["OBJTYPE"][i]).strip(), fiberstatus=flags,
                      window_mask_or=int(np.bitwise_or.reduce(mask)),
                      median_flux_density17=float(np.nanmedian(flux)) if np.isfinite(flux).any() else None)
        fit = forced_line(wave[bins], flux, hdu["IVAR"].data[i,bins], mask, profile)
        if flags or fit is None:
            result["unavailable"] = "flagged_fiber_or_unusable_local_profile"
        else:
            result.update(fit)
            result["snr"] = fit["flux"] / fit["sigma"]
        return result

    # Compute flux units using the cframe BUNIT and measured pixel widths.
    sky_rows, sky_distances, skipped = [], [], []
    for f, sep in zip(controls["sky_fibers"], controls["sky_separation_arcsec"]):
        row = measurement(f)
        if "sigma" not in row:
            skipped.append(dict(separation_arcsec=sep, **row))
            continue
        sky_rows.append(row)
        sky_distances.append(sep)
        if len(sky_rows) == 3:
            break
    output = dict(targetid=targetid, fiber=fiber, expid=expid, arm=arm,
                  barycentric_A=center, stellar_rest_A=hit["rest_A"],
                  lsf_fwhm_A=hit["rest_fwhm_A"]*center/hit["rest_A"],
                  topo_to_bary=float(hdu[0].header["HELIOCOR"]),
                  residual_flux17=hit["flux17"], residual_snr=hit["snr"],
                  bunit=hdu[0].header["BUNIT"], trace_fibers=controls["trace_fibers"],
                  sky_fibers=[s["fiber"] for s in sky_rows], sky_separation_arcsec=sky_distances,
                  nearer_unusable_sky=skipped)
    output["topocentric_A"] = center/output["topo_to_bary"]
    output["target_raw"] = measurement(fiber)
    output["sky_controls"] = sky_rows
    output["trace_controls"] = [measurement(f) for f in controls["trace_fibers"]]
    return output


def do_frame(job):
    (night, expid, camera), entries = job
    name = f"cframe-{camera}-{expid:08d}.fits.gz"
    prefix = f"{BASE}{night}/{expid:08d}/"
    checksum_name = f"redux_iron_exposures_{night}_{expid:08d}.sha256sum"
    checksum_bytes = download(prefix + checksum_name)
    checks = {row.split()[1].lstrip("*"): row.split()[0] for row in checksum_bytes.decode().splitlines() if row.strip()}
    compressed = download(prefix + name)
    digest = hashlib.sha256(compressed).hexdigest()
    if digest != checks[name]:
        raise ValueError(f"cframe checksum mismatch: {name}")
    manifest = dict(url=prefix+name, bytes=len(compressed), sha256=digest,
                    checksum_url=prefix+checksum_name, checksum_sha256=hashlib.sha256(checksum_bytes).hexdigest())
    with fits.open(io.BytesIO(gzip.decompress(compressed)), memmap=False) as hdu:
        data = hdu["FIBERMAP"].data
        out = [extract(hdu, data, entry) for entry in entries]
    return manifest, out


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("scans", type=Path)
    p.add_argument("output", type=Path)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--reuse", action="store_true",
                   help="Reuse stored measures when the complete target-centroid set agrees.")
    args = p.parse_args()
    args.output.mkdir(exist_ok=True, parents=True)
    jobs = defaultdict(list)
    for scan in sorted(args.scans.glob("scan-*.jsonl")):
        for line in scan.read_text().splitlines():
            s = json.loads(line)
            m = s["meta"]
            epochs = {e["expid"]: e for e in m["exposures"]}
            for arm, result in s["arms"].items():
                for candidate in result["candidates"]:
                    for hit in candidate["hits"]:
                        e = epochs[hit["expid"]]
                        jobs[e["night"], e["expid"], arm+str(e["petal_loc"])].append((m["targetid"], e, hit, arm))
    todo = []
    for key, group in jobs.items():
        _, expid, camera = key
        file = args.output / f"cframe-{camera}-{expid:08d}.json"
        if args.reuse and file.exists():
            saved = json.loads(file.read_text())["measurements"]
            current = {(t, h["observed_A"]) for t, _, h, _ in group}
            old = {(s["targetid"], s["barycentric_A"]) for s in saved if "nearer_unusable_sky" in s}
            if current == old:
                continue
        todo.append((key, group))
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        for manifest, rows in pool.map(do_frame, todo):
            data = dict(source=manifest, measurements=rows)
            file = args.output / (manifest["url"].split("/")[-1].split(".")[0]+".json")
            tmp = file.with_suffix(".partial")
            tmp.write_text(json.dumps(data, allow_nan=False, indent=2)+"\n")
            tmp.replace(file)
            print(file.name, len(rows), flush=True)


if __name__ == "__main__":
    main()
