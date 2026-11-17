"""Download and join DESI DR1 single-exposure spectra and matched RVSpecFit models.

Use the HEALPix spectra as an exposure inventory within one (survey, program).
A cross-program inventory is required before interpreting the search as all DR1.
DESI TARGETID and EXPID form a composite key. Every exposure must have a
successful, unflagged RVSpecFit stellar model and a usable fibre assignment.
"""
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from contextlib import ExitStack
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess
import tempfile

import fitsio
import numpy as np

BASE = "https://desidata.s3.amazonaws.com/"
ARMS = ("b", "r", "z")


def product_groups(survey, program, hp):
    tag = f"{survey}-{program}-{hp}"
    common = f"healpix/{survey}/{program}/{hp // 100}/{hp}"
    redux = f"dr1/spectro/redux/iron/{common}"
    vac = f"dr1/vac/dr1/mws/iron/v1.0/rv_output/240521/{common}"
    return [(redux, [f"spectra-{tag}.fits.gz", f"redrock-{tag}.fits", f"hpixexp-{tag}.csv"],
             f"redux_iron_{common.replace('/', '_')}.sha256sum"),
            (vac, [f"rvtab_spectra-{tag}.fits", f"rvmod_spectra-{tag}.fits"],
             f"{vac.replace('/', '_')}.sha256sum")]


def digest(path):
    with open(path, "rb") as f:
        return sha256_file(f)


def sha256_file(f):
    h = sha256()
    for block in iter(lambda: f.read(2**20), b""):
        h.update(block)
    return h.hexdigest()


def obtain(url, path, expected=None):
    if path.is_file() and (not expected or digest(path) == expected):
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as f:
        tmp = Path(f.name)
    try:
        subprocess.run(["curl", "--fail", "--silent", "--show-error", "--location",
                        "--retry", "3", "--connect-timeout", "20", "--max-time", "900",
                        "--output", str(tmp), url], check=True)
        if expected and digest(tmp) != expected:
            raise ValueError(f"DESI checksum mismatch: {url}")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def acquire(cache, survey, program, hp):
    """Download the five scientific inputs and verify DESI's published SHA-256."""
    cache = Path(cache)
    jobs, metadata = [], []
    for folder, names, checksum in product_groups(survey, program, hp):
        url = BASE + folder + "/"
        obtain(url + checksum, cache / checksum)
        checks = {line.split()[1].lstrip("*"): line.split()[0]
                  for line in (cache / checksum).read_text().splitlines() if line.strip()}
        for name in names:
            jobs.append((url + name, cache / name, checks[name]))
        metadata.append({"url": url + checksum, "sha256": digest(cache / checksum),
                         "bytes": (cache / checksum).stat().st_size})
    with ThreadPoolExecutor(max_workers=4) as pool:
        list(pool.map(lambda args: obtain(*args), jobs))
    metadata += [{"url": url, "sha256": check, "bytes": path.stat().st_size}
                 for url, path, check in jobs]
    return sorted(metadata, key=lambda item: item["url"])


def keys(table):
    return [(int(t), int(e)) for t,e in zip(table["TARGETID"], table["EXPID"])]


def valid_rv(row):
    return (row['SUCCESS'] and row['RVS_WARN'] == 0
            and np.isfinite(row['VRAD']) and np.isfinite(row['VRAD_ERR'])
            and row['VRAD_ERR'] >= 0)


def _targets_from_fits(cache, survey, program, hp, audit=None):
    """Yield target metadata and {arm: [Exposure, ...]} for all eligible stars.

    The wavelength arrays must match exactly; the model is an observed-frame
    output of the individual-exposure stellar fit. Do not convolve it again.
    """
    from laser_search import Exposure
    cache = Path(cache)
    tag = f"{survey}-{program}-{hp}"
    count = audit if audit is not None else Counter()
    with ExitStack() as stack:
        def opened(name):
            return stack.enter_context(fitsio.FITS(str(cache / name)))
        spectra = opened(f"spectra-{tag}.fits.gz")
        redrock = opened(f"redrock-{tag}.fits")["REDSHIFTS"].read()
        rv = opened(f"rvtab_spectra-{tag}.fits")["RVTAB"].read()
        model = opened(f"rvmod_spectra-{tag}.fits")
        fib = spectra["FIBERMAP"].read()
        spec_keys, rv_keys = keys(fib), keys(rv)
        mod_keys = keys(model["FIBERMAP"].read())
        if len(spec_keys) != len(set(spec_keys)) or len(rv_keys) != len(set(rv_keys)) or len(mod_keys) != len(set(mod_keys)):
            raise ValueError("duplicate TARGETID/EXPID in joined products")
        rv_rows = {key: i for i, key in enumerate(rv_keys)}
        model_rows = {key: i for i, key in enumerate(mod_keys)}
        stars = {int(x['TARGETID']) for x in redrock if x['SPECTYPE'] == 'STAR' and x['ZWARN'] == 0}
        waves = {}
        for arm in ARMS:
            ext = arm.upper() + '_WAVELENGTH'
            w = spectra[ext].read()
            if not np.array_equal(w, model[ext].read()):
                raise ValueError(f'{arm}-arm RVSpecFit and DESI grids differ')
            waves[arm] = w
        groups = defaultdict(list)
        for i, (target, expid) in enumerate(spec_keys):
            groups[target].append((expid, i))
        count['all_spectra_fibermap_rows'] = len(fib)
        count['all_unique_targets'] = len(groups)
        for target, rows in sorted(groups.items()):
            if target not in stars:
                count['nonstar_or_zwarn'] += 1
                continue
            if len(rows) < 2:
                count['star_single_exposure'] += 1
                continue
            rows.sort()
            joined = [(expid, spec_i, rv_rows.get((target, expid)), model_rows.get((target, expid)))
                      for expid, spec_i in rows]
            if any(ri is None or mi is None or fib[si]['FIBERSTATUS'] != 0 or not valid_rv(rv[ri])
                   for _,si,ri,mi in joined):
                count['repeat_star_incomplete_fits_or_fiber'] += 1
                continue
            x = fib[rows[0][1]]
            meta = {k.lower(): x[k].item() for k in ('TARGETID','TARGET_RA','TARGET_DEC','REF_ID','REF_CAT','PARALLAX')}
            meta.update(survey=survey, program=program, healpix=hp, exposures=[])
            by_arm = {arm: [] for arm in ARMS}
            for expid, i, j, m in joined:
                meta['exposures'].append({k.lower(): fib[i][k].item() for k in (
                    'EXPID','NIGHT','TILEID','PETAL_LOC','FIBER','MJD','EXPTIME')})
                for arm in ARMS:
                    def row(name, file=spectra, index=i):
                        hdu = file[arm.upper()+'_'+name]
                        ndim = len(hdu.get_info()['dims'])
                        return hdu[(slice(index, index + 1),) + (slice(None),) * (ndim - 1)][0]
                    by_arm[arm].append(Exposure(expid, waves[arm], row('FLUX'), row('IVAR'),
                                                row('MASK'), row('RESOLUTION'), row('MODEL', model, m),
                                                float(rv[j]['VRAD']), float(rv[j]['VRAD_ERR'])))
            count['eligible_stars'] += 1
            count['eligible_exposures'] += len(joined)
            yield meta, by_arm


def targets(cache, survey, program, hp, audit=None):
    """Read original FITS inputs, or locally packed eligible rows when available."""
    from laser_search import Exposure
    packed = Path(cache) / f"packed-{survey}-{program}-{hp}"
    index = packed / "index.json"
    if not index.is_file():
        yield from _targets_from_fits(cache, survey, program, hp, audit)
        return
    info = json.loads(index.read_text())
    if audit is not None:
        audit.update(info["audit"])
    for item in info["targets"]:
        path = packed / item["file"]
        if digest(path) != item["sha256"]:
            raise ValueError(f"packed row checksum mismatch: {path}")
        with np.load(path, allow_pickle=False) as data:
            by_arm = {}
            for arm in ARMS:
                by_arm[arm] = []
                for i, visit in enumerate(item["meta"]["exposures"]):
                    prefix = f"{arm}_{i}_"
                    by_arm[arm].append(Exposure(visit["expid"], *[
                        data[prefix + key] for key in ("wave", "flux", "ivar", "mask",
                        "resolution", "template")], *data[prefix+"velocity"]))
            yield item["meta"], by_arm


def pack(cache, survey, program, hp):
    """Losslessly cache all selected FITS rows with source and output hashes."""
    import shutil
    cache = Path(cache)
    name = f"packed-{survey}-{program}-{hp}"
    if (cache / name / "index.json").exists():
        return
    inputs = acquire(cache, survey, program, hp)
    staging = Path(tempfile.mkdtemp(dir=cache, prefix=name))
    index, audit = [], Counter()
    try:
        for meta, arms in _targets_from_fits(cache, survey, program, hp, audit):
            arrays = {}
            for arm, values in arms.items():
                for i, e in enumerate(values):
                    prefix = f"{arm}_{i}_"
                    for key in ("wave","flux","ivar","mask","resolution","template"):
                        arrays[prefix+key] = getattr(e, key)
                    arrays[prefix+"velocity"] = [e.velocity_kms, e.velocity_error_kms]
            path = staging / f"{meta['targetid']}.npz"
            np.savez_compressed(path, **arrays)
            index.append({"meta":meta, "file":path.name, "sha256":digest(path)})
        (staging / "index.json").write_text(json.dumps(
            {"inputs":inputs, "audit":audit, "targets":index}, sort_keys=True))
        os.replace(staging, cache / name)
    finally:
        if staging.exists():
            shutil.rmtree(staging)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pack', action='store_true', help='lossless eligible-row cache')
    parser.add_argument('cache', type=Path)
    parser.add_argument('survey')
    parser.add_argument('program')
    parser.add_argument('healpix', type=int)
    args = parser.parse_args()
    manifest = acquire(args.cache, args.survey, args.program, args.healpix)
    if args.pack:
        pack(args.cache, args.survey, args.program, args.healpix)
    print(json.dumps(manifest, indent=2))
