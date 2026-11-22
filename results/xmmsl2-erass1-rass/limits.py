"""Measure position-specific upper limits on observed 0.2–2 keV point-source flux.

Run only on the Actions runner, after catalogs.py. Rejected measurements are
retained, including off-footprint targets, invalid exposure and archive errors.
"""
from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import hashlib
import json
import math
import os
import re
import threading
import time
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import pandas as pd
import requests
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.stats import gamma

SKYVIEW = "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl"
CONFIDENCE = 0.99865
RASS_EEF = 0.5
RASS_APERTURE = 180.0  # arcsec, point-source assumption
THREAD = threading.local()
INPUTS = {}
MAP_HEADERS = {}
INPUT_LOCK = threading.Lock()
COLUMNS = [
    "detection_uid", "name", "ra", "dec", "mission", "status", "limit_flux_cgs",
    "confidence", "method", "counts", "exp_s", "eef", "ecf", "aperture_arcsec",
    "map_ref", "reason", "model_max_relative", "native_ul_b", "bkg_counts", "flag_pos", "field",
]


def session():
    if not hasattr(THREAD, "session"):
        THREAD.session = requests.Session()
        THREAD.session.headers["User-Agent"] = "Open-Space XMMSL2 transient study/1.0"
    return THREAD.session


def fetch(url, work, *, timeout=40, maxbytes=700_000_000, payload=None):
    """Cache exact immutable input bytes locally; record URLs and hashes."""
    request_key = url if payload is None else url + " " + json.dumps(payload, sort_keys=True)
    token = hashlib.sha256(request_key.encode()).hexdigest()
    dest = work / token
    if not dest.exists():
        for attempt in range(3):
            tmp = dest.with_suffix(f".tmp.{threading.get_ident()}")
            try:
                method = "GET" if payload is None else "POST"
                with session().request(method, url, json=payload, timeout=(12, timeout), stream=True) as r:
                    r.raise_for_status()
                    with tmp.open("wb") as f:
                        for block in r.iter_content(262144):
                            if block:
                                f.write(block)
                            if f.tell() > maxbytes:
                                raise ValueError("map exceeds configured single-file limit")
                os.replace(tmp, dest)
                break
            except (requests.RequestException, OSError) as exc:
                tmp.unlink(missing_ok=True)
                if ((isinstance(exc, requests.HTTPError) and
                     exc.response.status_code not in (408, 429, 500, 502, 503, 504))
                        or attempt == 2):
                    raise
                time.sleep(2 ** attempt)
            except Exception:
                tmp.unlink(missing_ok=True)
                raise
    with INPUT_LOCK:
        if token not in INPUTS:
            with dest.open("rb") as f:
                ingest = hashlib.file_digest(f, "sha256").hexdigest()
            INPUTS[token] = (url, dest.stat().st_size, ingest)
    return dest, token


def skyview_url(survey, ra, dec):
    # Return the native broad-band integer counts via nearest-neighbour sampling.
    # Each 20" output pixel samples one native 45" pixel. Exposure uses exactly
    # the same WCS. No intensity images, catalogue deblending or background fits.
    params = dict(Survey=survey, Position=f"{ra:.8f},{dec:.8f}", Coordinates="J2000",
                  Projection="Tan", Pixels="35", Size=f"{35 * 20 / 3600:.10f}",
                  Sampler="NN", Return="FITS")
    return SKYVIEW + "?" + urlencode(params)


def get_patch(survey, ra, dec, work):
    path, ref = fetch(skyview_url(survey, ra, dec), work, maxbytes=2_000_000)
    with fits.open(path, memmap=False) as hdus:
        hdu = next((h for h in hdus if h.data is not None and h.data.ndim == 2), None)
        if hdu is None or hdu.data.shape != (35, 35):
            raise ValueError("response is not the requested SkyView image")
        wcs = WCS(hdu.header).celestial
        if not np.allclose(proj_plane_pixel_scales(wcs) * 3600., 20., atol=0.02):
            raise ValueError("unexpected SkyView pixel scale")
        yy, xx = np.indices(hdu.data.shape)
        coords = wcs.pixel_to_world(xx, yy)
        sep = coords.separation(SkyCoord(ra, dec, unit="deg")).arcsec
        if np.nanmin(sep) > 15:
            raise ValueError("position offset from requested centre")
        # Include a 65" edge around the 3' physical aperture. A 45" square
        # native pixel subtends at most ~64" across. With grid spacing 20"
        # at least one output centre samples every included native pixel.
        keep = sep < RASS_APERTURE + 65.
        return np.array(hdu.data[keep], dtype=float), ref


def rass_cutout(row, work, model):
    """An upper bound allowing all counts in the aperture to be source counts.

    Counting repeated nearest-neighbour samples increases the count total. This
    conservatively bounds the Poisson upper limit constructed from the actual
    aperture counts. No subtraction of background or neighbouring sources.
    """
    ecf = model.get("rass_ecf")
    if not ecf or not np.isfinite(ecf) or ecf <= 0:
        return {"status": "unavailable", "reason": "no calibrated ROSAT count-to-flux conversion"}
    exp, eref = get_patch("RASS Exposure", row.ra, row.dec, work)
    if not np.all(np.isfinite(exp)) or np.min(exp) <= 0:
        return {"status": "no_exposure_or_masked", "map_ref": eref}
    cnt, cref = get_patch("RASS-Cnt Broad", row.ra, row.dec, work)
    if not np.all(np.isfinite(cnt) & (cnt >= 0) & np.isclose(cnt, np.rint(cnt), atol=1e-4, rtol=0)):
        raise ValueError("invalid or noninteger count pixels")
    total, exposure = int(np.rint(cnt).sum()), float(np.min(exp))
    high = gamma.ppf(CONFIDENCE, total + 1.0)
    return dict(status="ok", limit_flux_cgs=ecf * high / (exposure * RASS_EEF),
                confidence=CONFIDENCE, counts=total, exp_s=exposure, eef=RASS_EEF,
                ecf=ecf, aperture_arcsec=RASS_APERTURE, map_ref=f"{eref} {cref}",
                method="RASS native broad maps via SkyView; oversampled counts upper bound")


def write_inputs(out):
    with (out / "limits_inputs.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["map_ref", "url", "bytes", "sha256"])
        w.writerows((k, *v) for k, v in sorted(INPUTS.items()))


def note_map_header(url, hdus):
    # First observed header per archive product, including units and WCS for QA.
    suffix = url.rsplit("/",1)[-1].split("_")[-1]
    with INPUT_LOCK:
        if suffix not in MAP_HEADERS:
            keys = ["XTENSION","NAXIS","NAXIS1","NAXIS2","TELESCOP","INSTRUME","CONTENT","EXTNAME","HDUCLAS1","HDUCLAS2","HDUCLAS3","BUNIT","CTYPE1","CTYPE2","CUNIT1","CUNIT2","CDELT1","CDELT2","EQUINOX"]
            MAP_HEADERS[suffix] = {"url":url, "headers":[{k:h.header[k] for k in keys if k in h.header} for h in hdus]}


def image_hdu(hdus):
    for h in hdus:
        if h.data is not None and h.data.ndim == 2 and WCS(h.header).has_celestial:
            return h
    raise ValueError("no two-dimensional celestial FITS image")


def physical_native_limit(cnt, exp, wcs, coord, ecf):
    """Upper bound on point-source counts using the native integer image."""
    x, y = wcs.world_to_pixel(coord)
    scales = proj_plane_pixel_scales(wcs) * 3600.
    # RASS fields are small gnomonic projections. Include an entire edge pixel
    # for every pixel touched by the requested aperture.
    radius = RASS_APERTURE + math.hypot(*scales)
    half = int(np.ceil(radius / min(scales))) + 2
    i, j = int(np.rint(x)), int(np.rint(y))
    if (i - half < 0 or j - half < 0 or
            i + half >= cnt.shape[1] or j + half >= cnt.shape[0]):
        # A bounding *square* can cross a tile edge while the actual circular
        # aperture is still fully covered. Check its perimeter in the native WCS.
        # The 2" radial margin exceeds the 0.18" chord sagitta at 5° spacing.
        from astropy import units as u
        ring = coord.directional_offset_by(np.arange(0,360,5)*u.deg, (RASS_APERTURE+2)*u.arcsec)
        rx, ry = wcs.world_to_pixel(ring)
        if not np.all((rx >= -.5) & (ry >= -.5) & (rx <= cnt.shape[1]-.5) & (ry <= cnt.shape[0]-.5)):
            return {"status":"edge"}
    ys = slice(max(0,j-half), min(cnt.shape[0],j+half+1))
    xs = slice(max(0,i-half), min(cnt.shape[1],i+half+1))
    yy, xx = np.mgrid[ys, xs]
    aperture = wcs.pixel_to_world(xx, yy).separation(coord).arcsec <= radius
    expvals = np.asarray(exp[ys, xs][aperture], dtype=float)
    values = np.asarray(cnt[ys, xs][aperture], dtype=float)
    if not np.all(np.isfinite(expvals)) or np.min(expvals) <= 0:
        return {"status": "no_exposure_or_masked"}
    if not np.all(np.isfinite(values) & (values >= 0) & np.isclose(values, np.rint(values), atol=1e-4, rtol=0)):
        return {"status": "invalid_count_pixels"}
    n, t = int(np.rint(values).sum()), float(np.min(expvals))
    return dict(status="ok", counts=n, exp_s=t, eef=RASS_EEF, ecf=ecf,
                confidence=CONFIDENCE, aperture_arcsec=RASS_APERTURE,
                limit_flux_cgs=gamma.ppf(CONFIDENCE, n+1.0) * ecf / (t * RASS_EEF),
                method="RASS native broad count and exposure maps; unsubtracted aperture")


def list_urls(url, work):
    from html.parser import HTMLParser
    from urllib.parse import urljoin
    class Links(HTMLParser):
        links = None

        def handle_starttag(self, tag, attrs):
            if tag == "a":
                href = dict(attrs).get("href", "")
                if href and not href.startswith(("?", "#", "../", "/icons/")):
                    self.links.append(urljoin(url, href))
    parser = Links()
    parser.links = []
    path, _ = fetch(url, work, maxbytes=12_000_000)
    parser.feed(path.read_text(errors="replace"))
    return list(dict.fromkeys(parser.links))


def likely_maps(urls):
    counts = sorted(x for x in urls if re.search(r'[_\-](?:im1|imec|image1|counts|im_broad)(?:\.fits?|\.fit\.)(?:\.gz|\.Z)?$', x, re.I))
    exposure = sorted(x for x in urls if re.search(r'[_\-](?:mex|ime|ex1|exp|exposure|exposure_broad)(?:\.fits?|\.fit\.)(?:\.gz|\.Z)?$', x, re.I))
    if len(counts) == 1 and exposure:
        # MEX and IME are different archive products. Inspect their FITS units
        # and celestial grids before treating either as an exposure measurement.
        exposure.sort(key=lambda x: ("_mex." not in x.lower(),x))
        return counts[0], tuple(exposure)
    return None


def native_rass(rows, work, model, workers):
    """Search released native all-sky fields and choose by centre distance.

    Decisions among overlapping fields depend on geometry and coverage. An
    independent fallback can measure unresolved directions via SkyView.
    """
    from scipy.spatial import cKDTree
    base = "https://heasarc.gsfc.nasa.gov/FTP/rosat/data/pspc/processed_data/rass/release/"
    field = "rs930101n00"
    for doc in ("public_contents","public_contents_mpe"):
        try:
            url = f"{base}{field}/{field}.{doc}"
            path, _ = fetch(url, work, maxbytes=64000)
            MAP_HEADERS[doc] = {"url":url, "text":path.read_text(errors="replace")[:9000]}
        except (OSError, ValueError) as exc:
            MAP_HEADERS[doc] = {"error":str(exc)[:240]}
    directories = native_inventory(base, work)
    if not directories:
        return {}, ["RASS release index contained no identified field maps"]
    coords = SkyCoord([r.ra for r in rows], [r.dec for r in rows], unit="deg")
    tree = cKDTree(coords.cartesian.xyz.value.T)
    best, failure_counts = {}, {}
    def inspect(directory):
        maps = directory if isinstance(directory, tuple) else likely_maps(list_urls(directory, work))
        if maps is None:
            return [], "missing map pair"
        countpath, cref = fetch(maps[0], work)
        with open_map(countpath) as cc:
            note_map_header(maps[0],cc)
            ch = image_hdu(cc)
            wcs, shape = WCS(ch.header).celestial, ch.shape
            xmid, ymid = (shape[1] - 1)/2, (shape[0] - 1)/2
            center = wcs.pixel_to_world(xmid, ymid)
            corners = wcs.pixel_to_world(np.array([0,shape[1]-1,shape[1]-1,0]), np.array([0,0,shape[0]-1,shape[0]-1]))
            sep = center.separation(corners).max().radian
            if not 0 < sep < math.radians(15):
                raise ValueError("unexpected RASS field angular dimensions")
            nearby = tree.query_ball_point(center.cartesian.xyz.value, 2*math.sin(sep/2))
            if not nearby:
                return [], None
            from astropy import units as u
            problems = []
            for candidate in maps[1]:
                try:
                    exppath, eref = fetch(candidate, work)
                    with open_map(exppath) as ee:
                        note_map_header(candidate,ee)
                        eh = image_hdu(ee)
                        ewcs = WCS(eh.header).celestial
                        if ch.shape != eh.shape or not np.allclose(
                            ewcs.world_to_pixel(corners),
                            wcs.world_to_pixel(corners), atol=0.03, rtol=0):
                            raise ValueError("exposure and count WCS grids differ")
                        unit = str(eh.header.get("BUNIT", "")).strip().lower()
                        labels = " ".join(str(eh.header.get(k,"")) for k in
                                          ["CONTENT","EXTNAME","HDUCLAS1","HDUCLAS2","HDUCLAS3"]).lower()
                        if unit in ("s","sec","secs","seconds","second"):
                            seconds = 1.
                        elif unit:
                            seconds = float(u.Unit(unit).to(u.s))
                        elif "expos" in labels and candidate.lower().split("_")[-1].startswith(("ime.","mex.")):
                            # The archive's survey exposure-map product uses seconds.
                            seconds = 1.
                        else:
                            raise ValueError(f"missing exposure units; FITS classification {labels!r}")
                        exposure = np.asarray(eh.data, dtype=float)*seconds
                        result = []
                        for k in nearby:
                            d = physical_native_limit(ch.data, exposure, wcs, coords[k], model["rass_ecf"])
                            d["map_ref"] = f"{eref} {cref}"
                            rank = (d["status"] != "ok", center.separation(coords[k]).deg, str(directory))
                            result.append((k, rank, d))
                        return result, None
                except Exception as exc:
                    problems.append(candidate.rsplit("/",1)[-1].split("_")[-1]+": "+str(exc))
            raise ValueError("; ".join(problems))

    def safe(directory):
        try:
            return inspect(directory)
        except Exception as e:
            return [], f"{type(e).__name__}: {str(e)[:170]}"
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        for result, error in executor.map(safe, directories):
            if error:
                failure_counts[error] = failure_counts.get(error, 0) + 1
            for k, rank, d in result:
                if k not in best or rank < best[k][0]:
                    best[k] = (rank, d)
    return {rows[k].detection_uid: v[1] for k,v in best.items()}, failure_counts


def open_map(path):
    with path.open("rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x9d":
        import io
        import subprocess
        data = subprocess.run(["gzip", "-dc", str(path)], capture_output=True, check=True).stdout
        return fits.open(io.BytesIO(data), memmap=False)
    return fits.open(path, memmap=False)


def native_inventory(base, work):
    from collections import defaultdict
    from urllib.parse import quote, urlsplit
    from xml.etree import ElementTree
    prefix = "rosat/data/pspc/processed_data/rass/release/"
    bucket = "https://nasa-heasarc.s3.amazonaws.com/"
    urls = []
    try:
        params = {"list-type": "2", "prefix": prefix+"rs"}
        while True:
            path, _ = fetch(bucket + "?" + urlencode(params), work, timeout=35, maxbytes=8_000_000)
            tree = ElementTree.parse(path)
            keys = [n.text for n in tree.findall(".//{*}Key") if n.text]
            urls.extend(bucket + quote(key) for key in keys)
            token = tree.findtext(".//{*}NextContinuationToken")
            if not token:
                break
            if token == params.get("continuation-token"):
                raise ValueError("repeated S3 continuation token")
            params["continuation-token"] = token
    except (ValueError, OSError, ElementTree.ParseError):
        pass
    def group(links):
        groups = defaultdict(list)
        for url in links:
            basename = urlsplit(url).path.split("/")[-1]
            m = re.match(r"(rs\d{6}[a-z]\d{2})[_\-]", basename, re.I)
            if m:
                groups[m.group(1)].append(url)
        return [m for files in groups.values() if (m := likely_maps(files))]
    complete = group(urls)
    if complete:
        return complete
    # The RASS field table gives the release directory names.
    from astropy.table import Table
    params = {"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"votable","MAXREC":10000,
              "QUERY":"SELECT ror,index_id FROM rassmaster"}
    try:
        table_path, _ = fetch("https://heasarc.gsfc.nasa.gov/xamin/vo/tap/sync?" + urlencode(params), work)
        table = Table.read(table_path, format="votable")
        pairs = []
        for row in table:
            field = f"rs{str(row['ror']).strip()}{str(row['index_id']).strip().lower()}"
            if re.fullmatch(r"rs9[0-9]{5}[a-z][0-9]{2}",field):
                stem = base + field + "/" + field
                pairs.append((stem+"_im1.fits.Z", (stem+"_mex.fits",stem+"_ime.fits.Z")))
        if pairs:
            return list(dict.fromkeys(pairs))
    except (ValueError, OSError, KeyError):
        pass
    options = list_urls(base, work)
    complete = group(options)
    if complete:
        return complete
    return [s for s in options if s.startswith(base) and s.endswith("/") and s != base]


ERO_API = "https://erosita.mpe.mpg.de/erodat/upperlimit/service_multi"


def erass_limits(rows, work, model, out):
    results = {}
    report = out / "ero_ul_queries.jsonl.gz"
    with gzip.open(report, "wt", encoding="utf-8") as request_log:
        for start in range(0, len(rows), 500):
            batch = rows[start:start+500]
            payload = [{"ra":r.ra, "dec":r.dec, "band":"024", "dr_survey":"DR1_eRASS1"} for r in batch]
            try:
                path, ref = fetch(ERO_API, work, payload=payload, timeout=150, maxbytes=10_000_000)
                request_log.write(json.dumps({"map_ref": ref, "positions": payload})+"\n")
                reply = json.loads(path.read_bytes())
                if reply.get("error") not in (None, "None") or len(reply["limits"]) != len(batch):
                    raise ValueError("failed or truncated eRODat multiple-position result")
                for src, rec in zip(batch, reply["limits"]):
                    result = dict(map_ref=ref, field=rec.get("field", ""))
                    if (rec.get("band") != "024" or rec.get("dr_survey") != "DR1_eRASS1" or
                            abs(float(rec["ra"]) - src.ra) > 1e-5 or abs(float(rec["dec"]) - src.dec) > 1e-5):
                        raise ValueError("eRODat survey, band or coordinate validation failed")
                    if rec.get("de_sky") is False:
                        result["status"] = "outside_dr1"
                    elif rec.get("de_sky") is not True:
                        result["status"] = "invalid_coverage_flag"
                    elif rec.get("Flag_pos") != 0:
                        result.update(status="near_detected_source", flag_pos=rec.get("Flag_pos"))
                    elif all(np.isfinite(float(rec.get(k) or 0)) and float(rec.get(k) or 0)>0 for k in ["UL_B","Exposure"]):
                        result.update(status="ok", limit_flux_cgs=float(rec["UL_B"])*model["ero_band_factor"],
                                      confidence=0.9987, method="eRODat DR1_eRASS1 band 024 UL_B",
                                      model_max_relative=model["ero_relative"][1], native_ul_b=rec["UL_B"],
                                      counts=rec.get("Count"), bkg_counts=rec.get("Bkg_counts"),
                                      exp_s=rec["Exposure"], flag_pos=rec["Flag_pos"], eef=0.75)
                    else:
                        result["status"] = "invalid_upper_limit_or_exposure"
                    results[src.detection_uid] = result
            except (ValueError, KeyError, TypeError, OSError) as exc:
                for src in batch:
                    results[src.detection_uid] = dict(status="archive_error", reason=str(exc)[:220])
    return results


def measure_rass(rows, work, model, workers):
    if not rows:
        return {}, {}
    ecf = model.get("rass_ecf")
    if not ecf or not np.isfinite(ecf) or ecf <= 0:
        return {r.detection_uid: dict(status="no_response_calibration") for r in rows}, {}
    try:
        result, errors = native_rass(rows, work, model, workers)
    except Exception as exc:
        result, errors = {}, {"release_inventory":str(exc)[:240]}
    missed = [r for r in rows if r.detection_uid not in result or result[r.detection_uid]["status"]=="edge"]
    # The lower-flux detections remain eligible. Stop remote recovery if repeated
    # service failures show that this endpoint is unusable; preserve every row.
    missed.sort(key=lambda r: (-r.xmm_flux if np.isfinite(r.xmm_flux) else math.inf,r.detection_uid))
    baseline = []
    def get(r):
        try:
            return r.detection_uid, rass_cutout(r, work, model)
        except Exception as exc:
            return r.detection_uid, dict(status="archive_error", reason=str(exc)[:200])
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for start in range(0, len(missed), 2*workers):
            batch = missed[start:start+2*workers]
            for uid, d in pool.map(get, batch):
                result[uid] = d
                baseline.append(d["status"]=="archive_error")
            if len(baseline) >= 2*workers and all(baseline[-2*workers:]):
                for r in missed[start+len(batch):]:
                    result[r.detection_uid] = dict(status="archive_unavailable",reason="repeated SkyView failures; native maps unavailable")
                break
    for r in result.values():
        if r.get("status") == "ok":
            r["model_max_relative"] = model["rass_relative"][1]
    return result, errors


def main():
    arg = argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--work", required=True, type=Path)
    arg.add_argument("--out", required=True, type=Path)
    arg.add_argument("--workers", type=int, default=6)
    args = arg.parse_args()
    work = args.work / "map-cache"
    work.mkdir(parents=True, exist_ok=True)
    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    cat = pd.read_csv(out/"catalog_sources.csv.gz", low_memory=False)
    if cat.detection_uid.isna().any() or not cat.detection_uid.is_unique:
        raise ValueError("missing or duplicate slew detection key")
    model = json.loads((out/"spectral_model.json").read_text())
    all_rows = list(cat.itertuples(index=False))
    def missing(src, tag):
        return pd.isna(getattr(src,tag+"_name")) or not str(getattr(src,tag+"_name")).strip()
    def boolcol(src, col):
        return str(getattr(src,col,False)).lower().strip() in ("true","1","1.0","t")
    def needed(src, tag):
        return missing(src, tag) and not boolcol(src,tag+"_ambiguous") and np.isfinite(src.ra) and np.isfinite(src.dec)
    ero_rows = [r for r in all_rows if needed(r,"ero")]
    rass_rows = [r for r in all_rows if needed(r,"rass")]
    ero = erass_limits(ero_rows, work, model, out)
    rass, errors = measure_rass(rass_rows, work, model, args.workers)
    from collections import Counter
    counts = Counter()
    with (out/"limits.csv").open("w",newline="") as f:
        w = csv.DictWriter(f,fieldnames=COLUMNS)
        w.writeheader()
        for r in all_rows:
            for tag, data in (("ero",ero),("rass",rass)):
                d = data.get(r.detection_uid)
                if d is None:
                    status = "catalogued_or_ambiguous" if boolcol(r,tag+"_ambiguous") or not missing(r,tag) else "invalid_coordinates"
                    d = {"status":status}
                counts[(tag,d["status"])] += 1
                w.writerow(dict(d,detection_uid=r.detection_uid,name=r.name,ra=r.ra,dec=r.dec,mission=tag))
    write_inputs(out)
    (out/"limits_summary.json").write_text(json.dumps({
        "counts":{tag:{status:n for (m,status),n in counts.items() if m==tag} for tag in ("ero","rass")},
        "rass_native_map_errors": errors,
        "rass_native_example_headers": MAP_HEADERS,
        "model_dependence":"model_max_relative applies only to the catalog spectral model grid. The eROSITA grid varies the bandpass ratio at fixed released ECF. The RASS grid includes the folded PSPC ECF. ROSAT throughput assumes ≥50% encircled energy for a point source within 3 arcmin. The limits do not constrain arbitrary spectra, extended emission, shorter flares or all-sky simultaneity.",
        "confidence":"One-sided: eROSITA release UL_B 99.87%; ROSAT background-unsubtracted Poisson 99.865%.",
    },indent=2)+"\n")
    print("Map upper limits:",dict(counts),flush=True)


if __name__ == "__main__":
    main()
