"""Re-measure the two anomalous ZTF catalog points of Hephaistos D on science cutouts."""
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import sys
from urllib.parse import urlencode
from urllib.request import urlopen
import zipfile

from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "wise"))
from fetch import position, targets
DAYS = (58735, 58745, 58748, 58749)
RADII = (4, 6, 8)
PIXELS = ROOT / "d_science_cutouts.zip"


def fetch(row, product, center):
    stamp = row["filefracday"]
    ccd, qid = int(row["ccdid"], 0), int(row["qid"], 0)
    name = f"ztf_{stamp}_{int(row['field']):06d}_zr_c{ccd:02d}_o_q{qid}_{product}.fits"
    root = "https://irsa.ipac.caltech.edu/ibe/data/ztf/products/sci/"
    url = f"{root}{stamp[:4]}/{stamp[4:8]}/{stamp[8:]}/{name}?" + urlencode(
        dict(center=f"{center[0]:.7f},{center[1]:.7f}", size="240arcsec", gzip="false"))
    with urlopen(url, timeout=150) as response:
        data = response.read()
    with fits.open(io.BytesIO(data)) as hdu:
        if hdu[0].data is None or hdu[0].data.ndim != 2:
            raise ValueError(f"Invalid image at {url}")
    return name, data, dict(expid=row["expid"], product=product, url=url, name=name, bytes=len(data),
                            sha256=hashlib.sha256(data).hexdigest())


def cutouts(exposures, center):
    receipt_path = ROOT / "d_science_receipts.json"
    if not PIXELS.exists():
        tasks = [(r, product, center) for r in exposures for product in ("sciimg", "mskimg")]
        with ThreadPoolExecutor(4) as pool:
            results = list(pool.map(lambda t: fetch(*t), tasks))
        with zipfile.ZipFile(PIXELS, "w", zipfile.ZIP_DEFLATED) as archive:
            for name, data, _ in results:
                archive.writestr(name, data)
        receipt_path.write_text(json.dumps([receipt for _, _, receipt in results], indent=2) + "\n")
    out = {}
    meta = json.loads(receipt_path.read_text())
    with zipfile.ZipFile(PIXELS) as archive:
        for entry in meta:
            data = archive.read(entry["name"])
            if hashlib.sha256(data).hexdigest() != entry["sha256"]:
                raise ValueError(entry["name"])
            with fits.open(io.BytesIO(data)) as hdus:
                out[entry["expid"], entry["product"]] = (hdus[0].data.copy(), hdus[0].header.copy())
    return out


def controls(image, wcs, target_radec):
    safe = np.where(np.isfinite(image), image, np.nanmedian(image))
    smooth = gaussian_filter(safe, 1.)
    mid = np.median(smooth)
    scatter = 1.4826 * np.median(np.abs(smooth - mid))
    y, x = np.where((smooth == maximum_filter(smooth, size=7)) & (smooth > mid + 15 * scatter))
    # Reference-image peaks resolve isolated stars away from cutout edges and D.
    if len(x) < 2:
        raise ValueError("Too few comparison peaks")
    points = np.column_stack([x, y])
    nearest = cKDTree(points).query(points, k=2)[0][:, 1]
    tx, ty = wcs.all_world2pix(*target_radec, 0)
    acceptable = ((nearest > 20) & (np.hypot(x-tx, y-ty) > 20) & (x>18) & (y>18) &
                  (x<image.shape[1]-18) & (y<image.shape[0]-18) & (safe[y, x] < 40000))
    return np.asarray(wcs.all_pix2world(points[acceptable], 0))


def aperture(image, mask, x, y, radius):
    xx, yy = int(np.round(x)), int(np.round(y))
    patch = image[yy-15:yy+16, xx-15:xx+16]
    mm = mask[yy-15:yy+16, xx-15:xx+16].astype(np.int64)
    if patch.shape != (31,31):
        return np.nan, -1
    sy, sx = np.indices(patch.shape)
    rr = np.hypot(xx-15+sx-x, yy-15+sy-y)
    sky = np.nanmedian(patch[(rr>10) & (rr<=15)])
    inside = rr<=radius
    return np.sum(patch[inside]-sky), np.bitwise_or.reduce(mm[rr<2], initial=0)


def main():
    target = next(t for t in targets() if t["label"]=="D")
    center = position(target, 2000+(58745-51544.5)/365.25)
    with zipfile.ZipFile(ROOT / "lightcurve_inputs.zip") as archive:
        body = gzip.decompress(archive.read("cache/D.csv.gz"))
    exposures = [r for r in csv.DictReader(io.StringIO(body.decode())) if r["field"]=="1544"
                 and int(float(r["mjd"])) in DAYS and int(r["catflags"])==0]
    assert len(exposures) == 8
    frames = cutouts(exposures, center)
    img, header = frames[exposures[0]["expid"], "sciimg"]
    sky = np.vstack([center, controls(img, WCS(header), center)])
    base = np.array([int(float(r["mjd"])) != 58745 for r in exposures])
    summaries, curves = [], []
    for radius in RADII:
        flux = np.zeros((len(exposures), len(sky)))
        flags = np.zeros_like(flux, dtype=int)
        for i, r in enumerate(exposures):
            data, header = frames[r["expid"], "sciimg"]
            mask, _ = frames[r["expid"], "mskimg"]
            for j, (x,y) in enumerate(WCS(header).all_world2pix(sky, 0)):
                flux[i,j], flags[i,j] = aperture(data, mask, x, y, radius)
        use = np.all(np.isfinite(flux) & (flux>0) & (flags==0), axis=0)
        use[0] = False
        if use.sum() < 3:
            raise ValueError(f"Only {use.sum()} control stars")
        normalized = flux / np.median(flux[base], axis=0)
        transparency = np.median(normalized[:,use], axis=1)
        relative = -2.5 * np.log10(normalized / transparency[:,None])
        relative -= np.median(relative[base], axis=0)
        for i,r in enumerate(exposures):
            residuals = relative[i,use]
            summary = dict(expid=r["expid"], mjd=r["mjd"], radius_px=radius, n_controls=int(use.sum()),
                           catalog_mag=r["mag"], catalog_chi=r["chi"], target_dmag=round(relative[i,0],5),
                           control_mad_mag=round(1.4826*np.median(np.abs(residuals-np.median(residuals))),5),
                           control_p95_abs_mag=round(np.percentile(abs(residuals),95),5), target_core_mask=flags[i,0])
            summaries.append(summary)
            for j in np.r_[0,np.flatnonzero(use)]:
                curves.append(dict(expid=r["expid"],radius_px=radius,object_index=j,ra=sky[j,0],dec=sky[j,1],
                                   flux=flux[i,j],core_mask=flags[i,j],relative_mag=relative[i,j]))
    for name, rows in [("d_aperture_summary.csv",summaries),("d_aperture_controls.csv",curves)]:
        with (ROOT/name).open("w",newline="") as stream:
            writer=csv.DictWriter(stream,rows[0].keys(),lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    for row in summaries:
        if int(float(row["mjd"]))==58745:
            print(row,flush=True)


if __name__ == "__main__":
    main()
