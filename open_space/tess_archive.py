"""Plan two-sector samples from NASA/MAST SPOC inventories and cache verified FITS."""

from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
import time
from urllib.error import HTTPError, URLError
from urllib.request import urlopen


PRODUCT = re.compile(r"tess\d{13}-s(\d{4})-(\d{16})-\d{4}-s_lc\.fits")
INVENTORY = "https://archive.stsci.edu/missions/tess/download_scripts/sector/tesscurl_sector_{}_lc.sh"
MAST = "https://mast.stsci.edu/api/v0.1/Download/file/?uri=mast:TESS/product/"
CONTROL_TIC = 100100827
SEED = "open-space-tess-v1"
ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PLAN = ROOT / "data/tess-plan.json"


def atomic_bytes(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as staged:
        staged.write(data)
        tmp = Path(staged.name)
    try:
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def download(url, cap):
    for attempt in range(3):
        try:
            with closing(urlopen(url, timeout=75)) as response:
                raw = response.read(cap + 1)
            if len(raw) > cap:
                raise ValueError("Archive response exceeds size bound")
            return raw
        except HTTPError as exc:
            if exc.code not in {408, 429, 500, 502, 503, 504} or attempt == 2:
                raise
        except (URLError, TimeoutError):
            if attempt == 2:
                raise
        time.sleep(2 ** attempt)


def inventory(sector, cached_dir):
    if not 1 <= sector <= 9999:
        raise ValueError("Sector must be between 1 and 9999")
    path = cached_dir / f"s{sector}.sh"
    url = INVENTORY.format(sector)
    if not path.exists():
        atomic_bytes(path, download(url, 16_000_000))
    raw = path.read_bytes()
    products = {}
    # Treat the published cURL script as text. Never execute archive-supplied commands.
    for line in raw.decode().splitlines():
        tokens = line.split()
        if not tokens or tokens[0] != "curl" or "-o" not in tokens:
            continue
        name = tokens[tokens.index("-o") + 1]
        match = PRODUCT.fullmatch(name)
        if not match or int(match[1]) != sector:
            continue
        tic = int(match[2])
        if tic in products and name != products[tic]:
            raise ValueError(f"Multiple reductions for TIC {tic} in sector {sector}")
        products[tic] = name
    if not products:
        raise ValueError(f"No standard SPOC light curves in sector {sector}")
    return products, {"sector": sector, "url": url, "sha256": hashlib.sha256(raw).hexdigest(),
                      "products": len(products)}


def make_plan(search_sector, test_sector, limit, offset, seed, catalog_dir, control=False):
    if search_sector >= test_sector or limit < 1 or offset < 0:
        raise ValueError("Require increasing sectors, limit >= 1, and a nonnegative offset")
    sectors = sorted({search_sector, test_sector, *([2, 3] if control else [])})
    catalogs, metadata = {}, []
    for sector in sectors:
        catalogs[sector], meta = inventory(sector, catalog_dir)
        metadata.append(meta)
    shared = catalogs[search_sector].keys() & catalogs[test_sector].keys()
    ordered = sorted(shared - {CONTROL_TIC}, key=lambda tic: (
        hashlib.sha256(f"{seed}:{tic}".encode()).hexdigest(), tic))

    def target(tic, role, sectors):
        return {"tic_id": tic, "role": role, "products": [
            {"sector": s, "filename": catalogs[s][tic],
             "uri": "mast:TESS/product/" + catalogs[s][tic]} for s in sectors]}

    targets = [target(tic, "exploration", (search_sector, test_sector)) for tic in ordered[offset:offset+limit]]
    if control:
        targets.insert(0, target(CONTROL_TIC, "known_control", (2, 3)))
    if not targets:
        raise ValueError("The requested page contains no targets")
    return {"schema_version": 1, "track": "tess", "selection": {
        "search_sector": search_sector, "validation_sector": test_sector, "shared_tics": len(shared),
        "seed": seed, "offset": offset, "limit": limit,
        "order": "Ascending SHA-256 of UTF-8 'seed:TIC', excluding the known control"},
        "inventories": metadata, "targets": targets}


def cached_product(product, cache):
    from astropy.io import fits

    cache = Path(cache)
    name = product["filename"]
    match = PRODUCT.fullmatch(name)
    if not match or int(match[1]) != product["sector"] or product["uri"] != "mast:TESS/product/" + name:
        raise ValueError("Invalid TESS product identity")
    path = cache / name
    if path.exists():
        raw = path.read_bytes()
    else:
        # The STScI public mirror names the same SPOC product by sector and
        # four padded TIC directory components. The URI and hash stay constant.
        tic = match[2]
        keys = "/".join([f"s{match[1]}", *(tic[i:i+4] for i in (0,4,8,12)), name])
        mirror = "https://stpubdata.s3.amazonaws.com/tess/public/tid/" + keys
        try:
            raw = download(mirror, 32_000_000)
        except (HTTPError, URLError, TimeoutError):
            raw = download(MAST + name, 32_000_000)
    sha = hashlib.sha256(raw).hexdigest()
    if product.get("sha256") and product["sha256"] != sha:
        raise ValueError(f"Checksum mismatch for {name}")
    # Verify identity before making a new cache entry available to another worker.
    from io import BytesIO
    with fits.open(BytesIO(raw), memmap=False) as hdus:
        header = hdus[0].header
        if int(header["TICID"]) != int(match[2]) or int(header["SECTOR"]) != int(match[1]):
            raise ValueError(f"FITS identity mismatch for {name}")
        if not {"TIME", "PDCSAP_FLUX", "PDCSAP_FLUX_ERR", "QUALITY"}.issubset(hdus[1].columns.names):
            raise ValueError("Invalid SPOC light-curve table")
    if not path.exists():
        atomic_bytes(path, raw)
    return path, {**product, "sha256": sha, "bytes": len(raw)}


def fetch_plan(plan, cache, workers=3):
    if not 1 <= workers <= 16:
        raise ValueError("workers must be between 1 and 16")
    pinned = json.loads(json.dumps(plan))
    products = [p for target in pinned["targets"] for p in target["products"]]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for product, (_, record) in zip(products, pool.map(lambda p: cached_product(p, cache), products)):
            product.update(record)
    return pinned
