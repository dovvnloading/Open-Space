"""Reproducible acquisition of the two public HETDEX source spectra archives."""

from concurrent.futures import ThreadPoolExecutor
from functools import lru_cache
from hashlib import md5, sha256
import json
from pathlib import Path
import re
import time

from astropy.io import fits
import numpy as np
import requests


TACC = "https://web.corral.tacc.utexas.edu/hetdex/HETDEX"
MPE = "https://data.sciserver.mpe.mpg.de/public/hetdex"
SOURCES = {
    "sc1": (
        "v3.2", 7448504,
        ("c9c1be558d1e377334582fc38f8287b7", "4c117594446b8d9eb103a5e2c2a5e67b"),
        (TACC + "/catalogs/hetdex_source_catalog_1",),
    ),
    "sc2": (
        "v1.5", 19581262,
        ("d446718db26eedd5ac521c5f95ac9316", "11749c593829c42342e22e6e70b95deb"),
        (TACC + "/pdr/pdr1/hetdex_source_catalog_2", MPE + "/pdr/pdr1/hetdex_source_catalog_2"),
    ),
}
HEADERS = {"Accept-Encoding": "identity", "User-Agent": "Open-Space-HETDEX-research/1.0"}
TIMEOUT = (12, 50)
CHUNK = 16 * 1024 * 1024


def checksums(path):
    """Return publisher MD5 as an integrity check and archival SHA-256."""
    h_md5, h_sha = md5(), sha256()
    with open(path, "rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            h_md5.update(block)
            h_sha.update(block)
    return {"md5": h_md5.hexdigest(), "sha256": h_sha.hexdigest(), "bytes": path.stat().st_size}


def range_size(url):
    with requests.get(url, headers={**HEADERS, "Range": "bytes=0-0"},
                      stream=True, timeout=TIMEOUT) as response:
        response.raise_for_status()
        if response.status_code == 200:
            return None
        match = re.fullmatch(r"bytes 0-0/(\d+)", response.headers.get("Content-Range", ""))
        if not match or response.content != b"S":  # first SIMPLE card of a FITS file
            raise OSError("unexpected FITS byte-range probe/Content-Range")
        return int(match.group(1))


def download(url, path):
    """Bounded parallel range reads. Never accept a range returned at the wrong offset."""
    size = range_size(url)
    if size is None:
        with requests.get(url, headers=HEADERS, stream=True, timeout=TIMEOUT) as response:
            response.raise_for_status()
            with path.open("wb") as handle:
                for part in response.iter_content(2 * 1024 * 1024):
                    handle.write(part)
        return

    def fetch_interval(start):
        end = min(start + CHUNK, size) - 1
        for attempt in range(3):
            try:
                with requests.get(url, headers={**HEADERS, "Range": f"bytes={start}-{end}"},
                                  stream=True, timeout=TIMEOUT) as response:
                    response.raise_for_status()
                    expected = f"bytes {start}-{end}/{size}"
                    if response.status_code != 206 or response.headers.get("Content-Range") != expected:
                        raise OSError(f"range mismatch for {start}-{end}")
                    with path.open("r+b") as handle:
                        handle.seek(start)
                        for part in response.iter_content(1024 * 1024):
                            if handle.tell() + len(part) > end + 1:
                                raise OSError("oversized range")
                            handle.write(part)
                        if handle.tell() != end + 1:
                            raise OSError("truncated range")
                return
            except (OSError, requests.RequestException) as error:
                if attempt == 2:
                    raise OSError(f"{url}: {start}-{end}: {error}") from error
                time.sleep(2 ** attempt)

    with path.open("wb") as handle:
        handle.truncate(size)
    with ThreadPoolExecutor(max_workers=4) as pool:
        for _ in pool.map(fetch_interval, range(0, size, CHUNK)):
            pass


def prepare(release, cache, output):
    """Verify published MD5 for each file and persist endpoint and both checksums."""
    release = release.lower()
    version, record, expected, mirrors = SOURCES[release]
    cache, output = Path(cache), Path(output)
    cache.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    manifest_path = output / f"inputs-{release}.json"
    previous = {}
    if manifest_path.exists():
        previous = {r["file"]: r for r in json.loads(manifest_path.read_text())}
    manifest = []
    for kind, checksum in zip(("spec", "detinfo"), expected):
        name = f"hetdex_{release}_{kind}_{version}.fits"
        path = cache / name
        urls = [f"{base}/{name}" for base in mirrors]
        urls.append(f"https://zenodo.org/records/{record}/files/{name}?download=1")
        cached = checksums(path) if path.exists() else {}
        if cached.get("md5") == checksum:
            entry = {**previous.get(name, {}), **cached, "file": name, "published_md5": checksum,
                     "record": f"https://zenodo.org/records/{record}", "mirrors": urls}
        else:
            tmp = path.with_suffix(".part")
            failures = []
            for url in urls:
                try:
                    print(f"Retrieving {name} from {url}", flush=True)
                    download(url, tmp)
                    result = checksums(tmp)
                    if result["md5"] != checksum:
                        raise OSError(f"publisher MD5 mismatch: received {result['md5']}")
                    tmp.replace(path)
                    entry = {**result, "file": name, "url": url, "published_md5": checksum,
                             "record": f"https://zenodo.org/records/{record}", "mirrors": urls}
                    break
                except (OSError, requests.RequestException) as error:
                    tmp.unlink(missing_ok=True)
                    failures.append(str(error))
                    print(f"{name}: {error}", flush=True)
            else:
                raise RuntimeError(f"All mirrors failed for {name}: {failures}")
        manifest.append(entry)
        manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    return cache / manifest[0]["file"]


def plain(value):
    if isinstance(value, (np.generic,)):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("ascii", "replace").strip()
    if isinstance(value, str):
        return None if value.strip().lower() in ("", "n/a", "--") else value.strip()
    if isinstance(value, (int, float)) and (not np.isfinite(value) or value == -999):
        return None
    return value


class Spectra:
    """Mapped catalog. Metadata are consulted only after blind spectral measurement."""

    def __init__(self, path, release):
        self.path = Path(path)
        self.release = release.lower()

    def __enter__(self):
        self.hdus = fits.open(self.path, memmap=True, character_as_bytes=True)
        path_det = self.path.with_name(self.path.name.replace("_spec_", "_detinfo_"))
        self.det_hdus = fits.open(path_det, memmap=True, character_as_bytes=True)
        self.info = self.hdus["INFO"].data
        self.detinfo = self.det_hdus[1].data
        self.wave = np.asarray(self.hdus["WAVELENGTH"].data, dtype=float).ravel()
        self.flux = self.hdus["SPEC"].data
        self.error = self.hdus["SPEC_ERR"].data
        self.size = len(self.info)
        target_shape = (self.size, self.wave.size)
        if self.flux.shape == target_shape[::-1] and self.error.shape == target_shape[::-1]:
            # SC2 stores the image wavelength-major. These transposes are views.
            self.flux = self.flux.T
            self.error = self.error.T
        if self.flux.shape != target_shape or self.error.shape != target_shape:
            raise ValueError("Spectral rows and source observation table are not aligned")
        if len(self.wave) != 1036 or not np.allclose(self.wave, np.arange(3470, 5541, 2)):
            raise ValueError("Unexpected wavelength grid")
        self.det_ids = np.asarray(self.detinfo["detectid"], dtype="int64")
        self.det_order = np.argsort(self.det_ids)
        self.sorted_det_ids = self.det_ids[self.det_order]
        self.src_ids = np.asarray(self.info["source_id"], dtype="int64")
        self.src_order = np.argsort(self.src_ids)
        self.sorted_src_ids = self.src_ids[self.src_order]
        self.shot_ids = np.asarray(self.info["shotid"], dtype="int64")
        self.info_columns = {n.lower(): n for n in self.info.names}
        self.det_columns = {n.lower(): n for n in self.detinfo.names}
        self._info_vectors = {}
        self._det_vectors = {}
        return self

    def __exit__(self, typ, value, trace):
        self.metadata.cache_clear()
        self.hdus.close()
        self.det_hdus.close()

    def read(self, start, stop):
        flux = np.array(self.flux[start:stop], dtype=np.float32, copy=True, order="C")
        error = np.array(self.error[start:stop], dtype=np.float32, copy=True, order="C")
        invalid = ~np.isfinite(flux) | ~np.isfinite(error) | (error <= 0) | (flux == -999)
        flux[invalid] = np.nan
        error[invalid] = np.nan
        return flux, error

    def blocks(self, batch=256):
        for start in range(0, self.size, batch):
            yield start, *self.read(start, min(self.size, start + batch))

    def repeat_rows(self, index):
        """Exact shared source identifier in a different shot, within this release."""
        key = self.src_ids[index]
        lo = np.searchsorted(self.sorted_src_ids, key, side="left")
        hi = np.searchsorted(self.sorted_src_ids, key, side="right")
        rows = self.src_order[lo:hi]
        return rows[self.shot_ids[rows] != self.shot_ids[index]]

    def representative_detection(self, index):
        key = self._info_vectors["detectid"][index] if self._info_vectors else self.info["detectid"][index]
        pos = np.searchsorted(self.sorted_det_ids, key)
        if pos < self.sorted_det_ids.size and self.sorted_det_ids[pos] == key:
            return int(self.det_order[pos])
        return None

    @lru_cache(maxsize=512)
    def metadata(self, index):
        # Keep the joined detection columns coupled to the catalog primary line.
        # Do not confuse 'wave' with the wavelength found by the blind search.
        if not self._info_vectors:
            self._info_vectors = {k:self.info[n] for k,n in self.info_columns.items()}
            self._det_vectors = {k:self.detinfo[n] for k,n in self.det_columns.items()}
        result = {"release": self.release, "row": int(index)}
        det = self.representative_detection(index)
        if det is not None:
            result.update({k: plain(column[det]) for k, column in self._det_vectors.items()})
        result.update({k: plain(column[index]) for k, column in self._info_vectors.items()})
        if self.release == "sc1":
            for k in ("wave_group_a", "wave_group_b", "src_separation"):
                if result.get(k) is not None:
                    result[k] *= 3600
        if not result.get("ifuslot"):
            match = re.search(r"multi_\d+_(\d+)_", str(result.get("multiframe", "")))
            if match:
                result["ifuslot"] = match.group(1).zfill(3)
        return result
