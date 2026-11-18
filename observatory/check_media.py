"""Check observatory plate bytes and source metadata. Network replay is opt-in."""

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import subprocess
from urllib.parse import parse_qs, urlparse
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / "observatory/assets"
PLATES = {"mars", "sun", "quasar", "dwarf"}


def checked_bytes(directory, receipt):
    path = directory / receipt["file"]
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != receipt["sha256"]:
        raise ValueError(f"{receipt['file']} does not match its recorded digest")
    return raw


def archive_bytes(url, size):
    request = Request(url, headers={"User-Agent": "Open-Space-Observatory/1.0 (science image replay)"})
    with urlopen(request, timeout=45) as response:
        return response.read(size + 1)


def jpeg_dimensions(raw):
    if raw[:2] != b"\xff\xd8":
        raise ValueError("Expected a JPEG image")
    offset = 2
    while offset < len(raw) and raw[offset] == 255:
        while raw[offset] == 255:
            offset += 1
        marker = raw[offset]
        offset += 1
        length = int.from_bytes(raw[offset:offset + 2], "big")
        if marker in (0xc0, 0xc1, 0xc2):
            height = int.from_bytes(raw[offset+3:offset+5], "big")
            width = int.from_bytes(raw[offset+5:offset+7], "big")
            return width, height
        if length < 2 or marker in (0xda, 0xd9):
            break
        offset += length
    raise ValueError("Missing JPEG frame dimensions")


def check_solar_sequence(directory, sun, remote):
    frames = sun["sequence"]
    if len(frames) != 8:
        raise ValueError("Expected eight hourly COR2 frames")
    for hour, frame in enumerate(frames, 15):
        name = f"20210713_{hour:02}0845_n4c2A.jpg"
        url = f"https://stereo-ssc.nascom.nasa.gov/browse/2021/07/13/ahead/cor2/1024/{name}"
        if frame["archiveUrl"] != url or frame["observedAt"] != f"2021-07-13T{hour:02}:08:45Z":
            raise ValueError(f"COR2 spacecraft, instrument, or sequence time differs for {name}")
        if hour == 19:
            if any(frame[k] != sun[k] for k in ("archiveUrl", "observedAt", "file", "sha256")):
                raise ValueError("The 19:08:45 UTC frame must reuse the published investigation cover")
        elif frame["file"] != f"cor2/{name}":
            raise ValueError(f"Expected the named archive file for {name}")
        raw = checked_bytes(directory, frame)
        if jpeg_dimensions(raw) != (1024, 1024):
            raise ValueError(f"COR2 frame {name} differs from its credited 1024-pixel size")
        if remote and hour != 19 and archive_bytes(url, len(raw)) != raw:
            raise ValueError(f"COR2 frame {name} differs from the archive")


def fits_header(raw):
    """Read the primary header and locate its file-creation card for network comparisons."""
    values = {}
    created = None
    for offset in range(0, len(raw), 80):
        card = raw[offset:offset + 80]
        key = card[:8].decode("ascii").strip()
        if key == "END":
            if created is None:
                raise ValueError("PS1 FITS file lacks its creation-date card")
            return values, created
        if card[8:10] == b"= ":
            values[key] = card[10:].decode("ascii").split("/", 1)[0].strip().strip("'").strip()
            if key == "DATE":
                created = offset
    raise ValueError("Missing FITS primary header END")


def check_quasar_field(directory, quasar, pair, query, remote):
    receipt = quasar["astrometry"]
    origin = urlparse(receipt["archiveUrl"])
    fits_query = parse_qs(origin.query, strict_parsing=True)
    if (origin.scheme != "https" or origin.netloc != "ps1images.stsci.edu"
            or origin.path != "/cgi-bin/fitscut.cgi" or receipt["file"] != "j0826-r-wcs.fits"
            or fits_query != {"ra": query["ra"], "dec": query["dec"], "size": query["size"],
                              "red": query["green"], "format": ["fits"]}):
        raise ValueError("PS1 astrometry must use the same center, size and r-band stack as the displayed cutout")

    raw = checked_bytes(directory, receipt)
    h, created = fits_header(raw)
    if remote:
        # fitscut reissues DATE (file creation time). Compare every other header
        # card, padding byte and science pixel against the unmodified stored FITS.
        fresh = archive_bytes(receipt["archiveUrl"], len(raw))
        if (fresh[created:created + 10] != b"DATE    = " or len(fresh) != len(raw)
                or fresh[:created] != raw[:created] or fresh[created+80:] != raw[created+80:]):
            raise ValueError("The archived r-band FITS differs outside its creation-date card")
    if (h["CTYPE1"] != "RA---TAN" or h["CTYPE2"] != "DEC--TAN"
            or [float(h[k]) for k in ("PC001001", "PC001002", "PC002001", "PC002002")] != [-1, 0, 0, 1]):
        raise ValueError("Unsupported PS1 WCS orientation or projection")
    size = int(query["size"][0])
    for axis in ("1", "2"):
        if (int(h["NAXIS" + axis]) != size or not math.isclose(quasar["fieldArcsec"],
                size * abs(float(h["CDELT" + axis])) * 3600, rel_tol=0, abs_tol=1e-4)):
            raise ValueError("PS1 field scale differs from the native FITS pixel grid")

    # Gnomonic projection onto the native tangent plane. FITS CRPIX is
    # one-based. Normalized screen fractions are measured from pixel edges.
    ra0, dec0 = (math.radians(float(h[k])) for k in ("CRVAL1", "CRVAL2"))
    markers = quasar["markers"]
    if len(markers) != 2 or {m["label"] for m in markers} != {"A", "B"}:
        raise ValueError("Expected two labeled A/B markers")
    by_image = {r["image"]: r for r in pair}
    for marker in markers:
        r = by_image[marker["label"]]
        ra, dec = math.radians(float(r["ra_deg"])), math.radians(float(r["dec_deg"]))
        delta = ra - ra0
        z = math.sin(dec)*math.sin(dec0) + math.cos(dec)*math.cos(dec0)*math.cos(delta)
        east = math.degrees(math.cos(dec)*math.sin(delta) / z)
        north = math.degrees((math.sin(dec)*math.cos(dec0) - math.cos(dec)*math.sin(dec0)*math.cos(delta)) / z)
        fits_x = float(h["CRPIX1"]) - east / float(h["CDELT1"])
        fits_y = float(h["CRPIX2"]) + north / float(h["CDELT2"])
        for axis, expected in (("x", (fits_x-.5)/size), ("y", 1-(fits_y-.5)/size)):
            if not (0 < expected < 1 and math.isclose(marker[axis], expected, rel_tol=0, abs_tol=1e-8)):
                raise ValueError(f"PS1 marker {marker['label']} {axis} differs from the FITS projection")


def validate(directory, *, revision="HEAD", remote=False):
    plates = json.loads((directory / "plates.json").read_text(encoding="utf-8"))
    if set(plates) != PLATES:
        raise ValueError("Expected one plate for each featured dossier.")
    for key, plate in plates.items():
        for field in ("name", "tag", "field", "color", "file", "kind", "credit", "source", "description", "sha256"):
            if not isinstance(plate.get(field), str) or not plate[field].strip():
                raise ValueError(f"{key}: missing {field}")
        name = plate["file"]
        if Path(name).name != name or "/" in name or "\\" in name:
            raise ValueError(f"{key}: expected a single image filename")
        if urlparse(plate["source"]).scheme != "https":
            raise ValueError(f"{key}: missing source URL")
        raw = checked_bytes(directory, plate)
        if key == "dwarf":
            if raw[:4] != b"RIFF" or raw[8:12] != b"WEBP":
                raise ValueError("Expected the credited WebP illustration.")
            if plate.get("archiveUrl") or plate["kind"] != "Artist’s impression":
                raise ValueError("The generated dwarf illustration requires its artwork credit.")
            continue
        url = plate["archiveUrl"]
        origin = urlparse(url)
        hosts = {"mars": "assets.science.nasa.gov", "sun": "stereo-ssc.nascom.nasa.gov",
                 "quasar": "ps1images.stsci.edu"}
        if origin.scheme != "https" or origin.hostname != hosts[key] or raw[:3] != b"\xff\xd8\xff":
            raise ValueError(f"{key}: expected an archival JPEG and HTTPS download URL")
        if key == "quasar" and jpeg_dimensions(raw) != (512, 512):
            raise ValueError("PS1 JPEG dimensions differ from the displayed survey cutout")
        if remote and archive_bytes(url, len(raw)) != raw:
            raise ValueError(f"{key}: {name} differs from bytes currently returned by the archive")

    target_rows = subprocess.check_output(
        ["git", "-C", str(ROOT), "show", f"{revision}:results/j0826-lens-delay/targets.csv"], encoding="utf-8"
    )
    pair = [r for r in csv.DictReader(io.StringIO(target_rows)) if r["system"] == "J0826+7002"]
    if len(pair) != 2 or {r["image"] for r in pair} != {"A", "B"}:
        raise ValueError("Expected the A and B lens positions for J0826+7002")
    query = parse_qs(urlparse(plates["quasar"]["archiveUrl"]).query, strict_parsing=True)
    for axis in ("ra", "dec"):
        midpoint = sum(float(r[axis + "_deg"]) for r in pair) / 2
        if len(query[axis]) != 1 or not math.isclose(float(query[axis][0]), midpoint, rel_tol=0, abs_tol=1e-8):
            raise ValueError(f"PS1 request {axis} differs from the J0826+7002 midpoint")
    for channel, band in (("red", "i"), ("green", "r"), ("blue", "g")):
        if not query[channel][0].endswith(f".stk.{band}.unconv.fits"):
            raise ValueError(f"Unexpected PS1 band in the {channel} channel")
    if query.get("size") != ["160"] or query.get("output_size") != ["512"]:
        raise ValueError("PS1 field dimensions differ from the image credit")
    check_quasar_field(directory, plates["quasar"], pair, query, remote)
    check_solar_sequence(directory, plates["sun"], remote)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--directory", type=Path, default=ASSETS, help="Source or built assets directory")
    parser.add_argument("--revision", default="HEAD", help="Git revision containing the research target table")
    parser.add_argument("--remote", action="store_true", help="Also compare the source archives with the bundled files")
    args = parser.parse_args()
    validate(args.directory, revision=args.revision, remote=args.remote)
    print("4 covers, 8 COR2 epochs, PS1 field scale and both FITS-projected sky positions passed.")
    if args.remote:
        print("10 archive JPEGs and the r-band FITS (apart from creation time) match fresh downloads.")
