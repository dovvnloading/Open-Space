"""Render the two selected primary-HDU float32 cutouts without an astronomy runtime."""

from __future__ import annotations

import math
import struct
import zlib

# These original vector plots contain static drawing primitives, internal
# fragment references and embedded raster data only. Keep new/changed SVGs
# out of the same-origin asset directory until their source is reviewed.
VECTOR_PLOTS = {
    "docs/figures/calcium-k.svg": "3796676a5a604daccbe18d882377c7f11026cf7a97a990de6e23c9cfc949adac",
    "docs/figures/ffi-revisit.svg": "5838a018f495ef0348043b2ca3e60c0d2b0609f11bfd34029fac4ac2ce81c67e",
    "docs/figures/holdout.svg": "6b5666850791654817645e569692687c8bc8611e0f08f72edad74f6ed1ad7dc9",
    "results/at2023ugy-preflare-fading/fading.svg": "6219b34016225a97fbd23669b1296a2513ac2588cf23621ad44b6f137beef01a",
    "results/chime-high-dm-pair/frb-pair.svg": "c9e7540056d1dfc733a5f4e18a532c714e9061c852ebf45d0811e84dcf3ab694",
    "results/gaia-dr3-supernova-survivors/sensitivity.svg": "870414ff7171758b7414ee3558b1e0ba8d0ab86bbbea8900fcf487b90ccceb16",
    "results/gleamx-j1627-pulse-separation/spacing.svg": "be6494c542bf324eba62b286a716819efb72fc5508140c66bbed166e4b3f9940",
    "results/goes-pacific-double-bolide/light-curves.svg": "d5d80575e0cbba4753b29bc490b982f8237d3f890b53118176d043deea3c32e7",
    "results/magellanic-rr-lyrae-mode-transitions/mode-history.svg": "d85c7846c972bf29335c97b9bfb1964f5519a61aad4e65a05c2733c07a5981df",
    "results/nea-pair-search/residuals.svg": "bc87e2694b423e68744a1c9d215e776657ac82dbf0a1d98a03354adc20756bf1",
    "results/sims-corner-seismoacoustic/wavefronts.svg": "aa603a10a7edc16acd33f6bb9b4dd166c91c0f4e666a7f4213677c51cc92ba9b",
    "results/terrestrial-coorbitals/relative-longitude.svg": "23e68e17d4f12ced99c510397da77dfb60e17eb795f5ddb6c8d97a8f95c4a8a0",
    "results/tess-asymmetric-dips/sector105.svg": "1e03d44e8818bdaf4da16f458e806a1853fa425630937de4957247c4b5a35fd5",
    "results/tess-sector108/sector108.svg": "00e8345e2d309f553e3fff8575b4b18631a959cd916ef4fa3d1386f93a70591c",
    "results/xrism-ngc7172/deciles.svg": "20b47540e6593a53a353dc7d921ddd49ad188d52039eb17f693bcafe51cd4473",
}


def render_plate(raw: bytes) -> tuple[bytes, str]:
    """One display pixel per FITS sample. Axis 2 increases upward in the PNG."""
    header = {}
    for offset in range(0, min(len(raw), 128_000), 80):
        card = raw[offset:offset + 80]
        key = card[:8].decode("ascii").strip()
        if key == "END":
            start = (offset + 80 + 2879) // 2880 * 2880
            break
        if card[8:10] == b"= ":
            header[key] = card[10:].decode("ascii").split("/", 1)[0].strip().strip("'").strip()
    else:
        raise ValueError("Missing primary FITS END card")
    width, height = int(header["NAXIS1"]), int(header["NAXIS2"])
    axes = int(header["NAXIS"])
    if (header.get("SIMPLE") != "T" or header.get("BITPIX") != "-32" or not 2 <= axes <= 4
            or any(header.get(f"NAXIS{i}") != "1" for i in range(3, axes + 1))
            or float(header.get("BSCALE", 1)) != 1 or float(header.get("BZERO", 0)) != 0
            or not 1 <= width <= 1024 or not 1 <= height <= 1024):
        raise ValueError("Expected an unscaled, single-plane float32 FITS cutout")
    pixels = struct.unpack_from(f">{width * height}f", raw, start)
    valid = sorted(p for p in pixels if math.isfinite(p))
    if not valid:
        raise ValueError("Cutout has no finite pixels")
    low = valid[int((len(valid) - 1) * .01)]
    high = valid[int((len(valid) - 1) * .995)]
    if high <= low:
        raise ValueError("Cutout has no display range")
    # A linear stretch is reproducible directly against the native array.
    image = bytearray()
    for y in reversed(range(height)):
        image.append(0)  # PNG filter None
        for p in pixels[y * width:(y + 1) * width]:
            image.append(round(255 * min(1, max(0, (p - low) / (high - low)))) if math.isfinite(p) else 0)

    def chunk(tag, payload):
        data = tag + payload
        return struct.pack(">I", len(payload)) + data + struct.pack(">I", zlib.crc32(data))

    data = (b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">2I5B", width, height, 8, 0, 0, 0, 0))
            + chunk(b"IDAT", zlib.compress(image, 9)) + chunk(b"IEND", b""))
    return data, header["DATE-OBS"]
