"""Register the two original reference cutouts and measure both difference centroids."""
import gzip
from pathlib import Path
import sys

import numpy as np
from scipy.ndimage import shift
from scipy.optimize import least_squares
from scipy.signal import correlate2d

ROOT = Path(__file__).resolve().parent


def stamp(year, kind):
    raw = gzip.decompress((ROOT / "inputs" / f"{year}-{kind}.fits.gz").read_bytes())
    # All four supplied FITS files have one 2880-byte header and one 63x63
    # big-endian, 32-bit floating-point image, with no scaling or WCS.
    cards = [raw[i:i+80].decode() for i in range(0, 2880, 80)]
    header = {c[:8].strip(): c[10:30].strip() for c in cards if c[8] == "="}
    if [header.get(k) for k in ("BITPIX", "NAXIS", "NAXIS1", "NAXIS2")] != ["-32", "2", "63", "63"]:
        raise ValueError("Unexpected FITS representation")
    return np.frombuffer(raw, offset=2880, count=63*63, dtype=">f4").reshape(63, 63).astype(float)


def centroid(data, half):
    y, x = np.indices(data.shape)
    keep = (abs(y - 31) <= half) & (abs(x - 31) <= half)

    def residual(p):
        xc, yc, amp, sigma, base = p
        distance2 = (x[keep] - xc)**2 + (y[keep] - yc)**2
        return base + amp * np.exp(-distance2 / (2 * sigma**2)) - data[keep]

    return least_squares(residual, [31, 31, data[keep].max(), 1.5, 0]).x[:2]


def register(first, second, border):
    first = first - np.median(first)
    second = second - np.median(second)
    correlation = correlate2d(second, first, mode="same")
    start = np.array(np.unravel_index(correlation.argmax(), correlation.shape)) - 31
    inner = np.s_[border:-border, border:-border]

    def residual(p):
        dy, dx, gain, base = p
        shifted = shift(first, (dy, dx), order=1)
        return (gain * shifted[inner] + base - second[inner]).ravel()

    return least_squares(residual, [*start, 1, 0], loss="soft_l1", f_scale=3).x[:2]


def measure(figure=False):
    ref1, ref2 = (stamp(year, "template") for year in (2019, 2022))
    first, second = (stamp(year, "difference") for year in (2019, 2022))
    distances = []
    for border in (5, 10, 14):
        offset = register(ref1, ref2, border)
        for half in (3, 4, 5):
            p1, p2 = centroid(first, half), centroid(second, half)
            distance = np.linalg.norm(p2 - (p1 + offset[::-1]))
            distances.append(distance)
    dy, dx = register(ref1, ref2, 10)
    p1, p2 = centroid(first, 4) + [dx, dy], centroid(second, 4)
    print(f"Reference offset (dx, dy) = ({dx:.3f}, {dy:.3f}) pixels")
    print(f"Difference centroid separation = {np.linalg.norm(p2-p1):.3f} pixels, "
          f"window sensitivity {min(distances):.3f}–{max(distances):.3f}")
    if figure:
        import matplotlib.pyplot as plt
        from matplotlib import colors
        fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0), layout="constrained")
        frames = [shift(first, (dy, dx), order=1), second]
        for ax, im, title in zip(axes, frames, ("March 2019 · SN 2019bns", "June 2022 · SN 2022lxq")):
            center = im[22:42, 21:42]
            ax.imshow(im, origin="lower", cmap="gray_r",
                      norm=colors.AsinhNorm(linear_width=np.std(im[:12]), vmin=-np.std(im[:12]), vmax=center.max()))
            for pt, col, mark in ((p1, "#e87932", "o"), (p2, "#167bbb", "s")):
                ax.scatter(*pt, facecolors="none", edgecolors=col, s=190, linewidths=1.7, marker=mark)
            ax.set(xlim=(24, 37), ylim=(25, 39), title=title, xlabel="Registered detector x (pixels)")
        axes[0].set_ylabel("Registered detector y (pixels)")
        fig.savefig(ROOT / "difference-images.png", dpi=165)
        plt.close(fig)


if __name__ == "__main__":
    measure("--figure" in sys.argv[1:])
