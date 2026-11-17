from get_stamps import section
from catalog_screen import readcat
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
import numpy as np, pandas as pd, warnings
from scipy.ndimage import map_coordinates
from pathlib import Path
from io import BytesIO
from concurrent.futures import ThreadPoolExecutor

warnings.simplefilter("ignore", FITSFixedWarning)
site = pd.read_csv("onsets.csv").query("site==1004").iloc[0]
dx, dy = np.meshgrid(np.arange(-31, 32) * 0.04, np.arange(-31, 32) * 0.04)
out = Path("single_stamps")
out.mkdir(exist_ok=True)


def get(job):
    url = job.url
    f = out / Path(url).name.replace(".fits", ".npz")
    if f.exists():
        return
    with section(url, 0, 57599) as src:
        b = BytesIO(src.read())
    ph = fits.Header.fromfile(b)
    h = fits.Header.fromfile(b)
    off = b.tell()
    w = WCS(h, naxis=2)
    nx = h["NAXIS1"]
    ra = site.ra
    dec = site.dec
    x, y = w.all_world2pix(ra + dx / 3600 / np.cos(np.deg2rad(dec)), dec + dy / 3600, 0)
    lo = int(y.min()) - 1
    hi = int(y.max()) + 2
    with section(url, off + lo * nx * 4, off + hi * nx * 4 - 1) as src:
        a = np.frombuffer(src.read(), ">f4").reshape(-1, nx)
        im = map_coordinates(a, [y - lo, x], order=1)
    np.savez_compressed(
        f,
        im=im,
        date=ph["DATE-BEG"],
        filter=ph["FILTER"],
        x=x[31, 31],
        y=y[31, 31],
        url=url,
    )


if __name__ == "__main__":
    d = pd.read_csv("dither_fluxes.csv")
    j = d[(d.site == 1004) & (d.pid == 12577)]
    with ThreadPoolExecutor(6) as pool:
        list(pool.map(get, j.itertuples()))
    from matplotlib import pyplot as plt

    files = sorted(out.glob("*.npz"), key=lambda f: str(np.load(f)["date"]))
    fig, axes = plt.subplots(2, 3, figsize=(7, 5))
    for ax, p in zip(axes.flat, files):
        z = np.load(p)
        im = z["im"]
        xy = np.unravel_index(
            np.nanargmax(np.where(np.hypot(dx, dy) < 0.55, im, -1e30)), im.shape
        )
        print(str(z["date"]), float(z["x"]), float(z["y"]), dx[xy], dy[xy], im[xy])
        ax.imshow(
            im,
            origin="lower",
            vmin=-0.01,
            vmax=0.12,
            cmap="gray",
            extent=[-1.26, 1.26, -1.26, 1.26],
        )
        ax.add_patch(plt.Circle((0, 0), 0.15, fill=False, color="tomato", lw=1))
        ax.set_title(str(z["date"])[11:19] + " " + str(z["filter"]), fontsize=9)
    fig.tight_layout()
    fig.savefig("debug_dithers.png", dpi=110)
