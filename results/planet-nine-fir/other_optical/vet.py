"""Revisit each flagged sky position in the candidate's other epochs, including its stationary light."""

from chen_pixels import *
from query import acquire, urlencode
from concurrent.futures import ThreadPoolExecutor
from matplotlib import pyplot as plt
from matplotlib.patches import Circle


def examine(name):
    flags = json.loads((HERE / (name + "_flags.json")).read_text())
    rows = json.loads((HERE / "selected.json").read_text())[name]
    catalog = np.loadtxt(CACHE / ("ls_" + name + ".raw"), delimiter=",", skiprows=1)
    out = []
    for r in flags:
        coords = np.array(r["ra_dec"])
        fid = r["track"]
        ra, dc = coords.mean(0)
        url = "https://www.legacysurvey.org/viewer/fits-cutout?" + urlencode(
            dict(
                ra=round(float(ra), 6),
                dec=round(float(dc), 6),
                layer="ls-dr10",
                bands="r",
                size=360,
                pixscale=0.262,
            )
        )
        nm = f"static_{name}_{fid}"
        acquire(nm, url)
        stamp, sh = fits.getdata(CACHE / (nm + ".raw"), header=True)
        tw = WCS(sh)
        fig, axs = plt.subplots(1, 4, figsize=(10, 2.6), constrained_layout=True)
        data = []
        for k, q in enumerate(rows):
            md = q["md5sum"]
            files = fits.open(CACHE / (md + ".raw"))
            mask = fits.open(CACHE / (md + "_mask.raw"))
            vals = []
            for j, (a, b) in enumerate(coords):
                res = None
                for ext in range(1, len(files)):
                    h = files[ext].header
                    w = WCS(h)
                    x, y = w.world_to_pixel_values(a, b)
                    if not (40 < x < h["NAXIS1"] - 40 and 40 < y < h["NAXIS2"] - 40):
                        continue
                    im = files[ext].data
                    dq = mask[ext].data
                    cal = stars_on_ccd(im, dq, w, catalog)
                    if cal is None:
                        continue
                    ps, shift, zp, *_ = cal
                    res = stamp_measure(im, dq, x + shift[0], y + shift[1], ps)
                    if res is not None:
                        res = float(res[0] * 10 ** (0.4 * (22.5 - zp)))
                    break
                vals.append(res)
            files.close()
            mask.close()
            data.append(vals)
        # Sky-registered visual checks use common, contemporaneous 2014/2016 or independent DR10 reference.
        xs, ys = tw.world_to_pixel_values(coords[:, 0], coords[:, 1])
        bg = np.median(stamp)
        sp = mad(stamp)
        axs[3].imshow(
            np.arcsinh((stamp - bg) / sp), origin="lower", cmap="gray", vmin=-2, vmax=5
        )
        for j, (xp, yp) in enumerate(zip(xs, ys)):
            axs[3].add_patch(
                Circle(
                    (xp, yp), 5, fill=False, color=["cyan", "orange", "lime"][j], lw=0.9
                )
            )
        axs[3].set_title("Legacy DR10 r")
        for k, q in enumerate(rows):
            md = q["md5sum"]
            f = fits.open(CACHE / (md + ".raw"))
            im = None
            for hdu in f[1:]:
                w = WCS(hdu.header)
                x, y = w.world_to_pixel_values(ra, dc)
                if (
                    180 < x < hdu.header["NAXIS1"] - 180
                    and 180 < y < hdu.header["NAXIS2"] - 180
                ):
                    yy0, xx0 = np.indices(stamp.shape)
                    rra, dd = tw.pixel_to_world_values(xx0, yy0)
                    xx1, yy1 = w.world_to_pixel_values(rra, dd)
                    im = nd.map_coordinates(hdu.data, [yy1, xx1], order=1)
                    break
            if im is not None:
                axs[k].imshow(
                    np.arcsinh((im - np.median(im)) / mad(im)),
                    origin="lower",
                    cmap="gray",
                    vmin=-2,
                    vmax=5,
                )
                for j, (xp, yp) in enumerate(zip(xs, ys)):
                    axs[k].add_patch(
                        Circle(
                            (xp, yp),
                            5,
                            fill=False,
                            color=["cyan", "orange", "lime"][j],
                            lw=0.9,
                        )
                    )
            axs[k].set_title(q["dateobs_center"][:10])
            f.close()
        for ax in axs:
            ax.set(xticks=[], yticks=[])
        fig.savefig(HERE / f"{name}_{fid}_vet.jpg", dpi=100)
        plt.close(fig)
        out.append(
            dict(
                track=fid,
                position_epochs=[a["dateobs_center"] for a in rows],
                flux_nmgy_images_by_positions=data,
            )
        )
    (HERE / (name + "_stationarity.json")).write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    import sys

    for nm in sys.argv[1:]:
        examine(nm)
