from hostcat import read_range, url
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.cosmology import Planck18
import astropy.units as u
import numpy as np, pandas as pd, json, hashlib

exts = {x["name"]: x for x in json.load(open("host_exts.json"))}
memo = []


def get(name):
    x = exts[name]
    raw, _ = read_range(x["offset"], x["size"])
    memo.append(
        dict(
            url=url,
            extname=name,
            offset=x["offset"],
            bytes=len(raw),
            sha256=hashlib.sha256(raw).hexdigest(),
        )
    )
    h = fits.Header.fromstring(x["header"], sep="\n")
    raw += b"\0" * (-len(raw) % 2880)
    return fits.BinTableHDU.fromstring(h.tostring().encode() + raw).data


pos = get("SIZE")
zz = get("PHOTOZ")
assert np.array_equal(pos["ID"], zz["ID"])
site = pd.concat(
    [
        pd.read_csv("sites.csv"),
        pd.read_csv("onsets.csv")[["site", "ra", "dec"]],
        pd.read_csv("blind_pairs.csv")[["site", "ra", "dec"]],
    ]
)
q = SkyCoord(site.ra, site.dec, unit="deg")
c = SkyCoord(pos["RA"], pos["DEC"], unit="deg")
res = []
for rank in [1, 2, 3]:
    ix, sep, _ = q.match_to_catalog_sky(c, nthneighbor=rank)
    for i, j, d in zip(site.site, ix, sep.arcsec):
        rec = dict(
            site=i,
            neighbor_rank=rank,
            dr2_id=pos["ID"][j],
            host_ra=pos["RA"][j],
            host_dec=pos["DEC"][j],
            sep_arcsec=d,
        )
        rec.update(
            {
                s: zz[s][j]
                for s in [
                    "EAZY_z_a",
                    "EAZY_l68",
                    "EAZY_u68",
                    "EAZY_l95",
                    "EAZY_u95",
                    "EAZY_nfilt",
                    "EAZY_chisq_min",
                ]
            }
        )
        res.append(rec)
out = pd.DataFrame(res)
out.to_csv("host_neighbors.csv", index=False)
pd.DataFrame(memo).to_csv("host_input_ranges.csv", index=False)
print(
    out[(out.neighbor_rank == 1) & out.site.isin([93, 99, 185, 525, 829])]
    .round(3)
    .to_string(index=False)
)
