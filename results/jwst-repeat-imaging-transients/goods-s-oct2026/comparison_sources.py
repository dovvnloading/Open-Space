"""Acquire published comparison positions from their public primary data services."""

from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen
from io import StringIO
import hashlib, pandas as pd, numpy as np
from astropy.table import Table
from scipy.spatial import cKDTree
from lxml import html

urls = {
    "jades.html": "https://arxiv.org/html/2406.05060v2",
    "chandra.tsv": "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJS/228/2/sources&-out.all&-out.max=unlimited",
    "smiles.tsv": "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJ/966/229&-out.all&-out.max=unlimited",
    "poulain.tsv": "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/MNRAS/487/4285&-out.all&-out.max=unlimited",
    "lyu-mrt.txt": "https://content.cld.iop.org/journals/0004-637X/941/2/191/revision1/apjac9e5dt3_mrt.txt",
    "3dgoods.tsv": "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?-source=J/ApJS/214/24/3dhstall&Field=GOODS-S&-out=ID,RAJ2000,DEJ2000,zsp,zpk&-out.max=unlimited",
    "gaia.csv": "https://gea.esac.esa.int/tap-server/tap/sync?"
    + urlencode(
        {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "QUERY": "select source_id as name,ra,dec,pmra,pmdec from gaiadr3.gaia_source where 1=contains(point('ICRS',ra,dec),circle('ICRS',53.16,-27.79,0.12))",
        }
    ),
    "mirror0.raw": "https://transients.fyi/api/sky/markers?ra=53.16&dec=-27.79&r=12&window=1Y&source=tns&measure=new&mag=all",
}
if __name__ == "__main__":
    import sys

    provenance = []
    for name, url in urls.items():
        path = Path(name)
        if not path.exists():
            with urlopen(
                Request(
                    url,
                    headers={
                        "User-Agent": "Mozilla/5.0 (compatible; JWST academic archive analysis)"
                    },
                ),
                timeout=180,
            ) as r:
                path.write_bytes(r.read())
        provenance.append(
            dict(
                local_file=name,
                url=url,
                bytes=path.stat().st_size,
                sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            )
        )
    pd.DataFrame(provenance).to_csv("comparison_inputs.csv", index=False)

    def vizier(name):
        lines = [
            t
            for t in Path(name).read_text().splitlines()
            if t and not t.startswith("#")
        ]
        return pd.read_csv(StringIO("\n".join(lines)), sep="\t", skiprows=[1, 2])

    for name in ["smiles", "poulain", "chandra"]:
        d = vizier(name + ".tsv")
        d["ra"] = d.RAJ2000
        d["dec"] = d.DEJ2000
        d["name"] = (
            "Luo2017 " + d.C + " " + d.Seq.astype(str)
            if name == "chandra"
            else name + " " + d.ID.astype(str)
        )
        d.to_csv(("cdfs" if name == "chandra" else name) + ".csv", index=False)
    d = Table.read("lyu-mrt.txt", format="ascii.cds").to_pandas()
    d["name"] = "Lyu2022 " + d.iloc[:, 0].astype(str)
    d["ra"] = d.RAdeg
    d["dec"] = d.DEdeg
    d.to_csv("lyu2022.csv", index=False)
    d = vizier("3dgoods.tsv")
    s = pd.read_csv("sites.csv").query("site in [93,99,185,525,829]").copy()
    xy = lambda x, ra, dec: np.c_[x[ra] * np.cos(np.deg2rad(-27.79)), x[dec]] * 3600
    dist, j = cKDTree(xy(d, "RAJ2000", "DEJ2000")).query(xy(s, "ra", "dec"))
    s["dist"] = dist
    s["3DHST_ID"] = d.ID.values[j]
    s["redshift_spec"] = d.zsp.values[j]
    s["redshift_zpk"] = d.zpk.values[j]
    s.to_csv("3dhst_redshifts.csv", index=False)
