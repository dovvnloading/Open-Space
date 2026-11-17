"""Count primary stellar target spectra by number of integrated exposures and nights."""
from pathlib import Path
from catalog_audit import DESI, freeze

if __name__ == "__main__":
    folder = Path("catalog")
    folder.mkdir(exist_ok=True)
    q = ("SELECT coadd_numexp,coadd_numnight,count(*) AS n_stars FROM desi_dr1.zpix "
         "WHERE spectype='STAR' AND zwarn=0 AND zcat_primary='t' "
         "GROUP BY coadd_numexp,coadd_numnight ORDER BY coadd_numexp,coadd_numnight")
    rows = freeze(folder, "release_census", DESI, q, ["coadd_numexp","coadd_numnight","n_stars"])
    print("All clean STAR primary target records", sum(int(r["n_stars"]) for r in rows))
    for col in ("coadd_numexp", "coadd_numnight"):
        print(col, ">=2", sum(int(r["n_stars"]) for r in rows if int(r[col]) >= 2))
