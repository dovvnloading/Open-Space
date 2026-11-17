from query import *
from concurrent.futures import ThreadPoolExecutor


def cat(item):
    name, ra, dec = item
    q = f"SELECT ra,dec,flux_r FROM ls_dr10.tractor WHERE ra BETWEEN {ra-2.15} AND {ra+2.15} AND dec BETWEEN {dec-2.04} AND {dec+2.04} AND type='PSF' AND flux_r BETWEEN 22 AND 2200 AND maskbits=0"
    if name == "Chen1":
        q = "SELECT ra,dec,flux_r FROM ls_dr10.tractor_s WHERE ra BETWEEN 41.45 AND 43.90 AND dec BETWEEN -16.21 AND -13.82 AND type='PSF' AND flux_r BETWEEN 25 AND 1900 AND maskbits=0"
    url = "https://datalab.noirlab.edu/tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q, MAXREC=30000)
    )
    acquire("ls_" + name, url)


def masks(name):
    select = json.loads((HERE / "selected.json").read_text())[name]
    terms = [["instrument", "decam"], ["proc_type", "instcal"], ["prod_type", "dqmask"]]
    q = {
        "outfields": ["archive_filename", "md5sum", "dateobs_center"],
        "search": terms
        + [
            [
                "caldat",
                min(r["caldat"] for r in select),
                max(r["caldat"] for r in select),
            ],
            [
                "ra_center",
                min(r["ra_center"] for r in select) - 0.03,
                max(r["ra_center"] for r in select) + 0.03,
            ],
            [
                "dec_center",
                min(r["dec_center"] for r in select) - 0.03,
                max(r["dec_center"] for r in select) + 0.03,
            ],
        ],
    }
    acquire(
        "masks_" + name,
        "https://astroarchive.noirlab.edu/api/adv_search/find/?rectype=file&limit=20000",
        q,
    )
    rows = json.loads((CACHE / ("masks_" + name + ".raw")).read_text())[1:]
    keyed = {r["archive_filename"]: r for r in rows}
    todo = {}
    for r in select:
        m = keyed[r["archive_filename"].replace("_ooi_", "_ood_")]
        todo[r["md5sum"] + "_mask"] = (
            "https://astroarchive.noirlab.edu/api/retrieve/" + m["md5sum"] + "/"
        )
    with ThreadPoolExecutor(3) as ex:
        list(ex.map(lambda t: acquire(*t), todo.items()))


if __name__ == "__main__":
    with ThreadPoolExecutor(4) as e:
        jobs = [
            e.submit(cat, x)
            for x in [("Chen1", 42.676, -15.021), ("Chen2", 45.297, -16.711)]
        ] + [e.submit(masks, x) for x in ["Chen1", "Chen2"]]
        [j.result() for j in jobs]
