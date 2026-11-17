from ztf_meta import c, HERE, CACHE
from query import acquire
import numpy as np, json, csv
from concurrent.futures import ThreadPoolExecutor


def choose(name):
    st = c.circular_family(name)
    rs = {
        r["dateobs_center"]: r
        for r in json.loads((CACHE / f"{name}.raw").read_text())[1:]
        if r["release_date"] <= "2026-11-17"
        and r["dateobs_center"] >= ("2013" if name == "Chen1" else "2023")
        and r["ifilter"].split()[0] == "r"
        and r["exposure"] >= (15 if name == "Chen1" else 50)
    }
    pool = []
    for k, r in rs.items():
        t = c.epoch(k)
        a, d = c.other_sky(st, t)
        ang = np.hypot(
            (a - r["ra_center"]) * np.cos(np.deg2rad(d)), d - r["dec_center"]
        )
        cover = ang < 0.96
        pool.append((r, cover, t))
    choices = []
    available = pool.copy()
    history = np.zeros(len(st[0]), int)
    for step in range(3):
        if not available:
            break

        # Favor overlapping tracks at independent epochs in images of comparable depth.
        def score(o):
            r, ok, t = o
            points = ok.sum() if not choices else np.sum(ok * (0.2 + history))
            return (
                points
                * np.sqrt(min(r["exposure"], 200) / 90)
                * (1 + 0.04 * (int(r["caldat"][:4]) - 2023))
            )

        current = (
            [
                v
                for v in available
                if v[0]["dateobs_center"] >= "2024" and v[1].sum() > 200
            ]
            if name == "Chen1" and step == 0
            else available
        )
        q = max(current, key=score)
        r, ok, t = q
        choices.append(r)
        history += ok
        available = [x for x in available if abs(x[2] - t) > 2.5]
        print(
            name,
            r["dateobs_center"],
            r["exposure"],
            ok.sum(),
            np.sum(history >= 2),
            r["ra_center"],
            r["dec_center"],
            flush=True,
        )
    return choices


if __name__ == "__main__":
    select = {nm: choose(nm) for nm in ["Chen1", "Chen2"]}
    (HERE / "selected.json").write_text(json.dumps(select, indent=2) + "\n")
    jobs = {}
    for name, rr in select.items():
        for k, r in enumerate(rr):
            jobs[str(r["md5sum"])] = (
                "https://astroarchive.noirlab.edu/api/retrieve/" + r["md5sum"] + "/"
            )
    with ThreadPoolExecutor(5) as e:
        list(e.map(lambda t: acquire(*t), jobs.items()))
