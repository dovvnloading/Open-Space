import urllib.parse, urllib.request, concurrent.futures, time, threading, os
from pathlib import Path

ROOT = Path(__file__).resolve().parent
CACHE = ROOT / "cache"
CACHE.mkdir(exist_ok=True)


def get(name, url):
    p = CACHE / name
    if not p.exists():
        for i in range(3):
            try:
                with urllib.request.urlopen(url, timeout=90) as r:
                    data = r.read()
                if b'"error"' in data[:350] or data.startswith(b"API ERROR"):
                    raise ValueError(data[:300])
                tmp = p.with_suffix("." + str(threading.get_ident()) + ".tmp")
                tmp.write_bytes(data)
                os.replace(tmp, p)
                break
            except Exception:
                if i == 2:
                    raise
                time.sleep(1)
    return p.read_bytes()


base = "https://ssd.jpl.nasa.gov/api/horizons.api?"


def hurl(command, a, b, step):
    params = dict(
        COMMAND=command,
        EPHEM_TYPE="VECTORS",
        CENTER="500@10",
        START_TIME=a,
        STOP_TIME=b,
        STEP_SIZE=step,
        OUT_UNITS="AU-D",
        VEC_TABLE="2",
        VEC_CORR="NONE",
        REF_PLANE="ECLIPTIC",
        CSV_FORMAT="YES",
        OBJ_DATA="NO",
    )
    return base + urllib.parse.urlencode(
        {"format": "text", **{k: "'" + v + "'" for k, v in params.items()}}
    )


def vec(name, cmd, a, b, step):
    import numpy as np

    s = get(name + ".txt", hurl(cmd, a, b, step)).decode()
    if "$$SOE" not in s:
        raise ValueError(s[:1500])
    rows = [
        l.strip().split(",")
        for l in s.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()
    ]
    return np.array([[float(l[i]) for i in [0, 2, 3, 4, 5, 6, 7]] for l in rows])


craft = {
    "wind": ("-8", "1995-01-01", "2026-11-15"),
    "ace": ("-92", "2006-07-03", "2026-11-15"),
    "stereo-a": ("-234", "2006-10-27", "2027-02-25"),
    "stereo-b": ("-235", "2006-10-27", "2014-10-01"),
    "ulysses": ("-55", "1990-10-07", "2009-06-30"),
    "psp": ("-96", "2018-08-13", "2027-01-01"),
    "solo": ("-144", "2020-02-11", "2029-01-01"),
    "bepi": ("-121", "2018-10-21", "2026-11-01"),
    "cassini": ("-82", "1997-10-17", "2004-07-01"),
    "maven": ("499", "2014-09-22", "2029-01-01"),
    "juno": ("-61", "2011-08-06", "2016-07-05"),
    "vex": ("299", "2006-04-12", "2014-11-28"),
    "earth": ("399", "1995-01-01", "2031-01-01"),
}
sburl = "https://ssd-api.jpl.nasa.gov/sbdb_query.api?" + urllib.parse.urlencode(
    {"fields": "spkid,pdes,full_name,epoch,q,e,i,om,w,tp,M1,K1,per", "sb-kind": "c"}
)
if __name__ == "__main__":

    def do(k):
        if k == "sbdb":
            return k, len(get("sbdb.json", sburl))
        cmd, a, b = craft[k]
        return k, len(vec(k, cmd, a, b, "1 d"))

    with concurrent.futures.ThreadPoolExecutor(max_workers=7) as p:
        for f in concurrent.futures.as_completed(
            [p.submit(do, k) for k in ["sbdb", *craft]]
        ):
            try:
                print(f.result(), flush=True)
            except Exception as e:
                print(str(e)[:300], flush=True)
