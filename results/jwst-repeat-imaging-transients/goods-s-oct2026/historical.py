"""Recover additional 2022–2024 JADES/GOODS-S photometry at retained positions."""

from pathlib import Path
from s3 import keys
from concurrent.futures import ThreadPoolExecutor
import get_stamps as g, pandas as pd, sys
from collections import namedtuple

output = "history_stamps"
g.root = Path(output)
g.root.mkdir(exist_ok=True)
centers = pd.read_csv("sites.csv").query("site in [93,99,185,525,829]")
centers = pd.concat(
    [centers, pd.read_csv("onsets.csv").query("site==1004")[["site", "ra", "dec"]]]
)
centers.to_csv("historical_sites.csv", index=False)
g.center_file = Path("historical_sites.csv")
# Published JADES field observation numbers from DeCoursey et al.'s observing log.
q = {"01180": [7, 10, 11, 15, 17, 18], "03215": [], "06541": []}


def list_remote(item):
    pid, obs = item
    l = list(keys(f"jwst/public/jw{pid}/L3/t/"))
    return [
        k
        for k, n in l
        if (not obs or any("/o" + f"{i:03}/" in k for i in obs))
        and any(
            k.endswith(f"f{f}w_i2d.fits") for f in ["150", "200", "277", "356", "444"]
        )
    ]


K = namedtuple("K", "key")
allkeys = []
with ThreadPoolExecutor(3) as ex:
    for v in ex.map(list_remote, q.items()):
        allkeys.extend(v)
pd.DataFrame({"key": allkeys}).to_csv("historical_index.csv", index=False)
print("keys", len(allkeys), flush=True)


def go(k):
    for n in range(3):
        try:
            return g.download(K(k))
        except Exception as e:
            if n == 2:
                print(k, repr(e), flush=True)


with ThreadPoolExecutor(10) as ex:
    list(ex.map(go, allkeys))
