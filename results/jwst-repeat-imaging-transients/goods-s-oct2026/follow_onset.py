import get_stamps as g
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import pandas as pd
import sys

mode = "blind" if "--blind" in sys.argv else "follow"
centers = (
    pd.read_csv("blind_pairs.csv").query("site in [2010,2027,2028]")
    if mode == "blind"
    else pd.read_csv("onsets.csv").query("site==1004")
)
centers.to_csv(mode + "up_site.csv", index=False)
g.root = Path(mode + "_stamps")
g.root.mkdir(exist_ok=True)
g.center_file = Path(mode + "up_site.csv")
rows = pd.concat([pd.read_csv(f"{i}_index.csv") for i in [8060, 12577]])
rows = rows[rows.key.str.endswith("_i2d.fits")]


def run(x):
    for i in range(3):
        try:
            return g.download(x)
        except Exception as e:
            if i == 2:
                print(x.key, repr(e), flush=True)


with ThreadPoolExecutor(8) as ex:
    list(ex.map(run, rows.itertuples()))
