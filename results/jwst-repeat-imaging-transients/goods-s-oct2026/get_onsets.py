from pathlib import Path
import pandas as pd
import get_stamps as g
from concurrent.futures import ThreadPoolExecutor

g.root = Path("onset_stamps")
g.root.mkdir(exist_ok=True)
g.center_file = Path("onsets.csv")
rows = pd.concat([pd.read_csv(f"{i}_index.csv") for i in (8060, 12577)])
r = rows[rows.key.str.contains(r"/o00[12]/.*f(277|356)w_i2d\.fits")]


def run(x):
    for i in range(3):
        try:
            return g.download(x)
        except Exception as e:
            if i == 2:
                print(x.key, repr(e), flush=True)


with ThreadPoolExecutor(8) as ex:
    list(ex.map(run, r.itertuples()))
