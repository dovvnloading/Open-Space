from urllib.request import urlopen
from urllib.parse import urlencode
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

RA, DEC = 150.122033, 2.323386208


def query(level):
    q = f"SELECT TOP 20000 obs_id,t_min,t_max,access_url FROM ivoa.ObsCore WHERE calib_level={level} AND instrument_name='NIRCAM/IMAGE' AND dataproduct_type='image' AND 1=CONTAINS(POINT('ICRS',{RA},{DEC}),s_region) ORDER BY t_max DESC"
    u = "https://jwst.esac.esa.int/server/tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=q)
    )
    b = urlopen(u, timeout=90).read()
    p = Path(f"level{level}.csv")
    p.write_bytes(b)
    t = pd.read_csv(p)
    print(level, len(t), len(t.obs_id.unique()), t.obs_id.iloc[:3].tolist(), flush=True)


if __name__ == "__main__":
    with ThreadPoolExecutor(2) as ex:
        list(ex.map(query, [2, 3]))
