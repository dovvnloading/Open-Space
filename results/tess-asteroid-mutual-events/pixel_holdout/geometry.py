"""Show changes in target-to-Earth and target-to-Sun directions across TESS visits."""
import csv
from pathlib import Path
from urllib.parse import urlencode
import numpy as np

ROOT=Path(__file__).resolve().parent

def main():
    with (ROOT/"geometry.csv").open() as f:
        rows=list(csv.DictReader(f))
    for number in ("29489","3177"):
        for body in ("399","10"):
            record=[r for r in rows if r["mpc"]==number and r["toward_naif"]==body]
            v=np.array([[float(r[k]) for k in ("x_au","y_au","z_au")] for r in record])
            v/=np.linalg.norm(v,axis=1)[:,None]
            angles=np.rad2deg(np.arccos(np.clip(v@v[0],-1,1)))
            print(number,"toward NAIF",body,[(r["sector"],round(angle,1)) for r,angle in zip(record,angles)])
            times=" ".join(str(r["jd_tdb"]) for r in record)
            params=dict(COMMAND=f"'{number};'",CENTER=f"'500@{body}'", TLIST=f"'{times}'",
                        EPHEM_TYPE="'VECTORS'",VEC_TABLE="'1'",REF_PLANE="'FRAME'",
                        VEC_CORR="'NONE'",OUT_UNITS="'AU-D'",CSV_FORMAT="'YES'",format="text")
            print("https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode(params))
if __name__=="__main__":
    main()
