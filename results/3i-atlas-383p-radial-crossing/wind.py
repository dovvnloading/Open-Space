"""Rebuild the nearby Mars wind check, streaming the signed-off daily KP product."""
from collections import defaultdict
import hashlib
from urllib.request import urlopen
import numpy as np

from orbits import INPUTS

URL = ("https://lasp.colorado.edu/maven/sdc/public/data/sci/kp/insitu/2025/10/"
       "mvn_kp_insitu_20251004_v26_r01.tab")
EXPECTED = "45f0834f15cc092265daff3720ab4b296ca0565004186fae37d21c5f9a3f32a7"


def run():
    digest = hashlib.sha256()
    hourly = defaultdict(list)
    with urlopen(URL,timeout=120) as raw:
        for line in raw:
            digest.update(line)
            if line.startswith(b"#"):
                continue
            fields = line.split()
            vx,vy,vz,qx,qy,qz,x,y,z = [float(fields[i-1]) for i in (43,45,47,44,46,48,190,191,192)]
            x,y,z = x/3390.,y/3390.,z/3390.
            if x<=0 or qx!=1 or qy!=1 or qz!=1 or vx>=0 or not np.isfinite([vx,vy,vz,x,y,z]).all():
                continue
            focus_dist = np.linalg.norm([x-.64,y,z])
            bow_dist = 2.04/(1+1.03*(x-.64)/focus_dist)
            if focus_dist > bow_dist+.15:
                hourly[fields[0][:13].decode()].append(-vx)
    if digest.hexdigest()!=EXPECTED:
        raise ValueError("MAVEN raw product digest changed")
    rows = ["hour_utc,n,min_minus_vx,p10,median,p90,max"]
    for hour,values in sorted(hourly.items()):
        percentiles = np.percentile(values,[0,10,50,90,100])
        rows.append(f"{hour},{len(values)},"+",".join(map(str,percentiles)))
    assert "\n".join(rows)+"\n" == (INPUTS/"maven_hourly.csv").read_text()
    all_values = [x for values in hourly.values() for x in values]
    print("Sample count and [min,p10,median,p90,max] km/s:",len(all_values),np.percentile(all_values,[0,10,50,90,100]))


if __name__=="__main__":
    run()
