"""Check which sky directions point towards the off-axis satellite shadow."""
import argparse
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd
import spiceypy as sp
from scipy.spatial.transform import Rotation

from analysis import HERE
from ephemeris import load_kernels


def unit(v):
    return v / np.linalg.norm(v)


def projection(event, attitude, offset=0, tilt=True):
    t=sp.str2et(event.utc)
    sun=unit(sp.spkpos("SUN",t,"IAU_MARS","LT+S","-168")[0])
    moon=sp.spkpos("PHOBOS",t,"IAU_MARS","LT+S","-168")[0]
    cross_track=unit(moon-sun*np.dot(sun,moon))
    up=unit(sp.spkpos("-168",t,"IAU_MARS","NONE","499")[0])
    east=unit(np.cross([0,0,1],up))
    north=np.cross(up,east)
    # Local north-east-down. Figure 3E of Newman et al. independently identifies the azimuth zero.
    pitch,roll=(attitude.pitch_deg,attitude.roll_deg) if tilt else (0,0)
    orient=Rotation.from_euler("ZYX",[attitude.yaw_deg+offset,pitch,roll],degrees=True)
    values=[]
    for channel in range(2,9):
        az=np.deg2rad(-45*(channel-1))
        el=np.deg2rad(35 if channel==8 else 20)
        ned=orient.apply([np.cos(el)*np.cos(az),np.cos(el)*np.sin(az),-np.sin(el)])
        look=ned[0]*north+ned[1]*east-ned[2]*up
        values.append(np.dot(look,cross_track))
    return np.array(values)


def main(cache,output):
    load_kernels(cache)
    table=pd.read_csv(HERE/"geometry.csv")
    attitudes=pd.read_csv(HERE/"pointing.csv").set_index("sol")
    results=[]
    for event in table.itertuples():
        match=attitudes[attitudes.index==event.sol]
        if match.empty or int(match.iloc[0].sclk)!=event.sclk:
            continue
        att=match.iloc[0]
        nominal=projection(event,att)
        alternatives=np.array([projection(event,att,offset,tilt) for offset in (-5,0,5) for tilt in (True,False)])
        for i in range(7):
            results.append(dict(sol=event.sol,channel=i+2,projection=nominal[i],
                                range_min=alternatives[:,i].min(),range_max=alternatives[:,i].max(),
                                projection_rank=1+int((nominal>nominal[i]).sum())))
    pd.DataFrame(results).to_csv(output,index=False,float_format="%.6f")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache",type=Path,default=Path("kernels"))
    parser.add_argument("--output",type=Path,default=Path("directions.csv"))
    args=parser.parse_args()
    main(args.cache,args.output)
