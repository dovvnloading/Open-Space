"""Check whether the sparse Gaia archive independently recovers an eclipse signature."""
from pathlib import Path
import csv
import numpy as np
import pandas as pd
from astropy.time import Time
from fit_time import CLOCKS, HERE
def main():
    lines=[l for l in (HERE/'data/gaia-epochs.tsv').read_text().splitlines() if l and not l.startswith('#')]
    g=pd.DataFrame(csv.DictReader([lines[0],*lines[3:]],delimiter='\t'))
    g=g[(g.GrVFlag=='0') & g.TimeG.str.strip().ne('')].copy()
    time=Time(g.TimeG.astype(float)+2455197.5,format='jd',scale='tcb').tdb.jd-2457000
    mag=g.Gmag.astype(float)
    out=pd.DataFrame({'btjd':time,'G':mag,'delta_mag':mag-mag.median()})
    for c,(p,e,w) in CLOCKS.items():
        out[c]=((time-e)/p+.5)%1-.5
    out.to_csv(HERE/'gaia-check.csv',index=False,float_format='%.8f')
if __name__=='__main__':main()
