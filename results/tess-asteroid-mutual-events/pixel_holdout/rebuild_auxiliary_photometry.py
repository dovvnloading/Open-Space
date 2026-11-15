"""Reduce the two other visits separately on either side of the frozen date cut."""
import argparse
import csv
import gzip
from pathlib import Path
from auxiliary import VISITS
from locked.image import get,reduce

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir",type=Path)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    with gzip.open(args.out,"wt") as stream:
        writer=csv.writer(stream)
        writer.writerow(["mpc","sector","subset","jd_utc","reference_noise","prf_flux"])
        for key,(mpc,sector,split,*_) in VISITS.items():
            frames=get(args.workdir/("stamps"+key))
            for name,group in (("train",[x for x in frames if x[0]<split]),
                               ("test",[x for x in frames if x[0]>split+.5])):
                rows=reduce(group)[0]
                writer.writerows([mpc,sector,name,r[0],r[7],r[8]] for r in rows)
if __name__=="__main__":
    main()
