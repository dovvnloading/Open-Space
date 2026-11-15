"""Load the derived CSV photometry and execute the frozen scoring script."""
import argparse
import gzip
from pathlib import Path
import subprocess
import sys
import tempfile

import numpy as np


HERE=Path(__file__).resolve().parent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--photometry",type=Path,default=HERE)
    p.add_argument("--out",type=Path,default=HERE)
    p.add_argument("--locked-code",type=Path,default=HERE/"locked",help="Directory holding holdout.py, score.py and phase_lock.csv")
    args=p.parse_args()
    with tempfile.TemporaryDirectory() as name:
        local=Path(name)
        for mpc,sector in [(29489,20),(3177,39)]:
            with gzip.open(args.photometry/f"photometry_{mpc}.csv.gz","rt") as f:
                rows=np.loadtxt(f,delimiter=",",skiprows=1)
            np.savez_compressed(local/f"stamps{mpc}_{sector}.npz",rows=rows)
        subprocess.run([sys.executable,str(args.locked_code.resolve()/"holdout.py"),str(local),"--out",str(args.out)],check=True)


if __name__ == "__main__":
    main()
