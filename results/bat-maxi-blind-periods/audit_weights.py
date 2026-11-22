"""Describe the concentration of statistical weights in the original archive products."""
import csv
import io
from pathlib import Path
import sys
import tarfile

import numpy as np


def main(folder):
    with (folder/"weight-audit.csv").open("w",newline="") as handle:
        writer=csv.writer(handle)
        writer.writerow(["source_id","band","days","first_mjd","last_mjd","min_error","p10_error","median_error","max_error",
                         "highest_weight_mjd","rate_on_highest_weight_day","weight_share_top_day","weight_share_top_ten_days",
                         "original_effective_days"])
        for path in sorted(folder.glob("curves-*.tar.gz")):
            with tarfile.open(path,"r|*") as archive:
                for member in archive:
                    if not member.isfile() or not member.name.endswith(".npz"):
                        continue
                    source,band=Path(member.name).stem.rsplit(".",1)
                    with np.load(io.BytesIO(archive.extractfile(member).read()),allow_pickle=False) as curve:
                        t,y,e=curve["t"],curve["y"],curve["dy"]
                    k=e.argmin()
                    w=(e[k]/e)**2
                    w/=w.sum()
                    writer.writerow([source,band,len(t),t[0],t[-1],e[k],np.quantile(e,.1),np.median(e),np.max(e),
                                     t[k],y[k],w[k],np.sort(w)[-10:].sum(),1/(w@w)])


if __name__=="__main__":
    main(Path(sys.argv[1]))
