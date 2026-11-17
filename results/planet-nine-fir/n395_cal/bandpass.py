"""Extract the published dashed MAGIC system-response vector from Figure 2."""
import argparse
import csv
from pathlib import Path
import tarfile

import fitz
import numpy as np

ROOT = Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arxiv_source",type=Path,help="Source archive https://arxiv.org/src/2605.26581")
    args=parser.parse_args()
    with tarfile.open(args.arxiv_source) as tar:
        pdf=tar.extractfile("f4_filter.pdf").read()
    with fitz.open(stream=pdf,filetype="pdf") as doc:
        shapes=doc[0].get_drawings()
    path,=[s for s in shapes if s.get("dashes") not in (None,"[] 0") and len(s["items"])>50]
    # Native axes rectangle: 3850--4050 Angstrom, 0--0.5 relative transmission.
    axes=shapes[1]["rect"]
    assert np.isclose(axes.height,285.82,atol=.01)
    points=[line[1] for line in path["items"]]+[path["items"][-1][2]]
    wave=np.array([385+20*(x-axes.x0)/axes.width for x,y in points])
    throughput=np.array([max(0,.5*(axes.y1-y)/axes.height) for x,y in points])
    assert np.all(np.diff(wave)>0)
    with (ROOT/"bandpass.csv").open("w",newline="") as f:
        writer=csv.writer(f);writer.writerow(["wavelength_nm","relative_system_throughput"])
        writer.writerows(zip(wave,throughput))


if __name__=="__main__":main()
