"""Verify pixel orientation using predicted VIS light peaks and actual data/model correlation."""
import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from cutouts import CutoutArchive,image_extensions
from measure import HERE,table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("mge_csv",type=Path)
    p.add_argument("--index-cache",type=Path)
    args=p.parse_args()
    models={r["id_str"]:r for r in csv.DictReader(args.mge_csv.open())}
    archive=CutoutArchive("lens",args.index_cache)
    def check(obj):
        k1,k2="result/mge_lens_light.fits","result/mge_source_light.fits"
        d=archive.members(obj,[obj+".fits",k1,k2])
        raw=image_extensions(d[obj+".fits"])["VIS_FLUX"][1]
        l=image_extensions(d[k1])["VIS"][1];s=image_extensions(d[k2])["VIS"][1]
        ny,nx=l.shape
        cy,cx=(np.array(l.shape)-1)/2
        raw=raw[(raw.shape[0]-ny)//2:(raw.shape[0]+ny)//2,(raw.shape[1]-nx)//2:(raw.shape[1]+nx)//2]
        yy,xx=np.indices(l.shape)
        mask=((yy-cy)**2+(xx-cx)**2<40**2)&np.isfinite(raw)
        mod=l+s
        py,px=np.unravel_index(np.nanargmax(l),l.shape)
        ly=float(models[obj]["centre_0_median_pdf"]);lx=float(models[obj]["centre_1_median_pdf"])
        return dict(id_str=obj,light_y_arcsec=ly,light_x_arcsec=lx,peak_y_arcsec=(py-cy)*.1,peak_x_arcsec=(px-cx)*.1,
                    direct_peak_offset_arcsec=np.hypot((py-cy)*.1-ly,(px-cx)*.1-lx),
                    reversed_y_peak_offset_arcsec=np.hypot((py-cy)*.1+ly,(px-cx)*.1-lx),
                    data_model_corr=np.corrcoef(raw[mask],mod[mask])[0,1],
                    y_flipped_corr=np.corrcoef(raw[mask],mod[::-1][mask])[0,1])
    ids=(HERE/"priority.txt").read_text().split()
    with ThreadPoolExecutor(4) as pool:
        table(HERE/"registration.csv",list(pool.map(check,ids)))


if __name__=="__main__":main()
