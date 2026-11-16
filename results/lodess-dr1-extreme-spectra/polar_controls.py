"""Compare J100005's good-field flux with other isolated bright VLSSr emitters locally."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy import units as u

from release.photometry import RadioImage


def run(data, match_dir, parent, output):
    ref = pd.read_csv(parent).query("field == 'P139+79'").copy()
    target = SkyCoord(150.024625*u.deg,81.45075*u.deg)
    ref["distance_deg"] = target.separation(SkyCoord(ref.ra.to_numpy(),ref.dec.to_numpy(),unit="deg")).deg
    ref = ref[(ref.distance_deg < 3) & (ref.distance_deg > .02) & (ref.closest_as < 30)].copy()
    for survey in ("lotss","nvss","wenss"):
        tab=pd.read_csv(match_dir/f"match_{survey}.csv.gz",index_col="row")
        ref=ref.join(tab,on="closest_row")
    ref=ref[(ref.lotss_sep_as<15)&(ref.lotss_flux_jy/ref.lotss_sum135_jy>.9)&
            (ref.lotss_dc_maj<15)&(ref.nvss_sep_as<20)&(ref.wenss_sep_as<20)&
            (ref.wenss_source_type=="S")].copy()
    flux=[]
    with RadioImage(data/"images/P139+79-I.fits") as image:
        for row in ref.itertuples():
            m=image.measure(row.lotss_ra,row.lotss_dec)
            flux.append(np.nan if m is None else m["psf_jy"])
    ref["s23"]=flux
    ref["s74"]=ref.flux_74_mhz/1000
    ref["s23_over_144"]=ref.s23/ref.lotss_flux_jy
    ref["s74_over_144"]=ref.s74/ref.lotss_flux_jy
    ref[["name","closest_row","distance_deg","r_deg","lotss_name","s23","s74","lotss_flux_jy",
         "wenss_flux_jy","nvss_flux_jy","s23_over_144","s74_over_144"]].to_csv(output,index=False)
    print(ref[["s23_over_144","s74_over_144"]].describe())


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for arg in ["data","matches","vlssr_parent","output"]:
        p.add_argument(arg,type=Path)
    a=p.parse_args()
    run(a.data,a.matches,a.vlssr_parent,a.output)
