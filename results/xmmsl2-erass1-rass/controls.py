"""Check published slew TDE positions without using them to select thresholds."""
import argparse
import json
from pathlib import Path

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord

# Source coordinates are encoded in the published designations. NGC 3599 is
# represented by the slew designation from the original two-event discovery.
TARGETS = [
    ("NGC 3599", "XMMSL1 J111527.3+180638", "11:15:27.3", "+18:06:38",
     "Esquej et al. 2007, A&A 462, L49", "2007A&A...462L..49E"),
    ("SDSS J132341.97+482701.3", "SDSS J132341.97+482701.3", "13:23:41.97", "+48:27:01.3",
     "Esquej et al. 2007, A&A 462, L49", "2007A&A...462L..49E"),
    ("SDSS J120136.02+300305.5", "SDSS J120136.02+300305.5", "12:01:36.02", "+30:03:05.5",
     "Saxton et al. 2012, A&A 541, A106", "2012A&A...541A.106S"),
    ("2MASX J06192755-6553079", "2MASX J06192755-6553079", "06:19:27.55", "-65:53:07.9",
     "Saxton et al. 2014, A&A 572, A1", "2014A&A...572A...1S"),
    ("XMMSL1 J074008.2-853927", "XMMSL1 J074008.2-853927", "07:40:08.2", "-85:39:27",
     "Saxton et al. 2017, A&A 598, A29", "2017A&A...598A..29S"),
]


def main(out):
    cat = pd.read_csv(out/"catalog_sources.csv.gz",low_memory=False)
    pos = SkyCoord(cat.ra.to_numpy()*u.deg,cat.dec.to_numpy()*u.deg)
    pairs = pd.read_csv(out/"comparisons.csv.gz",low_memory=False)
    result=[]
    for host, designation,ra,dec,citation,bibcode in TARGETS:
        c=SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
        sep=pos.separation(c).arcsec
        near=np.flatnonzero(sep<=35)
        shared=dict(host=host,designation=designation,ra_deg=c.ra.deg,dec_deg=c.dec.deg,
                    reference=citation,bibcode=bibcode)
        if len(near)==0:
            result.append(dict(shared,match="no slew entry within 35 arcsec",nearest_arcsec=sep.min()))
        for i in near:
            row=cat.iloc[i]
            comp=pairs.loc[pairs.detection_uid.eq(row.detection_uid)]
            selected=comp.loc[comp.quality_ok].sort_values("guarded_ratio",ascending=False)
            record=dict(shared,match=row["name"],detection_uid=row.detection_uid,
                        separation_arcsec=sep[i],clean=row.clean,xmm_quality=row.xmm_quality,
                        xmm_flux=row.xmm_flux,ero_flux=row.ero_flux,rass_flux=row.rass_flux,
                        quality_comparisons=len(selected))
            if len(selected):
                b=selected.iloc[0]
                record.update(high=b.high,low=b.low,low_kind=b.low_kind,
                              ratio=b.ratio,guarded_ratio=b.guarded_ratio,
                              triggered=bool(b.photometric_lead))
            result.append(record)
    pd.DataFrame(result).to_csv(out/"positive_controls.csv",index=False)
    print(json.dumps({"controls":len(TARGETS),"rows":len(result)}))


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("out",type=Path)
    main(p.parse_args().out)
