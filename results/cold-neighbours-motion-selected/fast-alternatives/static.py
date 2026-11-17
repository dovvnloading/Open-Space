"""Test freely located stationary counterparts for two dense NEP W2 tracks."""
import gzip
import hashlib
import json
import argparse
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from astropy.table import Table
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parent
NAMES=["track2709p666_2","track2709p666_6"]
def main():
    args=argparse.ArgumentParser();args.add_argument("--l1b",type=Path,required=True);args.add_argument("--names",nargs="+",default=NAMES);args.add_argument("--out",type=Path,default=ROOT);a=args.parse_args()
    a.out.mkdir(parents=True,exist_ok=True)
    results={};epochs=[]
    for name in a.names:
        receipt=a.l1b/"retained"/(name+".json")
        if receipt.exists():
            z=json.loads(receipt.read_text())
        else:
            # The compact fast manifest preserves the archive URL, but its reduced
            # records need the separately published unTimely fit coordinates.
            refs={r["source"]:r for r in json.loads((a.l1b/"manifest.json").read_text())}
            z=refs[name].copy()
            track=pd.read_csv(a.l1b.parent/"tracks.csv").set_index("source_name").loc[name]
            z.update(ra=track.ra2015_deg,dec=track.dec2015_deg,pmra=track.pmra,pmdec=track.pmdec)
        path=a.l1b/"cache"/(name+".tbl.gz")
        raw=gzip.decompress(path.read_bytes()) if path.exists() else urlopen(z["archive_url"],timeout=240).read()
        science=b"\n".join(line for line in raw.splitlines() if not line.startswith(b"\\"))
        if hashlib.sha256(science).hexdigest()!=z["raw_science_rows_sha256"]:
            raise ValueError(f"Archive science rows changed for {name}")
        df=Table.read(raw.decode().splitlines(),format="ascii.ipac").to_pandas()
        good=(df.w2snr>=5)&(df.w2rchi2<3)&(df.sso_flg==0)&(df.qual_frame>=5)&df.cc_flags.astype(str).str.zfill(4).str[1].eq("0")
        t=df[good].copy()
        ra,dec,mx,my=z["ra"],z["dec"],z["pmra"],z["pmdec"]
        t["x"]=((t.ra-ra+180)%360-180)*np.cos(np.deg2rad(dec))*3600
        t["y"]=(t.dec-dec)*3600
        dt=(t.mjd-57170.)/365.25
        # Identify the dense held-out concentration which generated the initial hit counts.
        late=t.mjd>=float(z["heldout_mjd"])
        moving=t[late&(np.hypot(t.x-mx*dt,t.y-my*dt)<=2.)]
        centre=np.median(moving[["x","y"]],axis=0)
        # Localize a fixed counterpart *from pre-2019 data* within an 8-arcsec
        # field around the held-out hits. No W2 magnitude constraint.
        training=t[(t.mjd<58484.)&(np.hypot(t.x-centre[0],t.y-centre[1])<=8.)]
        pts=training[["x","y"]].to_numpy()
        if len(pts)<5:
            results[name]=dict(early_control_rows=len(pts),status="too few training rows for fixed counterpart");continue
        tree=cKDTree(pts)
        strength=tree.query_ball_point(pts,1.0,return_length=True)
        peak=pts[int(np.argmax(strength))]
        loc=np.median(pts[np.hypot(pts[:,0]-peak[0],pts[:,1]-peak[1])<=1.5],axis=0)
        fixed_dist=np.hypot(t.x-loc[0],t.y-loc[1])
        moving_dist=np.hypot(t.x-mx*dt,t.y-my*dt)
        selected=t[fixed_dist<=2].copy()
        selected['bin']=np.floor((selected.mjd-56600.)/182.625).astype(int)
        epoch_rows=[]
        for _,g in selected.groupby('bin'):
            if len(g)<5:continue
            epoch_rows.append(dict(track=name,n=len(g),mjd=g.mjd.median(),east=g.x.median(),north=g.y.median(),w2=g.w2mpro.median(),w1snr=g.w1snr.median()))
        e=pd.DataFrame(epoch_rows);epochs.extend(epoch_rows)
        if len(e)<3:
            results[name]=dict(early_control_rows=len(pts),status="too few well-sampled epochs");continue
        coeff=np.polyfit((e.mjd-58500.)/365.25,e[["east","north"]],1)
        results[name]=dict(fixed_offset_arcsec=loc.tolist(),fixed_ra=ra+loc[0]/(np.cos(np.deg2rad(dec))*3600),fixed_dec=dec+loc[1]/3600,
                           moving_rate=[mx,my],pre2019_fixed_within2=int(np.sum((t.mjd<58484.)&(fixed_dist<=2))),
                           heldout_moving_within2=len(moving),heldout_fixed_within2=int(np.sum(late&(fixed_dist<=2))),
                           heldout_matches_shared=int(np.sum(late&(moving_dist<=2)&(fixed_dist<=2))),
                           heldout_moving_within_3p5_of_fixed=int(np.sum(late&(moving_dist<=2)&(fixed_dist<=3.5))),
                           fixed_match_epoch_bins=len(e),fixed_centroid_fit_arcsec_yr=coeff[0].tolist(),
                           epoch_span_mjd=[float(e.mjd.min()),float(e.mjd.max())],
                           fixed_w2_magnitude_median=float(e.w2.median()),source_l1b_url=z["archive_url"],
                           reference_motion_magnitude=float(np.hypot(mx,my)))
        print(name,results[name],flush=True)
    pd.DataFrame(epochs).to_csv(a.out/"fixed-centroids.csv",index=False)
    (a.out/"results.json").write_text(json.dumps(results,indent=2)+"\n")
if __name__=="__main__":main()
