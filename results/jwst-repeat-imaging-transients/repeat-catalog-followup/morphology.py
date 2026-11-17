"""Measure centroid shifts and opposite-sign residuals for catalog coincidences."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from pixel_check import HERE,XX,YY


def centroid(im):
    r=np.hypot(XX,YY)
    sky=im[(r>.38)&(r<.65)]
    if np.isfinite(sky).sum()<100 or not np.isfinite(im[r<.22]).all():
        return np.nan,np.nan
    src=np.maximum(im-np.nanmedian(sky),0)*(r<.22)
    src=np.nan_to_num(src)
    if src.sum()<=0:return np.nan,np.nan
    return np.sum(src*XX)/src.sum(),np.sum(src*YY)/src.sum()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path)
    args=p.parse_args();hits=pd.read_csv(HERE/'coincidences.csv')
    pairs=pd.read_csv(HERE/'pairs.csv').merge(pd.read_csv(HERE/'screen_counts.csv'),on=['pair','band'])
    rows=[]
    for h in hits.itertuples():
        for band in h.bands.split(','):
            r=pairs[(pairs.pair==h.pair)&(pairs.band==band)].iloc[0]
            ims=[]
            for key in (r.reference,r.search):
                with np.load(args.cache/f'{h.candidate}_{key}.npz') as f:
                    ims.append(f['pixels'].astype('float64'))
            a,b=ims;u,v=centroid(a),centroid(b)
            diff=h.sign*(b-a*r.flux_scale);rad=np.hypot(XX,YY)
            bg=diff[(rad>.38)&(rad<.70)]
            z=diff[rad<.22]-np.nanmedian(bg) if np.isfinite(bg).any() else np.array([np.nan])
            noise=np.nanmedian(abs(bg-np.nanmedian(bg)))*1.4826 if np.isfinite(bg).any() else np.nan
            opposite=np.maximum(-z-2*noise,0).sum()
            same=np.maximum(z-2*noise,0).sum()
            rows.append(dict(candidate=h.candidate,band=band,centroid_shift_east_arcsec=v[0]-u[0],
                             centroid_shift_north_arcsec=v[1]-u[1],centroid_shift_arcsec=np.hypot(v[0]-u[0],v[1]-u[1]),
                             opposite_to_same_excess=opposite/same if same>0 else np.nan))
    pd.DataFrame(rows).to_csv(HERE/'morphology.csv',index=False,float_format='%.7g')


if __name__=='__main__':main()
