"""Check every coincident trigger against its nearest similar cataloged controls."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from campaign import HERE,compare,project


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('cache',type=Path)
    args=parser.parse_args();triggers=pd.read_csv(HERE/'coincidences.csv');rows=[]
    for pair in pd.read_csv(HERE/'pairs.csv').itertuples():
        hits=triggers[(triggers.pair==pair.pair)&triggers.bands.str.split(',').map(lambda bands:pair.band in bands)]
        if hits.empty:continue
        _,frame,_=compare(pair,args.cache)
        lf=np.log(frame.flux_search_njy/frame.flux_ref_njy)
        for hit in hits.itertuples():
            dist=np.linalg.norm(project(frame,(hit.ra,hit.dec)),axis=1)
            j=np.argmin(dist)
            geom=np.sqrt(frame.flux_ref_njy*frame.flux_search_njy)
            ok=frame.eligible & (dist>1) & (dist<60) & (abs(np.log10(geom/geom[j]))<1.0) & (frame.extended_ref==frame.extended_ref[j]) & (frame.extended_search==frame.extended_search[j])
            ix=np.flatnonzero(ok); ids=ix[np.argsort(dist[ix])[:35]]
            rec=dict(candidate=hit.candidate,band=pair.band,n=len(ids),furthest_arcsec=max(dist[ids],default=np.nan))
            if len(ids)>=10:
                bias=np.median(lf[ids]);width=1.4826*np.median(abs(lf[ids]-bias))
                formal=np.hypot(frame.err_ref_njy[j]/frame.flux_ref_njy[j],frame.err_search_njy[j]/frame.flux_search_njy[j])
                rec.update(control_flux_ratio=np.exp(bias),observed_flux_ratio=np.exp(lf[j]),scatter_log=width,
                           signed_local_score=(lf[j]-bias)/max(width,formal,.025))
            rows.append(rec)
    pd.DataFrame(rows).to_csv(HERE/'local_controls.csv',index=False,float_format='%.7g')


if __name__=='__main__':main()
