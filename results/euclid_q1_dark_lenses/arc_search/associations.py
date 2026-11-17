"""Associate image knots to MER centroid photometric redshifts without assuming identity."""
import csv
import numpy as np
from photometry import save_csv
from search import HERE

def main():
    root=HERE/'validation'
    rows=list(csv.DictReader((root/'photometry.csv').open()))
    cone=list(csv.DictReader((HERE.parent/'inference/cone_matches.csv').open()))
    out=[]
    for row in rows:
        if row['band']!='VIS' or row['method'] not in ('matched_035','native_core'):continue
        name=row['id']
        neighbors=[v for v in cone if v['id']==name and v['survey']=='euclid']
        if not neighbors:continue
        ra,dec=float(row['ra']),float(row['dec'])
        cosd=np.cos(np.deg2rad(dec))
        def separation(v):
            return 3600*np.hypot(float(v['dec'])-dec,(float(v['ra'])-ra)*cosd)
        v=min(neighbors,key=separation);dist=separation(v)
        ok=dist<=.5
        out.append(dict(id=name,component=row['point'],projected_ra=round(ra,8),projected_dec=round(dec,8),
                        nearest_separation_arcsec=round(dist,4),match_within_05=int(ok),object_id=v['object_id'] if ok else '',
                        merged_phz=v['redshift'] if ok else '',interval_95=v['PHZ_95'] if ok else '',flags=v['flags'] if ok else ''))
    save_csv(root/'component_associations.csv',out)
if __name__=='__main__':main()
