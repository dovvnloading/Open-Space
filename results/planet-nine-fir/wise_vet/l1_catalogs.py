"""Look for NEOWISE single-exposure point-source detections at each predicted position."""
import csv,json
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode
import numpy as np
from catalogs import HERE,query

def main():
    rows=list(csv.DictReader((HERE/'individual.csv').open()))
    jobs=[]
    for v in range(3):
        group=[r for r in rows if int(r['visit'])==v]
        ras=[float(x['ra']) for x in group];des=[float(x['dec']) for x in group];ts=[float(x['mjd']) for x in group]
        a=float(np.median(ras));d=float(np.median(des))
        radius=3600*np.max(np.hypot((np.array(ras)-a)*np.cos(np.radians(d)),np.array(des)-d))+12
        u='https://irsa.ipac.caltech.edu/cgi-bin/Gator/nph-query?'+urlencode(dict(catalog='neowiser_p1bs_psd',
             spatial='Cone',objstr=f'{a} {d}',outfmt=1,radius=radius,radunits='arcsec',
             selcols='ra,dec,mjd,scan_id,frame_num,w2mpro,w2sigmpro,w2snr,w2rchi2,cc_flags',
             constraints=f'mjd BETWEEN {min(ts)-.001:.7f} AND {max(ts)+.001:.7f}'))
        jobs.append((f'l1visit{v}',u))
    with ThreadPoolExecutor(max_workers=3) as e:list(e.map(query,jobs))
if __name__=='__main__':main()
