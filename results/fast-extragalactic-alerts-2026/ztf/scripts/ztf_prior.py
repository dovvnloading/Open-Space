"""365-day prior positive ZTF detections, including neighbouring broker object IDs."""
from concurrent.futures import ThreadPoolExecutor
import csv
from ztf_tap import get,ROOT

def main():
    with (ROOT/'derived'/'ztf_pairs.csv').open() as f:rows=list(csv.DictReader(f))
    def run(r):
        ra,dec,first=map(float,(r['meanra'],r['meandec'],r['firstmjd']))
        # Small object-mean cone allows investigation of changing object IDs.
        q=(f"SELECT o.oid,o.meanra,o.meandec,d.candid,d.ra,d.dec,d.mjd,d.magpsf,d.sigmapsf,d.isdiffpos,d.drb,d.fid "
           f"FROM ztf.object AS o INNER JOIN ztf.detection AS d ON o.oid=d.oid "
           f"WHERE 1=CONTAINS(POINT('ICRS',o.meanra,o.meandec),CIRCLE('ICRS',{ra:.9f},{dec:.9f},{3/3600:.10f})) "
           f"AND d.mjd>={first-365:.9f} AND d.mjd<{first-30:.9f} AND d.isdiffpos=1 AND d.sigmapsf>0 AND d.sigmapsf<=0.21714724095")
        result=get(r['oid']+'_prior3arcsec',q)
        print(r['oid'],len(result),flush=True)
        if float(r['mjdstarthist']) < first-365:
            previous=q.replace(f'd.mjd>={first-365:.9f} AND d.mjd<{first-30:.9f}',f'd.mjd<{first-365:.9f}')
            before=get(r['oid']+'_older3arcsec',previous)
            print(r['oid'],'older than 365d',len(before),flush=True)
    with ThreadPoolExecutor(max_workers=2) as p:list(p.map(run,rows))
if __name__=='__main__':main()
