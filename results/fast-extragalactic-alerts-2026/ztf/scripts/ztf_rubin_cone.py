"""ZTF forced observations and alerts associated with catalogued objects near a Rubin source.

This is an archival cone on ZTF object positions, distinct from a fresh IPAC
forced-photometry request at the exact sky position.
"""
from ztf_tap import get
def main():
    geom="1=CONTAINS(POINT('ICRS',o.meanra,o.meandec),CIRCLE('ICRS',342.8550704346,-11.7294969583,0.000833333333))"
    sources=get('rubin171198361324486700_ztf_objects','SELECT o.oid,o.meanra,o.meandec,o.firstmjd,o.lastmjd,o.ndet FROM ztf.object AS o WHERE '+geom)
    print('Objects',sources,flush=True)
    for table in ['detection','forced_photometry','non_detection']:
        result=get('rubin171198361324486700_ztf_'+table,'SELECT d.* FROM ztf.object AS o INNER JOIN ztf.'+table+' AS d ON o.oid=d.oid WHERE '+geom)
        print(table,len(result),flush=True)
if __name__=='__main__': main()
