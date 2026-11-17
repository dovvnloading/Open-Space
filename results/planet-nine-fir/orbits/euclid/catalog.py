"""External multi-year DES DR2 catalogue for static-source vetoes."""
from urllib.parse import urlencode
import hashlib, json
from datetime import datetime, timezone
from urllib.request import urlopen

def probe(item):
    label,url=item
    with urlopen(url,timeout=90) as resp:
        data=resp.read()
        if data.startswith(b'<?xml'): raise ValueError(data[:800])
        metadata=dict(url=url,utc=datetime.now(timezone.utc).isoformat(),status=resp.status,
                      bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    (HERE/(label+'.response')).write_bytes(data)
    (HERE/(label+'_receipt.json')).write_text(json.dumps(metadata,indent=2)+'\n')

from dense import HERE

def main():
    query=("SELECT coadd_object_id,ra,dec,mag_auto_r,mag_auto_i FROM des_dr2.main "
           "WHERE ra BETWEEN 34.87 AND 35.03 AND dec BETWEEN -49.77 AND -49.61 "
           "AND (mag_auto_r<25 OR mag_auto_i<25)")
    probe(('field_des','https://datalab.noirlab.edu/tap/sync?'+urlencode(
        dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=100000,QUERY=query))))
    gaia="SELECT source_id,ra,dec,pmra,pmdec,phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',34.9447,-49.6886,0.10)) AND phot_g_mean_mag BETWEEN 16.5 AND 21"
    probe(('psf_gaia','https://gea.esac.esa.int/tap-server/tap/sync?'+urlencode(
        dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=gaia))))

if __name__=='__main__':main()
