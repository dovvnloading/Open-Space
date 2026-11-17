"""Cross-check the optical residuals with pre-Euclid and Gaia catalogues."""
from concurrent.futures import ThreadPoolExecutor
from urllib.parse import urlencode
from metadata import retrieve
from common import RA0,DE0

def main():
    a,d=float(RA0),float(DE0)
    urls=[]
    for label,base,sql in [
        ('gaia','https://gea.esac.esa.int/tap-server/tap/sync',
         f"SELECT source_id,ra,dec,pmra,pmdec,phot_g_mean_mag,phot_variable_flag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{a},{d},0.036))"),
        ('ls_dr10','https://datalab.noirlab.edu/tap/sync',
         f"SELECT ls_id,ra,dec,type,flux_r,flux_i,flux_z FROM ls_dr10.tractor WHERE ra BETWEEN {a-.033} AND {a+.033} AND dec BETWEEN {d-.033} AND {d+.033}"),
        ('des_dr2','https://datalab.noirlab.edu/tap/sync',
         f"SELECT coadd_object_id,ra,dec,mag_auto_r,mag_auto_i FROM des_dr2.main WHERE ra BETWEEN {a-.033} AND {a+.033} AND dec BETWEEN {d-.033} AND {d+.033}")]:
        urls.append((label,base+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=30000,QUERY=sql))))
    urls.append(('ps1','https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv?'+urlencode(dict(ra=a,dec=d,radius=.036,columns='[objID,raMean,decMean,nDetections,rMeanKronMag,iMeanKronMag]',**{'nDetections.gte':2}))))
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(retrieve,urls))
if __name__=='__main__':main()
