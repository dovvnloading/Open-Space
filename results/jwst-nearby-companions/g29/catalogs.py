"""Catalogue matches and independent control-field counts; run from the campaign result folder."""
import json,sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'stellar'))
import mcc124_queries as q
q.OUT=HERE/'catalogue'
q.OUT.mkdir(exist_ok=True)
ir='https://irsa.ipac.caltech.edu/TAP/sync'
base=" WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',352.193797,5.246759,0.00139))"
control=" WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',53.12,-27.8,0.06))"
gaia="https://gea.esac.esa.int/tap-server/tap/sync"
visits=[
 ('seip',ir,"SELECT objid,ra,dec,i1_f_ap1,i1_df_ap1,i2_f_ap1,i2_df_ap1,i3_f_ap1,i4_f_ap1,m1_f_psf,m1_df_psf,m1_f_psf_bf,m1_df_psf_bf,m1_f_psf_3siglim,nbands FROM slphotdr4"+base),
 ('twomass',ir,"SELECT designation,ra,dec,j_m,h_m,k_m,ph_qual FROM fp_psc"+base),
 ('catwise',ir,"SELECT source_name,ra,dec,w1mpro,w2mpro FROM catwise_2020"+base),
 ('mips_control',ir,"SELECT objid,ra,dec,m1_f_psf,m1_df_psf,m1_fluxflag FROM slphotdr4"+control+" AND m1_f_psf >= 80 AND m1_fluxflag=0"),
 ('gaia_control',gaia,"SELECT COUNT(*) AS n FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',352.1956,5.2465,0.1)) AND phot_g_mean_mag IS NOT NULL"),
 ('gaia',gaia,q.GAIA.replace('153.57985,21.07495,0.00834','352.19561,5.24648,0.012')),
 ('confirmed_planets','https://exoplanetarchive.ipac.caltech.edu/TAP/sync',"SELECT DISTINCT hostname,pl_name FROM ps WHERE gaia_dr3_id='Gaia DR3 2660358032257156736' OR hostname IN ('G 29-38','G29-38','WD 2326+049','ZZ Psc')")]
if __name__=='__main__':
    with ThreadPoolExecutor(4) as ex:
        receipts=dict(ex.map(q.get_one,visits))
    (q.OUT/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
    if any('error' in r for r in receipts.values()):raise RuntimeError('Incomplete catalog response')
