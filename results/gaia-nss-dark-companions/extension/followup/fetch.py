from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from urllib.request import urlopen
from urllib.parse import urlencode
import hashlib,json
p=Path(__file__).resolve().parent
ga='https://gea.esac.esa.int/tap-server/tap/sync'
sim='https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
sid='275582835635329024';ss='Gaia DR3 '+sid
ir='https://irsa.ipac.caltech.edu/TAP/sync'
q={'2mass':(ir,"SELECT designation,j_m,j_cmsig,h_m,h_cmsig,k_m,k_cmsig,ph_qual,cc_flg FROM fp_psc WHERE designation='04164135+5351288'"),'gsp_quality':(ga,'SELECT flags_gspspec,teff_gspspec_lower,teff_gspspec_upper,logg_gspspec_lower,logg_gspspec_upper,mh_gspspec_lower,mh_gspspec_upper FROM gaiadr3.astrophysical_parameters WHERE source_id='+sid),'source':(ga,'SELECT source_id,ra,dec,parallax,parallax_error,ruwe,ipd_frac_multi_peak,duplicated_source,phot_g_mean_mag,phot_bp_mean_mag,phot_rp_mean_mag,rv_amplitude_robust,radial_velocity,rv_nb_transits,has_xp_sampled FROM gaiadr3.gaia_source WHERE source_id='+sid),'star':(ga,'SELECT teff_gspphot,logg_gspphot,mh_gspphot,ag_gspphot,ebpminrp_gspphot,logg_gspspec,teff_gspspec,mh_gspspec,radius_flame,lum_flame,mass_flame,flags_flame,evolstage_flame FROM gaiadr3.astrophysical_parameters WHERE source_id='+sid),'names':(sim,"SELECT basic.main_id,basic.otype,basic.sp_type,ids.ids FROM basic JOIN ident ON ident.oidref=basic.oid JOIN ids ON ids.oidref=basic.oid WHERE ident.id='"+ss+"'"),'refs':(sim,"SELECT ref.bibcode,ref.title FROM ref JOIN has_ref ON has_ref.oidbibref=ref.oidbib JOIN ident ON ident.oidref=has_ref.oidref WHERE ident.id='"+ss+"'")}
def get(t):
 k,(url,s)=t
 try:
  u=url+'?'+urlencode(dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',QUERY=s))
  dat=urlopen(u,timeout=60).read(); (p/(k+'.csv')).write_bytes(dat);return dict(file=k+'.csv',url=u,sha256=hashlib.sha256(dat).hexdigest())
 except Exception: raise
with ThreadPoolExecutor(4) as ex: receipts=list(ex.map(get,q.items()))
(p/'inputs.json').write_text(json.dumps(receipts,indent=2)+'\n')
