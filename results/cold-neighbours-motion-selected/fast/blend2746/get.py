from pathlib import Path
import urllib.parse,urllib.request
from concurrent.futures import ThreadPoolExecutor
from time import sleep
ROOT=Path(__file__).resolve().parent
a=274.5676647347;d=66.42021838013
q=f"SELECT ra,dec,mjd,w2mpro,w2sigmpro,w2snr,w2rchi2,qual_frame,cc_flags,scan_id,frame_num,sigra,sigdec FROM neowiser_p1bs_psd WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{a},{d},0.012)) AND w2snr >= 4 ORDER BY mjd"
jobs={"l1b":"https://irsa.ipac.caltech.edu/TAP/sync?"+urllib.parse.urlencode(dict(REQUEST="doQuery",LANG="ADQL",FORMAT="csv",MAXREC=100000,QUERY=q)),
"frames":"https://irsa.ipac.caltech.edu/ibe/search/wise/neowiser/p1bm_frm?"+urllib.parse.urlencode(dict(POS=f"{a},{d}",SIZE="0.001",COLUMNS="scan_id,frame_num,band,qual_frame,mjd_obs,crval1,crval2"))}
def get(it):
 k,u=it
 for i in range(3):
  try:
   raw=urllib.request.urlopen(u,timeout=180).read(); break
  except Exception:
   if i==2:raise
   sleep(2*i+1)
 (ROOT/(k+".url")).write_text(u);(ROOT/(k+(".csv" if k=="l1b" else ".tbl"))).write_bytes(raw);print(k,len(raw),raw[:32],flush=True)
if __name__=="__main__":
 with ThreadPoolExecutor(2) as t:list(t.map(get,jobs.items()))
