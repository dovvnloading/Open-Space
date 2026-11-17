"""Rebuild four matched high-line-of-sight-speed rows plus the available Gaia covariances."""
from pathlib import Path
import numpy as np
import pandas as pd
root=Path(__file__).resolve().parent
c=pd.read_csv(root.parent/'control/losrv.csv',dtype={'targetid':str,'source_id':str})
l,b=np.radians(c.glon),np.radians(c.glat)
c['v_gsr_los']=c.vrad+(11.1*np.cos(l)+245.6*np.sin(l))*np.cos(b)+7.25*np.sin(b)
c['delta_rv']=abs(c.vrad-299792.458*c.rr_z)
c=c.query("rvs_warn==0 and 0<vrad_err<20 and rr_spectype=='STAR' and abs(v_gsr_los)>580 and delta_rv<60 and source_id!='999999'")
c.to_csv(root/'leads.csv',index=False)
(root/'ids.txt').write_text('\n'.join(c.targetid)+'\n')
g=pd.read_csv(root/'gaia_join_raw.csv',dtype={'targetid':str,'source_id':str})
g[g.source_id.isin(c.source_id)].to_csv(root/'gaia.csv',index=False)
