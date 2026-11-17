"""Phase-guarded FITS metadata validation before any EVENTS table decode."""
import numpy as np
from astropy.io import fits

def metadata(h,detector,trigger_id,trigger,start,stop):
 a=h[0].header
 if a.get('DETNAM')!=detector or not str(a.get('OBJECT','')).endswith(trigger_id):raise ValueError('detector or trigger identity')
 if a.get('TIMESYS')!='TT' or a.get('TIMEUNIT')!='s' or a.get('MJDREFI')!=51910 or not np.isclose(a.get('MJDREFF',0),.0007428703703703703,rtol=0,atol=1e-15) or not np.isclose(a.get('TRIGTIME',0),trigger,rtol=0,atol=1e-6):raise ValueError('time reference')
 ev=h['EVENTS'];eb=h['EBOUNDS'];gti=h['GTI']
 if not {'TIME','PHA'}<=set(ev.columns.names or []) or not {'CHANNEL','E_MIN','E_MAX'}<=set(eb.columns.names or []) or not {'START','STOP'}<=set(gti.columns.names or []):raise ValueError('missing columns')
 if ev.columns['TIME'].unit!='s' or gti.columns['START'].unit!='s' or gti.columns['STOP'].unit!='s':raise ValueError('time column units')
 if eb.columns['E_MIN'].unit!='keV' or eb.columns['E_MAX'].unit!='keV':raise ValueError('energy units')
 ch=np.asarray(eb.data['CHANNEL']);lo=np.asarray(eb.data['E_MIN']);hi=np.asarray(eb.data['E_MAX']);gs=np.asarray(gti.data['START']);ge=np.asarray(gti.data['STOP'])
 if len(set(ch.tolist()))!=len(ch) or not np.isfinite(lo).all() or not np.isfinite(hi).all() or np.any(lo>=hi) or not np.isfinite(gs).all() or not np.isfinite(ge).all() or np.any(gs>=ge):raise ValueError('bad calibration or GTI')
 if not any(a<=trigger+start and b>=trigger+stop for a,b in zip(gs,ge)):raise ValueError('window not covered by one GTI')
 return ch,ch[(lo>=50)&(hi<=300)]

def decode(h,channels,selected,trigger,start,stop):
 ev=h['EVENTS'];times=np.asarray(ev.data['TIME']);pha=np.asarray(ev.data['PHA'])
 if not np.isfinite(times).all() or np.any(np.diff(times)<0) or not np.isin(pha,channels).all():raise ValueError('invalid photon rows')
 keep=(times>=trigger+start)&(times<trigger+stop)&np.isin(pha,selected)
 return times[keep]-trigger,pha[keep]

def validate(path,trigger,start,stop,detector='NAI_00',trigger_id='000000000'):
 with fits.open(path,memmap=False) as h:
  ch,selected=metadata(h,detector,trigger_id,trigger,start,stop);t,p=decode(h,ch,selected,trigger,start,stop);return t,p,selected
