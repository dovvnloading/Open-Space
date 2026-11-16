"""Archive-based mask/HK checks; do not reinterpret undocumented quality bytes."""
import pathlib,json,numpy as np
from astropy.io import fits
from astropy.time import Time
p=pathlib.Path(__file__).resolve().parents[1]/'newer-pj74';h=fits.Header.fromstring((p/'header-2.txt').read_text(),sep='');types={'J':'>i4','E':'>f4','D':'>f8'};d=np.dtype([(h[f'TTYPE{i}'],types[h[f'TFORM{i}'].strip()]) for i in range(1,20)]);x=np.fromfile(p/'pilot.bin',dtype=d);et=x['EPHEMERIS_TIME'];m=np.load(p/'aux-8.npz');hk=np.load(p/'aux-6.npz');frames=np.load(p/'aux-1.npz')
# This pilot lies entirely inside planned mask row2. Assert rather than generalize.
def t(s):
 return (Time(s.strip(),format='yday',scale='utc').tdb.jd-2451545)*86400
# PDS uses yyyy-dddTHH rather than Astropy's yyyy:ddd:HH.
def t(s):
 a=s.strip().replace('-',':',1).replace('T',':');return (Time(a,format='yday',scale='utc').tdb.jd-2451545)*86400
idx=[i for i in range(len(m['MASK_NAME'])) if t(m['PLANNED_START_UTC'][i])<=et[0] and t(m['PLANNED_STOP_UTC'][i])>=et[-1]];assert len(idx)==1;i=idx[0]
masked=np.zeros(len(x),bool);rect=[]
for k in range(1,9):
 a,b,c,e=[int(m[f'MSK{k}_{z}'][i]) for z in ['X_START','X_STOP','Y_START','Y_STOP']];rect.append([a,b,c,e])
 if [a,b,c,e]==[0,0,0,0]:continue
 masked|=(x['DETECTOR_X']>=a)&(x['DETECTOR_X']<=b)&(x['DETECTOR_Y']>=c)&(x['DETECTOR_Y']<=e)
q=(hk['EPHEMERIS_TIME']>=et[0])&(hk['EPHEMERIS_TIME']<=et[-1]);f=(frames['ET_LAST_HACK']>=et[0])&(frames['ET_FIRST_HACK']<=et[-1])
fields=['HVPS_STATUS','HVPS_SAFE','HVPS_LIMITED','HOT_PIXEL_MASKING','SCI_OVERFLOW','CHKSUM_ERROR','SCI_QUALITY','AVR_RAW_RATE','AVR_EVENT_RATE']
r={'planned_mask_name':m['MASK_NAME'][i].strip(),'mask_rectangles':rect,'declared_unmasked_detector_area_fraction':float(m['FRAC_UM_AREA'][i]),'photons_inside_planned_mask_rectangles':int(masked.sum()),'housekeeping_samples_within_pilot':int(q.sum()),'housekeeping_unique_values':{k:np.unique(hk[k][q]).tolist() for k in fields},'frame_quality_unique_values':np.unique(frames['FRAME_QUALITY'][f]).tolist(),'frame_tags_unique_values':np.unique(frames['FRAME_TAG'][f]).tolist(),'limits':['Planned mask is not independently confirmed actual per-photon mask state','No preceding housekeeping sample in extracted product for initial photons','High quality-byte numeric values not assumed bad/good: SIS descriptions do not decode bits','Mask check does not remove solar-reflection/scan exposure background']}
(p/'instrument-summary.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
