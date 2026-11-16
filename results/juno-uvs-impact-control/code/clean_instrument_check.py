import pathlib,numpy as np,json,hashlib
p=pathlib.Path(__file__).resolve().parents[1]/'clean-screen';x=np.load(p/'clean-photons.npz')['photons'];a=json.load(open(p/'clean-acquisition.json'));name='UVS_S02_805506529_2025192_P74OBS_V01.FIT';hk=np.load(p/(name+'.hk.npz'));anc=np.load(p/'s02-ancillary.npz');f=np.load(p/'s02-frames.npz');m=np.load(p/'s02-masks.npz')
q=(hk['EPHEMERIS_TIME']>=a['et_bounds'][0])&(hk['EPHEMERIS_TIME']<=a['et_bounds'][1]);fq=(f['ET_LAST_HACK']>=a['et_bounds'][0])&(f['ET_FIRST_HACK']<=a['et_bounds'][1]);tags=np.unique(f['FRAME_TAG'][fq]);assert len(tags)==1;mi=[i for i,v in enumerate(m['DATA_TAG']) if int(v.strip(),16)==tags[0]];assert len(mi)==1;i=mi[0]
masked=np.zeros(len(x),bool)
for k in range(1,9):
 a1,b1,c,d=[int(m[f'MSK{k}_{z}'][i]) for z in ['X_START','X_STOP','Y_START','Y_STOP']]
 if [a1,b1,c,d]!=[0,0,0,0]:masked|=(x['DETECTOR_X']>=a1)&(x['DETECTOR_X']<=b1)&(x['DETECTOR_Y']>=c)&(x['DETECTOR_Y']<=d)
on=(abs(x['JUPITER_LAT_400km'])<=90)&(x['JUPITER_LON_400km']>=0)&(x['JUPITER_LON_400km']<360)&(x['EMISSION_ANGLE']>=0)&(x['EMISSION_ANGLE']<90);night=on&(x['INCIDENCE_ANGLE']>90);long=(x['WAVELENGTH']>=170)&(x['WAVELENGTH']<200)
aq=(anc['EPHEMERIS_TIME_SLOW']>=a['et_bounds'][0])&(anc['EPHEMERIS_TIME_SLOW']<=a['et_bounds'][1]);centers=[]
for j in [1,2]:
 xs=np.unique(anc[f'STIM_{j}_X'][aq]);ys=np.unique(anc[f'STIM_{j}_Y'][aq]);assert len(xs)==len(ys)==1;centers.append([int(xs[0]),int(ys[0])])
stim=[]
for width in [0,4,8,16]:
 z=np.zeros(len(x),bool)
 for cx,cy in centers:z|=(abs(x['DETECTOR_X']-cx)<=width)&(abs(x['DETECTOR_Y']-cy)<=width)
 stim.append({'half_width_pixels':width,'all_photons':int(z.sum()),'nightside_photons':int((z&night).sum()),'nightside_longwave_photons':int((z&night&long).sum())})
keys=['SCI_OVERFLOW','HVPS_STATUS','HVPS_LIMITED','CHKSUM_ERROR','SAFETY_ACTIVE','DETECTOR_STIM','SCI_QUALITY'];r={'housekeeping_rows':int(q.sum()),'hk_values':{k:np.unique(hk[k][q]).tolist() for k in keys},'frame_quality':np.unique(f['FRAME_QUALITY'][fq]).tolist(),'actual_frame_tag':int(tags[0]),'actual_tag_matched_mask':m['MASK_NAME'][i].strip(),'photons_inside_mask':int(masked.sum()),'stim_centers':centers,'stim_box_sensitivity':stim,'limits':['Sampled overflow-clear state is not proof of perfectly lossless telemetry','Quality-byte bit semantics not independently decoded','Mask selected by actual acquisition tag, not planned wall-clock boundary','No survey completeness or physically calibrated false-alarm probability']};(p/'clean-instrument-summary.json').write_text(json.dumps(r,indent=2));print(json.dumps(r,indent=2))
