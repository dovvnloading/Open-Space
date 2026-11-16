"""Exposure-support diagnostics; occupied photon bins are NOT geometric exposure."""
import pathlib,json,re,math,io,numpy as np
from astropy.io import fits
p=pathlib.Path(__file__).resolve().parents[1];o=p/'systematic';results=[]
def load_masks(name):
 f=p/'clean-screen'/(name+'.tail');b=f.read_bytes();http=f.with_suffix('.http').read_text();start=int(re.findall(r'content-range: bytes (\d+)-',http,re.I)[-1])
 for match in re.finditer(b'XTENSION=',b):
  j=match.start()
  if (start+j)%2880:continue
  e=next(k+80 for k in range(j,len(b)-80,80) if b[k:k+8]==b'END     ');hs=math.ceil((e-j)/2880)*2880;h=fits.Header.fromstring(b[j:j+hs].decode(),sep='')
  if 'mask' not in h.get('EXTNAME','').lower():continue
  size=h['NAXIS1']*h['NAXIS2'];d=fits.open(io.BytesIO(fits.PrimaryHDU().header.tostring().encode()+b[j:j+hs+math.ceil(size/2880)*2880]))[1].data
  return {k:d[k] for k in d.names}
 raise RuntimeError('No mask table in verified source tail')
for s in [2,3]:
 frames=np.load(o/f's{s}-frames.npz');anc=np.load(o/f's{s}-ancillary.npz')
 for a in json.load(open(o/f's{s}-acquisition.json')):
  lo,hi=a['et_bounds'];ha,hb=a['hk_indices'];x=np.load(o/f's{s}-{ha}-{hb}.npz')['photons'];t=x['EPHEMERIS_TIME']-lo;m=load_masks(a['file']);fq=(frames['ET_LAST_HACK']>=lo)&(frames['ET_FIRST_HACK']<hi);tags=np.unique(frames['FRAME_TAG'][fq]);masked=np.zeros(len(x),bool);masknames=[]
  for tag in tags:
   ids=[i for i,v in enumerate(m['DATA_TAG']) if int(v.strip(),16)==tag];assert len(ids)==1;i=ids[0];masknames.append(str(m['MASK_NAME'][i]).strip())
   # Conservative union over actual active tags; intervals here are tiny.
   for k in range(1,9):
    x1,x2,y1,y2=[int(m[f'MSK{k}_{z}'][i]) for z in ['X_START','X_STOP','Y_START','Y_STOP']]
    if [x1,x2,y1,y2]!=[0,0,0,0]:masked|=(x['DETECTOR_X']>=x1)&(x['DETECTOR_X']<=x2)&(x['DETECTOR_Y']>=y1)&(x['DETECTOR_Y']<=y2)
  aq=(anc['EPHEMERIS_TIME_SLOW']>=lo)&(anc['EPHEMERIS_TIME_SLOW']<=hi);stim=np.zeros(len(x),bool);centers=[]
  assert aq.any()
  for k in [1,2]:
   pairs=np.unique(np.stack([anc[f'STIM_{k}_X'][aq],anc[f'STIM_{k}_Y'][aq]],axis=1),axis=0)
   for cx,cy in pairs:
    centers.append([int(cx),int(cy)]);stim|=(abs(x['DETECTOR_X']-cx)<=16)&(abs(x['DETECTOR_Y']-cy)<=16)
  lat,lon,inc,em,w=[x[k] for k in ['JUPITER_LAT_400km','JUPITER_LON_400km','INCIDENCE_ANGLE','EMISSION_ANGLE','WAVELENGTH']];on=(abs(lat)<=90)&(lon>=0)&(lon<360)&(em>=0)&(em<90);night=on&(inc>90)&(inc<=180);long=(w>=170)&(w<200);valid=~masked&~stim;bn=np.floor(t/.02).astype(int);cells=np.floor((np.clip(lat,-90,90)+90)/2).astype(int)*180+np.floor(np.clip(lon,0,359.999)/2).astype(int)
  def scan(q):
   c=np.bincount(bn[q],minlength=int(np.ceil((hi-lo)/.02)));tr=[]
   for j in np.flatnonzero(c>=20):
    v,n=np.unique(cells[q&(bn==j)],return_counts=True)
    if n.max()/c[j]>=.6:tr.append({'bin':int(j),'counts':int(c[j]),'cell':int(v[n.argmax()]),'fraction':float(n.max()/c[j])})
   return {'max_counts':int(c.max()),'count_cut_bins':int((c>=20).sum()),'triggers':tr,'occupied_bins':int((c>0).sum())}
  groups=np.floor(t/.1).astype(int)*10000+(x['DETECTOR_Y']//16)*10+night.astype(int);idx=[np.flatnonzero(groups==g) for g in np.unique(groups)];null=[]
  for seed in range(20):
   rng=np.random.default_rng(seed);lab=long.copy()
   for ids in idx:lab[ids]=rng.permutation(lab[ids])
   z=scan(night&valid&lab);null.append({'seed':seed,'max_counts':z['max_counts'],'triggers':len(z['triggers'])})
  nt=t[night&valid];r={**a,'mask_names':masknames,'actual_tags':tags.tolist(),'frame_quality':np.unique(frames['FRAME_QUALITY'][fq]).tolist(),'masked_photons':int(masked.sum()),'stim_centers':centers,'stim_box_half_width':16,'stim_removed_photons':int(stim.sum()),'stim_removed_night_photons':int((stim&night).sum()),'on_disc_photons':int(on.sum()),'nightside_valid_photons':int((night&valid).sum()),'nightside_longwave_valid_photons':int((night&long&valid).sum()),'night_photon_time_support':[float(nt.min()),float(nt.max())] if len(nt) else [],'max_gap_between_night_photons':float(np.diff(nt).max()) if len(nt)>1 else None,'night_any_wavelength_occupied_20ms_bins':int(len(np.unique(bn[night&valid]))),'nightside':scan(night&long&valid),'dayside':scan(on&~night&long&valid),'off_disc':scan(~on&long&valid),'conditional_diagnostics':null};results.append(r)
(o/'systematic-results.json').write_text(json.dumps({'intervals':results,'limits':['Nominal time and photon-occupied bins are not effective exposure','No injection completeness or physical false-alarm calibration','Conditional wavelength labels break exact detector-X mapping','Masks unioned conservatively across actual tags; STIM unioned across sampled centers','Sampled HK does not prove loss-free observing']},indent=2));print([(r['file'][:7],r['hk_indices'],r['nightside_valid_photons'],r['nightside_longwave_valid_photons'],r['nightside']['max_counts'],len(r['nightside']['triggers'])) for r in results])
