"""Fixed DES peak diagnostics. Read only authorized, hash-verified frozen inputs."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
from scipy import ndimage as nd
from astropy.io import fits
from astropy.wcs import WCS
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'archives'))
import decam_difference as dd
LABELS=('DECam_605285','DECam_605571','mask_605285','mask_605571','template_605571')
CASES={'dec1':{'pixel':[608,1312],'offset':[-17,-11],'prior':'pr143-des-dec1-f6b97d8-296a210-20261117','score':5.692525918881558},'dec2':{'pixel':[609,1312],'offset':[-17,-11],'prior':'pr143-des-dec2-f6b97d8-b000d53-20261117','score':5.717945263780821}}
CONTROLS=[(-100,0),(100,0),(0,-100),(0,100)]

def sha(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()

def geometry(first,offset,w1,w2,shift):
 x,y=map(float,first);dx,dy=offset;sx,sy=shift
 # dd.nd.shift(second, [-dy,-dx]) means output q samples registered input q+shift.
 def to_native2(q):
  ra,de=w1.pixel_to_world_values(q[0]+sx,q[1]+sy)
  return tuple(map(float,w2.world_to_pixel_values(ra,de)))
 detector_sky=w2.pixel_to_world_values(x,y);rx,ry=w1.world_to_pixel_values(*detector_sky)
 return {'first_native':(x,y),'shifted_registered':(x+dx,y+dy),'sky_fixed_native2':to_native2((x,y)),
         'shifted_native2':to_native2((x+dx,y+dy)),'detector_fixed_native2':(x,y),
         'detector_fixed_registered':(float(rx-sx),float(ry-sy))}

def cut(arr,x,y,r=12):
 ix,iy=np.rint([x,y]).astype(int)
 if ix-r<0 or iy-r<0 or ix+r>=arr.shape[1] or iy+r>=arr.shape[0]:return None
 return arr[iy-r:iy+r+1,ix-r:ix+r+1].copy()

def summarize(stamp,mask):
 if stamp is None or mask is None:return {'status':'outside_image','expected_aperture_pixels':49,'valid_aperture_pixels':0,'expected_sky_pixels':192,'valid_sky_pixels':0,'minimum_valid_fraction':0.8,'plane_fit_rank':None}
 yy,xx=np.indices(stamp.shape);r=np.hypot(xx-12,yy-12);good=(mask==0)&np.isfinite(stamp);ap=(r<=4);ann=(r>=9)&(r<=12);use=good&ap;sky_mask=good&ann;sky=stamp[sky_mask]
 counts={'expected_aperture_pixels':int(ap.sum()),'valid_aperture_pixels':int(use.sum()),'expected_sky_pixels':int(ann.sum()),'valid_sky_pixels':int(sky.size),'minimum_valid_fraction':0.8,'plane_fit_rank':None}
 if use.sum()<.8*ap.sum() or sky.size<.8*ann.sum():return dict(counts,status='inconclusive_insufficient_coverage')
 med=float(np.median(sky));mad=float(1.4826*np.median(np.abs(sky-med)));flux=float(np.sum((stamp-med)[use]));w=np.maximum(stamp-med,0)*use
 design=np.column_stack((np.ones(sky_mask.sum()),xx[sky_mask]-12,yy[sky_mask]-12));coef,_,rank,_=np.linalg.lstsq(design,stamp[sky_mask],rcond=None)
 counts['plane_fit_rank']=int(rank)
 if rank<3:return dict(counts,status='inconclusive_background_rank')
 plane=coef[0]+coef[1]*(xx-12)+coef[2]*(yy-12);alternate={'radius3_median_flux':float(np.sum((stamp-med)[good&(r<=3)])),'radius4_plane_flux':float(np.sum((stamp-plane)[use])),'plane_coefficients':coef.tolist()}
 if not w.sum():return dict(counts,status='nonpositive_centroid',flux_aperture=flux,sky=med,noise_mad=mad)
 cx=float((xx*w).sum()/w.sum());cy=float((yy*w).sum()/w.sum());variance=float(((xx-cx)**2+(yy-cy)**2)[use].dot(w[use])/w.sum());sn=abs(flux)/(mad*np.sqrt(max(use.sum(),1))) if mad else 0
 return dict(counts,status='ok',flux_aperture=flux,sky_median=med,sky_mean=float(sky.mean()),noise_mad=mad,centroid_xy=[cx,cy],radial_second_moment=variance,centroid_uncertainty_uncalibrated_proxy_pixels=float(np.sqrt(variance)/max(sn,1e-12)),background_crop_sensitivity=alternate)

def reconstruction_gate(case,red,stats,score,expected,out):
 checks={'shift_max_abs':float(np.max(np.abs(np.array(red['shift'])-np.array(stats['empirical_shift_xy_pix'])))),'scale_abs':abs(red['scale']-stats['scale_b_to_a']),'star_count_equal':len(red['stars'])==stats['n_psf_stars'],'score_abs':abs(score-expected)}
 passed=(checks['shift_max_abs']<=1e-6 and checks['scale_abs']<=1e-6 and checks['star_count_equal'] and np.isfinite(checks['score_abs']) and checks['score_abs']<=1e-5)
 receipt={'case':case,'status':'passed' if passed else 'failed','checks':checks,'tolerances':{'shift_pixels':1e-6,'scale':1e-6,'score':1e-5,'star_count':'exact'}}
 if not passed:
  (out/'failed_reconstruction.json').write_text(json.dumps(receipt,indent=2)+'\n')
  raise ValueError('reconstruction mismatch: '+case)
 return receipt

def reduced(case,images,dq,hdr,w1,w2):
 from timing_scenarios import SCENARIOS
 inputs=json.loads((ROOT/'orbits/orbit_input.json').read_bytes());inputs['akari_midpoint_utc']=SCENARIOS[case][0];states=dd.p.family(inputs)
 coords=[]
 for h in hdr:
  t=dd.p.epoch(h['DATE-OBS'])+h['EXPTIME']/172800
  coords.append(dd.p.sky(states,t))
 xy=[np.array(w1.world_to_pixel_values(*c)) for c in coords];box=np.hstack(xy)
 xs=np.arange(int(np.floor(box[0].min()))-370,int(np.ceil(box[0].max()))+371);ys=np.arange(int(np.floor(box[1].min()))-370,int(np.ceil(box[1].max()))+371)
 xs=xs[(xs>=0)&(xs<images[0].shape[1])];ys=ys[(ys>=0)&(ys<images[0].shape[0])]
 a=images[0][np.ix_(ys,xs)].astype('f4');b=dd.register(images[1].astype('f4'),w2,w1,xs,ys,1)
 invalid=(dq[0][np.ix_(ys,xs)]!=0)|(dd.register((dq[1]!=0).astype('f4'),w2,w1,xs,ys,1)!=0)|~np.isfinite(a+b)
 finite=~np.isfinite(a+b);a=dd.sky_subtract(np.where(invalid,np.nan,a));b=dd.sky_subtract(np.where(invalid,np.nan,b));a=np.nan_to_num(a);b=np.nan_to_num(b)
 stars,psfs=dd.psf_stars(a,b,invalid);dx,dy=np.median(stars[:,4:],axis=0)
 b=nd.shift(b,[-dy,-dx],order=1,mode='nearest');invalid|=nd.shift(invalid.astype('f4'),[-dy,-dx],order=1,mode='nearest')>0
 stars,psfs=dd.psf_stars(a,b,invalid);scale=np.median(stars[:,2]/stars[:,3])
 bright=nd.maximum_filter(np.minimum(a,b*scale),size=27)>220
 def sharp(im):return im-nd.median_filter(im,3)>8*dd.mad(im[~invalid])+0.35*np.maximum(im,0)
 particle=nd.maximum_filter(sharp(a)|sharp(b),size=37);dilated=nd.maximum_filter(invalid,size=29)
 return {'a':a,'b':b,'invalid':invalid,'finite':finite,'bright':bright,'particle':particle,'dilated':dilated,'psfs':psfs,'stars':stars,'scale':float(scale),'shift':[float(dx),float(dy)],'origin':[int(xs[0]),int(ys[0])]}

def validate(cache,out):
 cache=cache.resolve();out=out.resolve()
 if cache.drive.upper()!='D:' or out.drive.upper()!='D:' or out.exists() or cache==out or cache in out.parents or out in cache.parents:raise ValueError('invalid D: input/output')
 receipts=json.loads((ROOT/'archives/decam_inputs.json').read_bytes())
 for label in LABELS:
  row=next(r for r in receipts if r['label']==label);p=cache/(label+'.response')
  if not p.is_file() or p.stat().st_size!=row['size'] or sha(p)!=row['sha256']:raise ValueError('input hash mismatch '+label)

def run(cache,out,prior_root):
 validate(cache,out)
 # Output path is exclusive and all native inputs have passed hashes before FITS decoding.
 out.mkdir(parents=True,exist_ok=False)
 images=[];dq=[];hdr=[]
 for ident in (605285,605571):
  im,h=fits.getdata(cache/f'DECam_{ident}.response',header=True);mask,mh=fits.getdata(cache/f'mask_{ident}.response',header=True)
  if h['EXPNUM']!=ident or mh['EXPNUM']!=ident or h['CCDNUM']!=mh['CCDNUM']:raise ValueError('exposure/CCD mismatch')
  images.append(im);dq.append(mask);hdr.append(h)
 w1,w2=map(WCS,hdr);result={};arrays={}
 for case,conf in CASES.items():
  prior=prior_root/conf['prior']/'run';stats=json.loads((prior/'difference_result.json').read_bytes());peaks=json.loads((prior/'dense_peaks.json').read_bytes())
  if len(peaks)!=1 or abs(peaks[0]['score']-conf['score'])>1e-10:raise ValueError('prior feature mismatch')
  red=reduced(case,images,dq,hdr,w1,w2)
  g=geometry(conf['pixel'],conf['offset'],w1,w2,stats['empirical_shift_xy_pix']);info={'geometry':g,'original_score':conf['score'],'original_threshold':stats['threshold'],'original_noise':stats['empirical_control_noise_adu'],'original_nearest_track':209,'nearest_track_status':'masked','stamps':{}}
  arrays[case+'_psf1']=red['psfs'][0];arrays[case+'_psf2']=red['psfs'][1];arrays[case+'_psf_stars']=red['stars']
  fx,fy=conf['pixel'];ox,oy=red['origin'];rx,ry=fx-ox,fy-oy
  for layer in ('a','b','invalid','finite','bright','particle','dilated'):
   stamp=cut(red[layer],rx,ry)
   if stamp is not None:arrays[case+'_registered_'+layer]=stamp
  scoremap=dd.filtered(red['a'],red['b'],red['psfs'],red['scale']);vx,vy=conf['offset']
  for layer in ('a','b','invalid','finite','bright','particle','dilated'):
   stamp=cut(red[layer],rx+vx,ry+vy)
   if stamp is not None:arrays[case+'_registered_shifted_'+layer]=stamp
  f0=float(scoremap[ry,rx]);f1=float(scoremap[ry+vy,rx+vx]);info['fixed_lobes']={'first':f0,'second':f1,'noise_adu':stats['empirical_control_noise_adu'],'center_adu':stats['control_median_adu'],'score':(f0-f1-stats['control_median_adu'])/stats['empirical_control_noise_adu']}
  info['reconstruction']=reconstruction_gate(case,red,stats,info['fixed_lobes']['score'],conf['score'],out)
  for hypothesis,(hx,hy) in {'sky_fixed':(rx,ry),'shifted':(rx+vx,ry+vy),'detector_fixed':(g['detector_fixed_registered'][0]-ox,g['detector_fixed_registered'][1]-oy)}.items():
   hx,hy=map(int,np.rint((hx,hy)))
   if 0<=hy<scoremap.shape[0] and 0<=hx<scoremap.shape[1]:
    second=float(scoremap[hy,hx]);info.setdefault('hypothesis_lobes',{})[hypothesis]={'first':f0,'second':second,'score_fixed_noise':(f0-second-stats['control_median_adu'])/stats['empirical_control_noise_adu']}
   else:info.setdefault('hypothesis_lobes',{})[hypothesis]={'status':'outside_range'}
  for i,(a,b) in enumerate([(0,0)]+CONTROLS):
   xx,yy=rx+a,ry+b;qx,qy=xx+vx,yy+vy
   if not (0<=xx<scoremap.shape[1] and 0<=qx<scoremap.shape[1] and 0<=yy<scoremap.shape[0] and 0<=qy<scoremap.shape[0]):
    info.setdefault('fixed_controls',[]).append({'offset':[a,b],'status':'outside_range'});continue
   masks={layer:[bool(red[layer][yy,xx]),bool(red[layer][qy,qx])] for layer in ('invalid','finite','bright','particle','dilated')}
   ff=float(scoremap[yy,xx]);gg=float(scoremap[qy,qx]);bad=any(any(masks[k]) for k in ('dilated','bright','particle'))
   info.setdefault('fixed_controls',[]).append({'offset':[a,b],'status':'masked' if bad else 'eligible','masks':masks,'lobes':[ff,gg],'score':(ff-gg-stats['control_median_adu'])/stats['empirical_control_noise_adu']})
  for name,pos1,pos2 in [('sky_fixed',g['first_native'],g['sky_fixed_native2']),('shifted',g['first_native'],g['shifted_native2']),('detector_fixed',g['first_native'],g['detector_fixed_native2'])]+[(f'control_{i}',(conf['pixel'][0]+a,conf['pixel'][1]+b),geometry((conf['pixel'][0]+a,conf['pixel'][1]+b),conf['offset'],w1,w2,stats['empirical_shift_xy_pix'])['shifted_native2']) for i,(a,b) in enumerate(CONTROLS)]:
   for epoch,pos in enumerate((pos1,pos2)):
    st=cut(images[epoch],*pos);m=cut(dq[epoch],*pos);key=f'{case}_{name}_epoch{epoch+1}';summary=summarize(st,m)
    if summary['status']=='ok':summary['centroid_native_xy']=[float(np.rint(pos[0])-12+summary['centroid_xy'][0]),float(np.rint(pos[1])-12+summary['centroid_xy'][1])]
    info['stamps'][key]=summary
    if st is not None:arrays[key]=st;arrays[key+'_dq']=m
  result[case]=info
 np.savez_compressed(out/'native_stamps.npz',**arrays)
 (out/'diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')

def main():
 a=argparse.ArgumentParser();a.add_argument('--cache',type=Path,required=True);a.add_argument('--out',type=Path,required=True);a.add_argument('--prior-root',type=Path,required=True);a.add_argument('--manifest',type=Path,required=True);a.add_argument('--manifest-sha256',required=True);a.add_argument('--auth',type=Path,required=True);a.add_argument('--auth-sha256',required=True);q=a.parse_args()
 from feature_packet import verify
 verify(q.manifest,q.manifest_sha256,q.auth,q.auth_sha256,q.cache,q.prior_root,q.out)
 run(q.cache,q.out,q.prior_root)
if __name__=='__main__':main()
