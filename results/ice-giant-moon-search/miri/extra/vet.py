"""Compare MIRI Uranus screening tracks against native point sources and alternate motion."""
import argparse,csv,json
import numpy as np
from scipy import ndimage as ndi
from pathlib import Path
from search import SIGMA,maps
from coverage import OUT,get_states,pixels,frame
from fetch import products,read_product
from miri_audit import write_csv


def photometry(data,x,y,sigma):
    x,y=float(x),float(y);width=int(np.ceil(8*sigma))+6; ix,iy=round(x),round(y)
    if x<width or x>=data['SCI'].shape[1]-width or y<width or y>=data['SCI'].shape[0]-width:return [np.nan]*7
    gy,gx=np.mgrid[-width:width+1,-width:width+1];gx=gx+ix-x;gy=gy+iy-y;rr=np.hypot(gx,gy)
    cut=np.s_[iy-width:iy+width+1,ix-width:ix+width+1];d=data['SCI'][cut];dq=data['DQ'][cut];err=data['ERR'][cut]
    valid=((dq==0)|(dq==4))&np.isfinite(d)&np.isfinite(err)&(err>0);bg=(rr>=3.5*sigma)&(rr<6*sigma)&valid
    if bg.sum()<24:return [np.nan]*7
    # Local plane and positive response at fixed native sky position, then a diagnostic free peak.
    basis=np.stack([np.ones_like(gx),gx,gy],axis=-1);sel=bg.copy()
    for _ in range(3):
        co=np.linalg.lstsq(basis[sel],d[sel],rcond=None)[0];res=d-basis@co
        med=np.median(res[sel]);mad=max(1e-8,1.4826*np.median(abs(res[sel]-med)));sel=bg&(abs(res-med)<3*mad)
        if sel.sum()<12:return [np.nan]*7
    aperture=(rr<=2.2*sigma);ok=aperture&valid
    den=ndi.gaussian_filter(valid.astype(float),sigma)
    sm=ndi.gaussian_filter(np.where(valid,res,0),sigma)/np.maximum(.01,den)
    local=(rr<3.5*sigma)&(den>.8)
    if not local.any() or ok.sum()<.8*aperture.sum():return [np.nan]*7
    j,k=np.unravel_index(np.argmax(np.where(local,sm,-np.inf)),d.shape)
    flux=res[ok].sum()*aperture.sum()/ok.sum()
    err_flux=np.sqrt(np.sum(err[ok]**2))*aperture.sum()/ok.sum()
    return [flux,err_flux,gx[j,k],gy[j,k],sm[j,k]/mad,ok.sum()/aperture.sum(),mad]


def main(cache,head):
    states=get_states();source=list(csv.DictReader((OUT/'uranus_50000_650000_leads.csv').open()))
    origin=json.loads((OUT/'uranus_50000_650000_summary.json').read_text())['epoch_jd_tdb']
    rr=np.array([float(r['a_km']) for r in source]);ph=np.deg2rad([float(r['phase_deg']) for r in source])
    rows=[];control=[];maps_by=[];astro=[];stamps=[]
    for f,p in products():
        if f['planet']!='Uranus':continue
        g=frame(f,p,head,states[f['file']],{});idx,x,y=pixels(g,rr,ph,origin);v=f['visit'];sigma=SIGMA[f['filter']]
        data,_=read_product(p,cache,True);scores,good=maps(g,p,cache);g.update(score=scores,good=good)
        factor=1.e9*float(g['w'].proj_plane_pixel_area().value)*(np.pi/180.)**2 # MJy/sr * pixel solid angle -> mJy
        xy={int(i):np.array([xx,yy]) for i,xx,yy in zip(idx,x,y)};maps_by.append((g,xy))
        for i,pos in xy.items():
            measured=photometry(data,*pos,sigma);flux,e,dx,dy,pk,frac,mad=measured
            known=[]
            for k,state in g['states'].items():
                dra=(state[0]-g['state'][0]+180)%360-180;l,r,b,t=g['skybox']
                if k==799 or not(l-.04<dra<r+.04 and b-.04<state[1]<t+.04):continue
                kx,ky=g['w'].all_world2pix(*state[:2],0,quiet=True)
                if np.isfinite(kx+ky):known.append((np.hypot(kx-pos[0],ky-pos[1]),k))
            kdist,kcode=min(known) if known else (np.nan,0)
            rc=pos[::-1,None];in_mask=bool(ndi.map_coordinates(good.astype('f4'),rc,order=1)[0]>.995)
            z=ndi.map_coordinates(scores,rc,order=1)[0]
            rows.append(dict(template=i+1,a_km=rr[i],phase_deg=np.rad2deg(ph[i]),file=f['file'],filter=f['filter'],visit=v,mjd=float(f['exp_mid_mjd_utc']),
                x=pos[0],y=pos[1],usable=int(in_mask),normalized_score=z,aperture_mJy=flux*factor,pipeline_error_mJy=e*factor,
                peak_dx_pixels=dx,peak_dy_pixels=dy,peak_over_raw_MAD=pk,good_aperture=frac,nearest_body=kcode,nearest_body_distance_px=kdist))
            if i==0 and in_mask:
                xx,yy=np.meshgrid(np.arange(-12,13)+pos[0],np.arange(-12,13)+pos[1])
                s=ndi.map_coordinates(scores,[yy,xx],order=1);gs=ndi.map_coordinates(good.astype('f4'),[yy,xx],order=1)>.995
                stamps.append((v,np.where(gs,s,np.nan)))
        for body in (701,702,703,704,705):
            state=g['states'][body];dra=(state[0]-g['state'][0]+180)%360-180;l,r,b,t=g['skybox']
            if not l<=dra<=r or not b<=state[1]<=t:continue
            px,py=g['w'].all_world2pix(*state[:2],0,quiet=True)
            if not np.isfinite(px+py) or not (0<px<g['nx']-1 and 0<py<g['ny']-1):continue
            flux,e,dx,dy,pk,frac,mad=photometry(data,px,py,sigma)
            control.append(dict(file=f['file'],visit=v,filter=f['filter'],body=body,x=px,y=py,aperture_mJy=flux*factor,peak_dx=dx,peak_dy=dy,peak_over_raw_MAD=pk))
    for v in sorted({g['f']['visit'] for g,_ in maps_by}):
        visit=[pair for pair in maps_by if pair[0]['f']['visit']==v];anchor=visit[len(visit)//2]
        for i in range(len(rr)):
            if i not in anchor[1]:continue
            ga,pa=anchor[0],anchor[1][i]
            sky=ga['w'].all_pix2world([pa],0)[0];rel=sky-ga['state'][:2]
            for g,orbit_coords in maps_by:
                if i not in orbit_coords:continue
                band=g['f']['filter'];gv=g['f']['visit']
                models={'orbit':orbit_coords[i],'fixed_detector':pa,'fixed_planet':g['w'].all_world2pix([g['state'][:2]+rel],0,quiet=True)[0],
                    'fixed_sky':g['w'].all_world2pix([sky],0,quiet=True)[0]}
                for model,pos in models.items():
                    rc=np.asarray(pos)[::-1,None]
                    ok=ndi.map_coordinates(g['good'].astype('f4'),rc,order=1,mode='constant')[0]>.995
                    z=ndi.map_coordinates(g['score'],rc,order=1,mode='constant')[0]
                    astro.append(dict(template=i+1,anchor=ga['f']['file'],anchor_visit=v,visit=gv,file=g['f']['file'],model=model,
                                      support=int(ok),score=z,x=pos[0],y=pos[1]))
    for name,d in [('uranus_native_vet',rows),('uranus_control',control),('uranus_motion',astro)]:
        if d:write_csv(OUT/(name+'.csv'),list(d[0]),d)
    shape=[]
    for v in sorted({q[0] for q in stamps}):
        stack=np.nanmedian([s for vv,s in stamps if vv==v],axis=0)
        for j,k in np.ndindex(stack.shape):
            shape.append(dict(visit=v,dx=k-12,dy=j-12,median_normalized_score=stack[j,k]))
    write_csv(OUT/'uranus_stamps.csv',list(shape[0]),shape)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('cache',type=Path);p.add_argument('headers',type=Path);a=p.parse_args();main(a.cache,a.headers)
