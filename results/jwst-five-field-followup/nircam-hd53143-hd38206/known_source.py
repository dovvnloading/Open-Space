"""Independently centre the known bright source on each channel, rectify the SW registration error,
and test its diffraction pattern on a held-out roll. Archive the calibration input hashes.
"""
import argparse,csv,hashlib,json,sys
from datetime import datetime,timezone
from pathlib import Path
from urllib.request import urlopen
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi
from scipy.optimize import least_squares

def main(parent,cache,out):
    sys.path.insert(0,str(parent.resolve()/'nircam'))
    from nircam_reduce import median,reduce_cube,write_table
    from standard_photometry import band_jy
    from nircam_check_sources import geometry
    entries=[r for r in csv.DictReader((parent/'nircam/go8826_inputs.csv').open()) if r['filename'].startswith('jw01537026')]
    entries.append(dict(filename='F200W.vot',url='https://svo2.cab.inta-csic.es/theory/fps/fps.php?ID=JWST/NIRCam.F200W',sha256=''))
    def get(r):
        path=cache/r['filename']
        if not path.exists():
            with urlopen(r['url'],timeout=90) as h:path.write_bytes(h.read())
        with path.open('rb') as f:sha=hashlib.file_digest(f,'sha256').hexdigest()
        if r['sha256']:assert sha==r['sha256']
        return dict(filename=r['filename'],url=r['url'],sha256=sha,bytes=path.stat().st_size)
    with ThreadPoolExecutor(4) as pool:inputs=list(pool.map(get,entries))
    write_table(out/'short_calibration_inputs.csv',list(inputs[0]),inputs)
    flux_jy=band_jy(cache/'g191b2b_mod_012.fits',cache/'F200W.vot')
    def read(f):
        with fits.open(cache/f) as h:
            im=h['SCI'].data.astype(float);im[(h['DQ'].data&3)!=0]=np.nan
            return im,h['SCI'].header.copy()
    stamps=[]
    y,x=np.mgrid[-18:19,-18:19]
    for r in entries[:3]:
        data,h=read(r['filename']);im=median(data,axis=0)
        cy,cx=np.unravel_index(np.nanargmax(im),im.shape)
        core=im[cy-1:cy+2,cx-1:cx+2];dy,dx=ndi.center_of_mass(core-np.nanmin(core));cx+=dx-1;cy+=dy-1
        yy,xx=np.indices(im.shape);rad=np.hypot(xx-cx,yy-cy);im-=median(im[(rad>=35)&(rad<43)])
        val=flux_jy/(h['PIXAR_SR']*1e6)
        from nircam_reduce import fill
        stamps.append(ndi.map_coordinates(fill(im),[cy+y,cx+x],order=3)/val)
    psf=median(stamps,axis=0)
    rows=[];side=[]
    bright={r['visit']:r for r in csv.DictReader((out/'brightest_psf_fits.csv').open())}
    for name in ['hd53143_o001','hd53143_o008']:
        meta=json.loads((out/(name+'_measurement.json')).read_text());lw=bright[name]
        dat,h=read(meta['science_filename'].replace('nrcalong','nrca2'));mat,cen,scale,rad=geometry(h,dat.shape[1:])
        refs=[read(f.replace('nrcalong','nrca2'))[0] for f in meta['reference_filenames']]
        sm=median(reduce_cube(dat,refs,rad,sky_bounds=(3.9,4.5)),axis=0)
        pred=cen+np.linalg.solve(mat,[float(lw['east']),float(lw['north'])]);px,py=np.rint(pred).astype(int)
        iy,ix=np.mgrid[py-8:py+9,px-8:px+9];ok=np.isfinite(sm[iy,ix])&(np.hypot(ix-pred[0],iy-pred[1])<8)
        sel=sm[iy[ok],ix[ok]];i=np.argmax(sel);xx,yy=ix[ok][i],iy[ok][i]
        # Separate astrometry from flux fitting. The fitted constant removes residual halo.
        def residual(v):
            tx,ty,amp,base=v
            return base+amp*ndi.map_coordinates(psf,[18+iy[ok]-ty,18+ix[ok]-tx],order=3,mode='constant')-sel
        par=least_squares(residual,[xx,yy,np.max(sel)/psf.max(),np.nanmedian(sel)],bounds=([xx-1.5,yy-1.5,0,-np.inf],[xx+1.5,yy+1.5,np.inf,np.inf])).x
        short=par[2]*h['PIXAR_SR']*1e12;long=float(lw['fit_flux_uJy']);spos=mat@(par[:2]-cen);corr=np.corrcoef(sel,sel+residual(par))[0,1]
        rows.append(dict(visit=name,lw_x=lw['fit_x'],lw_y=lw['fit_y'],f444w_uJy=long,pred_sw_x=pred[0],pred_sw_y=pred[1],peak_sw_x=par[0],peak_sw_y=par[1],sw_east=spos[0],sw_north=spos[1],registration_offset_arcsec=np.linalg.norm(par[:2]-pred)*scale,
                    f200w_uJy=short,SW_fit_r=corr,model_template_jy=flux_jy,color_AB=-2.5*np.log10(short/long)))
        oth='hd53143_o008' if name.endswith('001') else 'hd53143_o001'
        bm=bright[oth];plist=list(csv.DictReader((out/(oth+'_peaks.csv')).open()));other=np.array([[int(r['x0']),int(r['y0'])] for r in plist])
        for r in csv.DictReader((out/(name+'_peaks.csv')).open()):
            pos=np.array([int(r['x0']),int(r['y0'])])
            pred=pos+np.array([float(bm['fit_x'])-float(lw['fit_x']),float(bm['fit_y'])-float(lw['fit_y'])])
            d=np.linalg.norm(other-pred,axis=1)
            side.append(dict(visit=name,x=pos[0],y=pos[1],other=oth,other_expected_x=pred[0],other_expected_y=pred[1],
                             nearest_other_crossing_pix=d.min(),x_other=other[np.argmin(d),0],y_other=other[np.argmin(d),1]))
    for name,recs in [('known_source_color.csv',rows),('sidelobe_roll_test.csv',side)]:write_table(out/name,list(recs[0]),recs)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('cache',type=Path);p.add_argument('out',type=Path);a=p.parse_args();main(a.parent,a.cache,a.out)
