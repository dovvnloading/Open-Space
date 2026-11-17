"""Cross-check every F444W annular threshold crossing using independent sky rolls/epochs and the simultaneous SW camera."""
import argparse,csv,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from scipy import ndimage as ndi
from reduce_fields import VISITS

def run(parent,out,cache):
    sys.path.insert(0,str(parent.resolve()/'nircam'))
    from nircam_check_sources import geometry,aperture,fit_primary
    from nircam_reduce import reduce_cube,offaxis_psf,median,matched,thresholds,write_table
    psfs=sorted(cache.glob('jw01537028*long_calints.fits'))
    psf,_=offaxis_psf(psfs[2:]);det,_=offaxis_psf(psfs[:2]);kernel=det[14:23,14:23]
    frames={};phot=[];fits_out=[];allpeaks={};arch=[]
    for visit,(_,__,___,____,host) in VISITS.items():
        info=json.loads((out/(visit+'_measurement.json')).read_text())
        with fits.open(cache/info['science_filename']) as hd:
            h=hd['SCI'].header.copy();prim=hd[0].header.copy()
        with np.load(out/(visit+'_rdi_maps.npz')) as f:
            med=f['residual_MJy_sr'].astype(float);amp=f['amplitude'].astype(float);cut=f['threshold'].astype(float);snr=f['snr'].astype(float)
        mat,cen,scale,rad=geometry(h,med.shape)
        tables=list(csv.DictReader((out/(visit+'_peaks.csv')).open()));allpeaks[visit]=tables
        bright=max(tables,key=lambda r:float(r['snr']))
        p,mod,corr=fit_primary(med,psf,int(bright['x0']),int(bright['y0']))
        model_amplitude=matched(np.where(np.isfinite(med),mod,np.nan),kernel)
        cleaned=matched(med-mod,kernel);ccut,crank,_=thresholds(cleaned,rad,cen,scale)
        fits_out.append(dict(visit=visit,fit_x=p[0],fit_y=p[1],east=(mat@(p[:2]-cen))[0],north=(mat@(p[:2]-cen))[1],
                             fit_flux_uJy=p[2]*h['PIXAR_SR']*1e12,core_r=corr))
        def read_short(n):
            with fits.open(cache/n.replace('nrcalong','nrca2')) as f:
                assert f[0].header['APERNAME']=='NRCA2_MASK335R'
                assert f[0].header['FILTER'] in ('F200W','F210M')
                img=f['SCI'].data.astype('f8');img[(f['DQ'].data&3)!=0]=np.nan
                return img,f['SCI'].header.copy(),f[0].header['FILTER']
        sw,h2,band=read_short(info['science_filename']);mat2,cen2,scale2,r2=geometry(h2,sw.shape[1:])
        sr=[read_short(k)[0] for k in info['reference_filenames']]
        sim=median(reduce_cube(sw,sr,r2,sky_bounds=(3.9,4.5)),axis=0)
        for row in tables:
            x,y=int(row['x0']),int(row['y0'])
            sky=np.array([float(row[k]) for k in ['east_arcsec','north_arcsec']]);xy=cen2+np.linalg.solve(mat2,sky)
            radius=np.linalg.norm(xy-cen2);phi=np.linspace(0,2*np.pi,int(2*np.pi*radius/8),endpoint=False)
            pts=(cen2+radius*np.array([np.cos(phi),np.sin(phi)]).T)
            controls=np.array([aperture(sim,*pt) for pt in pts if np.linalg.norm(pt-xy)>12])
            flux=(aperture(sim,*xy)-np.nanmedian(controls))*h2['PIXAR_SR']*1e12
            noise=np.nanstd(controls,ddof=1)*h2['PIXAR_SR']*1e12
            pars,model,r=fit_primary(med,psf,x,y)
            phot.append(dict(visit=visit,x=x,y=y,filter=band,psf_core_r=r,peak_modelsub_SNR=crank[y,x],
                             peak_modelsub_above_cut=int(cleaned[y,x]>ccut[y,x]),
                             calibrated_bright_source_model_fraction=model_amplitude[y,x]/amp[y,x],
                             dq_bitwise_or_at_peak=int(np.bitwise_or.reduce(fits.getdata(cache/info['science_filename'],'DQ')[:,y,x])),
                             short_aperture_uJy=flux,short_rms_uJy=noise,short_coverage=int(np.isfinite(flux)),
                             short_nctrl=int(np.isfinite(controls).sum()),long_flux_uJy=row['flux_uJy']))
        xyz=np.array([h['JWST_'+c] for c in 'XYZ'])/149597870.7
        frames[visit]=(mat,cen,amp,cut,snr,info['expstart_mjd'],xyz,h['PIXAR_SR']*1e12)
    xrows=[]
    for a,rows in allpeaks.items():
        b=next(k for k in VISITS if k!=a and k.split('_')[0]==a.split('_')[0])
        host=VISITS[a][-1];host_id={'HD38206':'2967287596507432704','HD53143':'5479222240596469632'}[host]
        cats=list(csv.DictReader((out/(host+'-gaia.csv')).open()));h=next(g for g in cats if g['source_id']==host_id)
        ra,de=np.deg2rad([float(h['ra']),float(h['dec'])])
        ax=np.array([[-np.sin(ra),np.cos(ra),0],[-np.cos(ra)*np.sin(de),-np.sin(ra)*np.sin(de),np.cos(de)]])
        a_mat,a_cen,_,__,___,a_mjd,a_xyz,____=frames[a]
        mat,cen,amp,cut,snr,mjd,xyz,unit=frames[b]
        mu=np.array([float(h['pmra']),float(h['pmdec'])]);par=float(h['parallax'])
        bg=(-mu*(mjd-a_mjd)/365.25+par*ax@(xyz-a_xyz))/1000.
        for row in rows:
            x,y=int(row['x0']),int(row['y0']);sky=a_mat@([x,y]-a_cen)
            approx_host_sky=(mu*(2000+(a_mjd-51544.5)/365.25-2016)-par*ax@a_xyz)/1000.
            neighbors=[]
            for g in cats:
                if g['source_id']==host_id:continue
                base=np.array([(float(g['ra'])-np.rad2deg(ra))*np.cos(de),float(g['dec'])-np.rad2deg(de)])*3600
                gm=np.array([float(g[k] or 0) for k in ['pmra','pmdec']])
                pos=base+(gm*(2000+(a_mjd-51544.5)/365.25-2016)-float(g['parallax'] or 0)*ax@a_xyz)/1000
                neighbors.append((np.linalg.norm(pos-approx_host_sky-sky),g['source_id']))
            near=min(neighbors) if neighbors else (np.nan,'')
            for label,pix in [('comoving',cen+np.linalg.solve(mat,sky)),('background',cen+np.linalg.solve(mat,sky+bg)),('detector',[x,y])]:
                xx,yy=np.rint(pix).astype(int)
                patch=snr[yy-1:yy+2,xx-1:xx+2]
                dist=min([np.linalg.norm(pix-np.array([int(e['x0']),int(e['y0'])])) for e in allpeaks[b]],default=np.nan)
                xrows.append(dict(origin=a,x=x,y=y,other=b,hypothesis=label,pred_x=pix[0],pred_y=pix[1],
                                  origin_snr=row['snr'],forced_snr=snr[yy,xx],max_3x3_snr=np.nanmax(patch),
                                  nearest_crossing_px=dist,forced_flux_uJy=amp[yy,xx]*unit,
                                  local_cut_uJy=cut[yy,xx]*unit,bg_east_shift_mas=1000*bg[0],bg_north_shift_mas=1000*bg[1],
                                  nearest_gaia_other_id=near[1],nearest_gaia_separation_arcsec=near[0]))
    for fn,rows in [('paired_rolls.csv',xrows),('paired_color_morphology.csv',phot),('brightest_psf_fits.csv',fits_out)]:write_table(out/fn,list(rows[0]),rows)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('cache',type=Path);p.add_argument('out',type=Path)
    a=p.parse_args();run(a.parent,a.out,a.cache)
