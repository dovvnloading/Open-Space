"""At both hypotheses, force apertures in the historical HST rolls and later full-frame NIRCam images.
Use Gaia-propagated apparent host centres + FITS WCS; quote optical/IR scores as diagnostics.
"""
import argparse,csv,hashlib,json,sys
from pathlib import Path
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from astropy.coordinates import get_body_barycentric
from scipy import ndimage as ndi

def main(parent,cache,out):
    sys.path.insert(0,str(parent.resolve()/'vetting'))
    from archival_images import test_locus
    tables=[];ground=[]
    for host,js,names in [('HD38206','hst_inputs.json',['hd38206_o141','hd38206_o150']),('HD53143','wide_inputs.json',['hd53143_o001','hd53143_o008'])]:
        sid={'HD53143':'5479222240596469632','HD38206':'2967287596507432704'}[host]
        astro=next(r for r in csv.DictReader((out/(host+'-gaia.csv')).open()) if r['source_id']==sid)
        ra,de=np.deg2rad([float(astro['ra']),float(astro['dec'])]);ax=np.array([[-np.sin(ra),np.cos(ra),0],[-np.sin(de)*np.cos(ra),-np.sin(de)*np.sin(ra),np.cos(de)]])
        mu=np.array([float(astro['pmra']),float(astro['pmdec'])]);par=float(astro['parallax'])
        def project(t,xyz):return (mu*(t.jyear-2016)-par*ax@xyz)/1000
        for j in json.loads((out/js).read_text()):
            if j['target_name'] not in ['HD-38206','HD-53143']:continue
            file=j['dataURL'].split('/')[-1];path=cache/file
            with path.open('rb') as f:assert hashlib.file_digest(f,'sha256').hexdigest()==j['sha256']
            with fits.open(path) as f:
                im=f['SCI'].data.astype(float);im[f['WHT'].data<=0]=np.nan
                h=f['SCI'].header.copy();w=WCS(h,fix=False);p=f[0].header
                t=Time(h.get('MJD-AVG',.5*(j['t_min']+j['t_max'])),format='mjd')
                obs_xyz=np.array([h['JWST_'+c] for c in 'XYZ'])/149597870.7 if host=='HD53143' else get_body_barycentric('earth',t).xyz.to_value('au')
                host_offset=project(t,obs_xyz);scale=np.sqrt(abs(np.linalg.det(w.pixel_scale_matrix)))*3600
                def sky(off):
                    return tuple(np.rad2deg([ra,de])+off/np.array([3600*np.cos(de),3600]))
                hostpx=w.world_to_pixel_values(*sky(host_offset))
                for visit in names:
                    m=json.loads((out/(visit+'_measurement.json')).read_text())
                    h0=fits.getheader(cache/m['science_filename'],1);t0=Time(m['expstart_mjd'],format='mjd')
                    xyz0=np.array([h0['JWST_'+c] for c in 'XYZ'])/149597870.7
                    d_host=host_offset-project(t0,xyz0)
                    for row in csv.DictReader((out/(visit+'_peaks.csv')).open()):
                        off=np.array([float(row['east_arcsec']),float(row['north_arcsec'])])
                        for label,xy in [('co-moving',host_offset+off),('stationary',host_offset+off-d_host)]:
                            coord=sky(xy);px=w.world_to_pixel_values(*coord)
                            radius=.10 if host=='HD38206' else .16
                            rec=test_locus(im,w,px,hostpx,radius,scale)
                            tables.append(dict(origin=visit,x0=row['x0'],y0=row['y0'],dataset=file,epoch=t.utc.isot,filter=j['filters'],hypothesis=label,ra=coord[0],dec=coord[1],aperture_radius_arcsec=radius,**rec))
                    if host=='HD53143':
                        epoch=Time(2009.03,format='jyear');start=4.479*np.array([np.sin(np.deg2rad(286.2)),np.cos(np.deg2rad(286.2))])
                        bg=start-project(t0,xyz0)+project(epoch,get_body_barycentric('earth',epoch).xyz.to_value('au'))
                        fit=next(r for r in csv.DictReader((out/'brightest_psf_fits.csv').open()) if r['visit']==visit)
                        meas=np.array([float(fit[k]) for k in ('east','north')])
                        ground.append(dict(visit=visit,nici_epoch=2009.03,reference='Wahhaj et al. 2013, Table 7, candidate 1',published_rho=4.479,published_pa=286.2,pred_bg_east=bg[0],pred_bg_north=bg[1],measured_east=meas[0],measured_north=meas[1],residual_bg=np.linalg.norm(meas-bg),residual_comoving=np.linalg.norm(meas-start)))
    def write(fn,rows):
        with (out/fn).open('w') as f:
            tab=csv.DictWriter(f,list(dict.fromkeys(k for r in rows for k in r)));tab.writeheader();tab.writerows(rows)
    write('HST_and_wide_loci.csv',tables);write('NICI_track.csv',list({r['visit']:r for r in ground}.values()))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('parent',type=Path);p.add_argument('cache',type=Path);p.add_argument('out',type=Path);a=p.parse_args();main(a.parent,a.cache,a.out)
