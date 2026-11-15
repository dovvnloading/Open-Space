"""Independent eclipse localization on KIC 8110758, a separate Q10 target and detector."""
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from catalog_screen import HERE, fetch
from pixel_localization import (EB, TARGET, NEIGHBOR, event_differences, fit_image,
                                fixed_location, local_prf, propagated_field)

TPF = "kplr008110758-2011271113734_lpd-targ.fits.gz"
PRF = "kplr09.3_2011265_prf.fits"
TPF_URL = "https://stpubdata.s3.amazonaws.com/kepler/public/target_pixel_files/0081/008110758/" + TPF
PRF_URL = "https://archive.stsci.edu/missions/kepler/fpc/prf/" + PRF


def main():
    tpf = fetch(TPF, TPF_URL, "7c00cba45d0ead3b23d188605c79c3e74ed5a3023e173088ee4b664e3fb739be")
    prf = fetch(PRF, PRF_URL, "6d2dd7337687c0eb2aa7d0f41a511e6113d82349c1e5cf46b280cef36fd43543")
    with fits.open(tpf) as hdus:
        table = hdus[1].data
        good = (table['QUALITY']==0) & np.isfinite(table['TIME'])
        stamp, time = table['FLUX'][good].astype(float), table['TIME'][good]
        aperture, header = hdus[2].data, hdus[2].header
        wcs, physical = WCS(header), WCS(header,key='P')
    finite = np.all(np.isfinite(stamp),axis=0)
    yy,xx = np.indices(aperture.shape)
    x,y = xx[finite],yy[finite]
    flux = stamp[:,finite]
    mask = ((aperture&2)>0)[finite]
    kernel = local_prf(prf, np.asarray(physical.pixel_to_world_values(x.mean(),y.mean())))
    field = propagated_field(wcs,Time(2454833+time.mean(),format='jd',scale='tdb').jyear)
    coords = field.set_index('gaia_id')[['x_pix','y_pix']]
    field.to_csv(HERE/'q10_gaia_field.csv',index=False)
    scale=3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
    # Predict the two distinct events directly from Gaia's geometric model.
    # Gaia TimeRef=2234.0540609516934 days after JD 2455197.5; BKJD zero is JD 2454833.
    period=1/1.9302825156283367480
    epoch=(2455197.5-2454833)+2234.0540609516934
    rows,full_maps=[],{}
    for label,phase in [('primary',.7466577),('secondary',.2377278)]:
        dates,diffs,baseline=event_differences(time,flux,phase,period=period,epoch=epoch)
        for split,use in [('all',np.ones(len(dates),bool)),('first',dates<np.median(time)),('second',dates>=np.median(time))]:
            arr=diffs[use]
            values=arr.mean(axis=0)
            error=arr.std(axis=0,ddof=1)/np.sqrt(len(arr))
            pars,chi2=fit_image(values,error,x,y,kernel,[x[np.argmax(values)],y[np.argmax(values)]])
            row=dict(eclipse=label,subset=split,n_events=len(arr),prf_x=pars[0],prf_y=pars[1],
                     chi2_free=chi2,n_cadences=len(time),start_bkjd=min(time),end_bkjd=max(time),
                     aperture_loss_e_s=values[mask].sum(),
                     ap_depth_ppm=1e6*values[mask].sum()/baseline[use][:,mask].sum(axis=1).mean())
            for key,sid in [('koi7865',TARGET),('eb',EB),('koi7866',NEIGHBOR)]:
                at=coords.loc[sid].to_numpy()
                row[f'{key}_sep_arcsec']=np.linalg.norm(at-pars[:2])*scale
                row[f'{key}_fixed_chi2']=fixed_location(values,error,x,y,kernel,at)
            rows.append(row)
            if split=='all': full_maps[label]=values
    pd.DataFrame(rows).to_csv(HERE/'q10_localizations.csv',index=False)
    pd.DataFrame(dict(x=x,y=y,optimal_aperture=mask,median_flux_e_s=np.median(flux,axis=0),
                      primary_loss_e_s=full_maps['primary'],secondary_loss_e_s=full_maps['secondary'])).to_csv(HERE/'q10_pixel_maps.csv',index=False)
    print(pd.DataFrame(rows)[['eclipse','subset','n_events','eb_sep_arcsec','koi7865_sep_arcsec','koi7866_sep_arcsec','ap_depth_ppm']].round(3).to_string(index=False))


if __name__=='__main__':
    main()
