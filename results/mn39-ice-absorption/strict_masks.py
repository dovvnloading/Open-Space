"""Test band depths on identical pixels with nominal flags in all relevant channels."""
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from analyze import BANDS, INPUT, TARGET, band_measurement, write_table


def main():
    cube = np.load(INPUT/'spherex.npz')
    wave = cube['wavelength'].astype(float)
    wcs = WCS(fits.Header.fromstring(cube['header'].item()))
    yy, xx = np.indices(cube['image'].shape[1:])
    output = []
    for name, coord in [('MN39',TARGET),('J161036.07-512153.8',(242.650318,-51.364950))]:
        x,y = wcs.all_world2pix([coord], 0)[0]
        dist = np.hypot(xx-x, yy-y)
        for label, band in BANDS.items():
            center,left,right = band
            indices = np.zeros(len(wave), dtype=bool)
            for lo,hi in [left,right]:
                indices |= (wave>=lo-.0005) & (wave<=hi+.0005)
            central_index = np.argmin(abs(wave-center))
            indices[central_index] = True
            # These are identical spatial masks for the flux and both continuum sides.
            nominal = (cube['flags'][indices].astype('int64') & ~(1<<21)) == 0
            common = (nominal & np.isfinite(cube['image'][indices]) & (cube['hits'][indices]>0)).all(axis=0)
            sky = (dist>=3) & (dist<=5) & common
            for radius in [1.5,2.0]:
                aperture = (dist<=radius) & common
                signal = np.full(len(wave), np.nan)
                if sky.any() and aperture.any():
                    for i in np.flatnonzero(indices):
                        plane = cube['image'][i]
                        signal[i] = np.sum(plane[aperture] - np.median(plane[sky]))
                ratio,_ = band_measurement(wave,signal,band)
                output.append(dict(source=name,feature=label,radius_pix=radius,
                                   complete_aperture_pixels=int((dist<=radius).sum()),
                                   nominal_aperture_pixels=int(aperture.sum()),nominal_sky_pixels=int(sky.sum()),
                                   signal_at_band_center=round(float(signal[central_index]),5),
                                   transmission='' if not np.isfinite(ratio) else round(float(ratio),5)))
    write_table('strict-mask.csv',output)

if __name__=='__main__':
    main()
