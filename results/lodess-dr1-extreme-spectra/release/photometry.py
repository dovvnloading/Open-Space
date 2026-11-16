"""Forced point-source photometry on a circular-beam, Jy/beam radio image.

The flux is a fixed-position PSF amplitude with a local background. The error is
the scatter of identical flux measurements at a grid of adjacent blank positions,
so spatially correlated image noise is measured at the same filtering scale.
This is appropriate for an isolated unresolved source. A separate aperture estimate
and local residual diagnostics are provided for inspecting this assumption.
"""
import math

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales


def mad(x):
    return 1.4826022185 * np.nanmedian(np.abs(x-np.nanmedian(x)))


class RadioImage:
    def __init__(self, filename):
        self.hdus = fits.open(filename, memmap=True)
        header = self.hdus[0].header
        if header['BUNIT'].replace(' ', '').upper() != 'JY/BEAM':
            raise ValueError('Expected flux density in Jy/beam')
        self.data = np.squeeze(self.hdus[0].data)
        if self.data.ndim != 2:
            raise ValueError('Expected singleton frequency/Stokes dimensions')
        self.wcs = WCS(header).celestial
        scale = proj_plane_pixel_scales(self.wcs)
        if abs(header['BMAJ']/header['BMIN']-1) > .01 or abs(scale[0]/scale[1]-1) > .01:
            raise ValueError('This estimator requires a circular beam and square pixels')
        self.fwhm = header['BMAJ'] / float(np.mean(scale))
        self.sigma = self.fwhm/np.sqrt(8*np.log(2))
        self.beam_area = 2*np.pi*self.sigma**2
        self.pixel_arcsec = float(np.mean(scale))*3600

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, traceback):
        self.hdus.close()

    def _sample(self, x, y, halfsize):
        cx, cy = int(round(x)), int(round(y))
        h = int(math.ceil(halfsize))
        if cx-h < 0 or cy-h < 0 or cx+h >= self.data.shape[1] or cy+h >= self.data.shape[0]:
            return None
        a = self.data[cy-h:cy+h+1, cx-h:cx+h+1]
        yy, xx = np.ogrid[cy-h:cy+h+1, cx-h:cx+h+1]
        return a, (xx-x)**2+(yy-y)**2

    def _estimate(self, x, y):
        patch = self._sample(x, y, 2.5*self.fwhm)
        if patch is None:
            return None
        pixels, r2 = patch
        sky = (r2 >= (1.5*self.fwhm)**2) & (r2 <= (2.5*self.fwhm)**2)
        core = r2 <= (1.25*self.fwhm)**2
        finite = np.isfinite(pixels)
        if finite[sky].mean() < .95 or finite[core].mean() < .95:
            return None
        background = np.median(pixels[sky & finite])
        core &= finite
        instr = np.exp(-.5*r2[core]/self.sigma**2)
        flux = np.sum((pixels[core]-background)*instr)/np.sum(instr**2)
        aperture = r2 <= self.fwhm**2
        encircled = 1 - np.exp(-.5*self.fwhm**2/self.sigma**2)
        aperture_flux = np.nansum(pixels[aperture]-background)/(self.beam_area*encircled)
        return float(flux), float(background), float(aperture_flux)

    def measure(self, ra, dec):
        x, y = map(float, self.wcs.world_to_pixel_values(ra, dec))
        if not np.isfinite([x,y]).all():
            return None
        target = self._estimate(x, y)
        patch = self._sample(x, y, 12*self.fwhm)
        if target is None or patch is None:
            return None
        pixels, r2 = patch
        annulus = (r2 >= (3*self.fwhm)**2) & (r2 <= (12*self.fwhm)**2)
        if np.isfinite(pixels[annulus]).mean() < .9:
            return None
        noise = mad(pixels[annulus])
        controls = []
        for iy in range(-4, 5):
            for ix in range(-4, 5):
                if 2 <= np.hypot(ix,iy) <= 4:
                    v = self._estimate(x+ix*3*self.fwhm, y+iy*3*self.fwhm)
                    if v is not None:
                        controls.append(v[0])
        if len(controls) < 25:
            return None
        flux, bg, ap = target
        err = mad(np.asarray(controls))
        cy, cx = np.unravel_index(np.nanargmax(np.where(r2 <= self.fwhm**2, pixels, np.nan)), pixels.shape)
        offset = math.sqrt(r2[cy,cx])*self.pixel_arcsec
        core = r2 <= (1.25*self.fwhm)**2
        residual = pixels[core]-bg-flux*np.exp(-.5*r2[core]/self.sigma**2)
        return dict(psf_jy=flux, psf_noise_jy=float(err),
                    controls=len(controls), control_median_jy=float(np.median(controls)),
                    background_jy_beam=bg, pixel_noise_jy_beam=float(noise),
                    aperture_jy=ap, peak_jy_beam=float(pixels[cy,cx]-bg),
                    peak_offset_arcsec=float(offset), core_resid_rms_jy_beam=float(np.sqrt(np.nanmean(residual**2))),
                    local_high_jy_beam=float(np.nanmax(pixels[annulus])-bg),
                    local_low_jy_beam=float(np.nanmin(pixels[annulus])-bg))
