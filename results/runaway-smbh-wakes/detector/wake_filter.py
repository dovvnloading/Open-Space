"""Compensated straight-line matched filters and independent profile measurements."""

import math

import numpy as np
from scipy import fft, ndimage

FWHM_SIGMA = 2 * np.sqrt(2 * np.log(2))


def background(image, valid, mesh=48):
    """Bilinear interpolation of robust mesh medians."""
    h, w = image.shape
    yy, xx = np.indices(image.shape, dtype=np.float32)
    ny, nx = math.ceil(h / mesh), math.ceil(w / mesh)
    # Equal-size bins keep the interpolation knots aligned near partial edges.
    cell = np.full((ny, nx), np.nan, dtype=np.float32)
    for iy, rows in enumerate(np.array_split(np.arange(h), ny)):
        for ix, cols in enumerate(np.array_split(np.arange(w), nx)):
            s = np.s_[rows[0]:rows[-1]+1, cols[0]:cols[-1]+1]
            values = image[s][valid[s]]
            if values.size:
                cell[iy, ix] = np.median(values)
    cell[~np.isfinite(cell)] = np.nanmedian(cell)
    return ndimage.map_coordinates(cell, [(yy+.5)*ny/h-.5, (xx+.5)*nx/w-.5],
                                   order=1, mode="nearest")


def axes(angle):
    return np.array([np.cos(angle), np.sin(angle)]), np.array([-np.sin(angle), np.cos(angle)])


def profile(data, valid, x, y, length, width, angle, pixel_scale, min_length):
    """Verify coverage and measure extent/FWHM independently of the filter dimensions."""
    u, v = axes(angle)
    along = np.arange(-math.ceil(.9*length), math.ceil(.9*length)+1, dtype=float)
    cross = np.arange(-math.ceil(4*width), math.ceil(4*width)+.25, .25)
    px = x + cross[:, None]*v[0] + along*u[0]
    py = y + cross[:, None]*v[1] + along*u[1]
    coverage = ndimage.map_coordinates(valid.astype(np.uint8, copy=False), [py, px],
                                       order=0, mode="constant", cval=0)
    data = ndimage.map_coordinates(data, [py, px], order=1, mode="constant", cval=0)
    core = np.abs(along) <= .5*length
    if coverage[:, core].mean() < .98 or coverage[np.abs(cross) <= width][:, core].min() == 0:
        return None
    sides = (np.abs(cross) >= 2.5*width) & (np.abs(cross) <= 4*width)
    data -= np.median(data[sides], axis=0)
    trace = data[np.abs(cross) <= .5*width].mean(axis=0)
    smooth = ndimage.gaussian_filter1d(trace, max(1.5, .06*length))
    level = np.median(smooth[np.abs(along) < .35*length])
    if level <= 0:
        return None
    # Require extended signal on both sides of the center, reject a central compact peak.
    parts = np.array_split(trace[core], 5)
    if sum(np.mean(p) > .15*level for p in parts) < 4:
        return None
    center = len(along) // 2
    above = (smooth > .35*level) & (coverage[np.abs(cross) <= width].min(axis=0) > 0)
    if not above[center]:
        return None
    labels, _ = ndimage.label(above)
    run = along[labels == labels[center]]
    measured_length = run[-1] - run[0]
    if measured_length < max(min_length / pixel_scale, .7*length):
        return None
    interior = (along >= run[0]+.12*measured_length) & (along <= run[-1]-.12*measured_length)
    transverse = np.mean(data[:, interior], axis=1)
    transverse = ndimage.gaussian_filter1d(transverse, 1)
    # A centered and resolved cross-section must fall below half-peak on both sides.
    near = np.flatnonzero(np.abs(cross) <= .7*width)
    peak = near[np.argmax(transverse[near])]
    if transverse[peak] <= 0:
        return None
    mask = transverse >= .5*transverse[peak]
    ids, _ = ndimage.label(mask)
    selected = np.flatnonzero(ids == ids[peak])
    i, j = selected[0], selected[-1]
    if i == 0 or j == len(cross)-1:
        return None
    half = .5*transverse[peak]
    left = np.interp(half, transverse[i-1:i+1], cross[i-1:i+1])
    right = np.interp(half, transverse[[j+1,j]], cross[[j+1,j]])
    measured_width = max(1., right-left)
    ratio = measured_length/measured_width
    if ratio <= 20:
        return None
    shift = cross[peak]*v
    p0 = np.array([x, y])+shift+u*run[0]
    p1 = np.array([x, y])+shift+u*run[-1]
    return dict(x=float((p0[0]+p1[0])/2), y=float((p0[1]+p1[1])/2),
                x0=float(p0[0]), y0=float(p0[1]), x1=float(p1[0]), y1=float(p1[1]),
                length_px=float(measured_length), fwhm_px=float(measured_width),
                length_arcsec=float(measured_length*pixel_scale), axis_ratio=float(ratio),
                pa_image_deg=float(np.degrees(angle)), profile_at_limit=bool(
                    run[0] == along[0] or run[-1] == along[-1]))


def compatible(a, b):
    """Match unoriented segments using angle, midpoint and length."""
    angle = np.deg2rad(a["pa_image_deg"])
    u, v = axes(angle)
    delta = np.array([a["x"]-b["x"], a["y"]-b["y"]])
    axis_diff = abs((a["pa_image_deg"]-b["pa_image_deg"]+90) % 180 - 90)
    length, width = a["length_px"], a["fwhm_px"]
    return bool(axis_diff <= min(12, np.degrees(np.arctan2(2*width+1, .5*length)))
                and abs(delta @ u) <= .35*length
                and abs(delta @ v) <= max(2.5, 1.5*width)
                and .55 <= b["length_px"]/length <= 1.8)


def detect(image, valid, lengths, widths, pixel_scale, threshold=6, min_length=3):
    if valid.sum() < 1024:
        return []
    bg = background(image, valid)
    data = np.where(valid, image-bg, 0).astype("f4")
    if data[valid].max() == data[valid].min():
        return []
    # Zero-padding prevents features on opposite boundaries from supporting a seed.
    pad = int(np.ceil(max(lengths)/2 + 4*max(widths)))
    h, w = data.shape
    shape = (fft.next_fast_len(h+2*pad), fft.next_fast_len(w+2*pad))
    spec = fft.rfft2(data, s=shape)
    kx = fft.rfftfreq(shape[1]).astype("f4")[None, :]
    ky = fft.fftfreq(shape[0]).astype("f4")[:, None]
    best = np.full((h,w), -np.inf, "f4")
    which = np.zeros((h,w), np.int32)
    bank = []
    # Near the edges the convolution is truncated; omit those pixels from the noise scale.
    margin = max(2, min(pad, (min(h,w)-32)//2))
    sample = np.s_[margin:h-margin:3, margin:w-margin:3]
    sample_valid = valid[sample]
    if sample_valid.sum() < 100:
        return []
    for length in lengths:
        for width in widths:
            ntheta = max(12, math.ceil(np.pi*length/(1.25*width)))
            for angle in np.arange(ntheta)*np.pi/ntheta:
                u, v = axes(angle)
                axial = kx*np.float32(u[0])+ky*np.float32(u[1])
                normal = kx*np.float32(v[0])+ky*np.float32(v[1])
                # Top-hat along the line, Gaussian across it, minus two parallel sidebands.
                kernel = (np.sinc(length*axial)
                          * np.exp(-2*np.pi**2*((width/FWHM_SIGMA*normal)**2+(.6*axial)**2))
                          * (1-np.cos(2*np.pi*(2.5*width)*normal))).astype("f4")
                response = fft.irfft2(spec*kernel, s=shape)[:h, :w]
                sample_values = response[sample][sample_valid]
                q25, median, q75 = np.percentile(sample_values, [25,50,75])
                noise = (q75-q25)/1.3489795
                if noise <= 0:
                    continue
                score = ((response-np.float32(median))/np.float32(noise))
                better = score > best
                best[better] = score[better]
                which[better] = len(bank)
                bank.append((length, width, angle))
    seeds = valid & (best >= threshold) & (best == ndimage.maximum_filter(best, size=7))
    ys, xs = np.nonzero(seeds)
    output = []
    valid8 = valid.astype(np.uint8)
    for idx in np.argsort(best[ys, xs])[::-1]:
        x, y = xs[idx], ys[idx]
        length, width, angle = bank[which[y,x]]
        row = profile(data, valid8, x, y, length, width, angle, pixel_scale, min_length)
        if row and not any(compatible(old, row) for old in output):
            row.update(score=float(best[y,x]), filter_length_px=float(length),
                       filter_fwhm_px=float(width))
            output.append(row)
    return output
