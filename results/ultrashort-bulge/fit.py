"""Fit a bank of one- and two-lens templates at the measured epochs of each alert.

Build the C++ library as described in README.md. Input streams share the lens
geometry and have independent baselines, source fluxes and native error scales.
"""
from pathlib import Path
import ctypes
import sys
import numpy as np
import VBMicrolensing as vbm

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / 'ogle-dark-lenses'))
import lens  # shared OGLE flux scale and point-lens implementation

HALF = 1024


def engine():
    vb = vbm.VBMicrolensing()
    vb.LoadESPLTable(str(Path(vbm.__file__).parent / 'data/ESPL.tbl'))
    vb.Tol = 0.0001
    vb.a1 = 0.0
    return vb


def templates():
    vb = engine()
    z = np.linspace(-np.log(33), np.log(33), 2 * HALF + 1)
    tau = np.sign(z) * np.expm1(abs(z))
    rows, labels = [], []
    for u0 in (0.08, 0.35, 0.8, 1.4):
        rows.append(lens.magnification(tau**2 + u0**2) - 1)
        labels.append(('ps', u0, 0., 0., 0., 0.))
    for rho in (0.25, 0.75, 1.25, 2.5, 5.0):
        for frac in (0.15, 0.85):
            u0 = frac * rho
            rows.append([vb.ESPLMag2(np.hypot(t, u0), rho) - 1 for t in tau])
            labels.append(('fs', u0, rho, 0., 0., 0.))
    # Small family of smooth and caustic-crossing binary controls; masses in
    # total-Einstein-radius units. Center of mass is the origin of coordinates.
    for s, q, angle in [(0.8,0.1,0.7),(1.2,0.25,0.4),(1.7,0.05,1.4),(2.0,0.3,2.5)]:
        u0, rho = 0.22, 0.04
        co, si = np.cos(angle), np.sin(angle)
        rows.append([vb.BinaryMag2(s, q, t*co-u0*si, t*si+u0*co, rho)-1 for t in tau])
        labels.append(('bin', u0, rho, s, q, angle))
    return np.ascontiguousarray(rows), labels


def setup(library_path=None):
    global bank, labels, native
    bank, labels = templates()
    native = ctypes.CDLL(str(library_path if library_path is not None else ROOT / 'grid.so')).evaluate
    f64 = np.ctypeslib.ndpointer(np.float64, flags='C_CONTIGUOUS')
    i32 = np.ctypeslib.ndpointer(np.int32, flags='C_CONTIGUOUS')
    integer = ctypes.c_int
    native.argtypes = [integer,f64,f64,f64,i32,integer,f64,integer,f64,
                       integer,f64,integer,integer,f64,f64]
    native.restype = None


def clean(t, f, e, g):
    t, f, e, g = map(np.asarray, (t, f, e, g))
    keep = (t > 2450000) & (t < 2470000) & np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e>0)
    ii = np.where(keep)[0]
    ii = ii[np.argsort(t[ii], kind='stable')]
    return t[ii], f[ii], e[ii], np.unique(g[ii], return_inverse=True)[1]


def reported_error_floor(e, g):
    """Protect difference-image fits against near-zero reported uncertainties."""
    e = e.copy()
    for stream in np.unique(g):
        mask = g == stream
        e[mask] = np.maximum(e[mask], 0.15*np.median(e[mask]))
    return e


def axes(t, f, e, g, t0, te, season):
    z = np.zeros(len(t))
    for stream in np.unique(g):
        use = g == stream
        med = np.median(f[use])
        mad = 1.4826*np.median(np.abs(f[use]-med))
        z[use] = (f[use]-med)/np.hypot(e[use], max(0., mad))
    # Data-driven peak with support in adjacent observations. Additionally use
    # peak nearest the catalog epoch, which can have stronger local coverage.
    good = np.where(season)[0]
    if not len(good):
        good = np.arange(len(t))
    zg = z[good]
    if len(good)>=3:
        peak = good[np.argmax(np.minimum(np.minimum(zg[:-2],zg[1:-1]),zg[2:]))+1]
    else:
        peak = good[np.argmax(zg)]
    epoch = t0 if np.isfinite(t0) else t[peak]
    centers = [epoch, t[peak]]
    near = good[abs(t[good]-epoch)<max(2.,min(abs(te),15.) if np.isfinite(te) else 2.)]
    if len(near):
        centers.append(t[near[np.argmax(z[near])]])
    for shift in (-0.16, -0.04, 0.04, 0.16):
        centers.append(epoch + shift)
    scales = [.03,.065,.14,.28,.45,.65,1.2,2.5,6.,15.,40.,100.,250.]
    if np.isfinite(te) and 0.01 < te < 2000:
        scales.extend([te*.5,te,te*2])
    return np.unique(centers), np.unique(scales)


def screen(t, f, e, g, t0=np.nan, te=np.nan, year=2025, difference=False, caps=None):
    t,f,e,g = clean(t,f,e,g)
    if difference:
        e=reported_error_floor(e,g)
    if len(t)<6:
        raise ValueError('fewer than six usable epochs')
    season_start = 2458849.5 + (np.datetime64(f'{year}-01-01')-np.datetime64('2020-01-01')).astype(int)
    centers,scales = axes(t,f,e,g,t0,te,(t>=season_start)&(t<season_start+366))
    ng = max(g)+1
    cap = np.array([1e30 if difference else max(1e-12,1.15*np.median(f[g==s])) for s in range(ng)])
    if caps is not None:
        cap=np.ascontiguousarray(caps,dtype="d")
        if len(cap)!=ng: raise ValueError("baseline cap alignment")
    output=np.empty((len(centers),len(scales),len(labels)))
    native(len(t),np.ascontiguousarray(t,dtype='d'),np.ascontiguousarray(f,dtype='d'),
           np.ascontiguousarray(1/e**2,dtype='d'),g.astype('i'),ng,cap,
           len(centers),centers,len(scales),scales,len(labels),HALF,bank,output)
    result={'n':len(t),'streams':ng,'first_hjd':t[0],'last_hjd':t[-1]}
    for family in ('ps','fs','bin'):
        inds=np.array([i for i,l in enumerate(labels) if l[0]==family])
        a,b,c=np.unravel_index(np.argmax(output[:,:,inds]),(len(centers),len(scales),len(inds)))
        m=inds[c]
        result.update({family+'_t0':centers[a],family+'_te':scales[b],family+'_gain':output[a,b,m],
                       family+'_u0':labels[m][1]})
        if family=='fs': result['fs_rho']=labels[m][2]
        if family=='bin': result.update(zip(('bin_s','bin_q','bin_angle'),labels[m][3:]))
    # Actual sampling near the best single-lens peak.
    key=max(('ps','fs'),key=lambda x:result[x+'_gain'])
    off=t-result[key+'_t0']
    for radius in (.125,.25,.5):
        result[f'n_within_{radius}d']=np.sum(abs(off)<=radius)
    result['n_before_0.5d']=np.sum((off<0)&(off>=-.5))
    result['n_after_0.5d']=np.sum((off>0)&(off<=.5))
    return result
