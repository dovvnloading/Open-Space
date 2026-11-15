"""Search both common heliocentric directions of the osculating orbital planes."""
import numpy as np

from orbits import AU_KM, K, catalog, day


def node_times(x, y, q, e, tp):
    """Time at true anomaly atan2(y, x) on a bound or unbound orbit."""
    t = np.empty_like(x)
    ell = e < 1
    anomaly = 2*np.arctan2(np.sqrt(1-e[ell])*y[ell], np.sqrt(1+e[ell])*(1+x[ell]))
    t[ell] = (anomaly-e[ell]*np.sin(anomaly))*(q[ell]/(1-e[ell]))**1.5/K
    hyp = e > 1
    # This function is called only for finite, positive heliocentric radii.
    arg = np.sqrt((e[hyp]-1)/(e[hyp]+1))*y[hyp]/(1+x[hyp])
    anomaly = 2*np.arctanh(arg)
    t[hyp] = (e[hyp]*np.sinh(anomaly)-anomaly)*(q[hyp]/(e[hyp]-1))**1.5/K
    par = e == 1
    d = y[par]/(1+x[par])
    t[par] = (d+d**3/3)*np.sqrt(2*q[par]**3)/K
    return tp + t


def run():
    names, tp, q, e, _, p, v = catalog()
    source = np.flatnonzero(names == "3I/ATLAS").item()
    k = np.where(np.arange(len(names)) != source)[0]
    normals = np.cross(p[k], v[k])
    counts = []
    for angle in range(360):
        a = np.deg2rad(angle)
        rotate = np.array([[np.cos(a),-np.sin(a),0], [np.sin(a),np.cos(a),0], [0,0,1]])
        u, w = rotate@p[source], rotate@v[source]
        rays = np.cross(np.cross(u,w),normals)
        rays /= np.linalg.norm(rays,axis=1)[:,None]
        matches = []
        for sign in (-1,1):
            d = sign*rays
            xa, ya = d@u, d@w
            xb, yb = np.sum(p[k]*d,axis=1), np.sum(v[k]*d,axis=1)
            ra = q[source]*(1+e[source])/(1+e[source]*xa)
            rb = q[k]*(1+e[k])/(1+e[k]*xb)
            sel = (ra>.1) & (rb>ra) & (rb<4) & (rb-ra<1.5)
            idx = k[sel]
            src = np.full(len(idx), source)
            ta = node_times(xa[sel],ya[sel],q[src],e[src],tp[src])
            tb = node_times(xb[sel],yb[sel],q[idx],e[idx],tp[idx])
            ell = e[idx] < 1
            per = 2*np.pi*(q[idx][ell]/(1-e[idx][ell]))**1.5/K
            tb[ell] += np.ceil((ta[ell]-tb[ell])/per)*per
            speed = (rb[sel]-ra[sel])*AU_KM/(86400*(tb-ta))
            ok = (speed>=250)&(speed<=850)&(ta>=day("2010-01-01"))&(tb<day("2046-01-01"))
            matches.extend(zip(names[idx[ok]], speed[ok]))
        if angle == 0:
            print(f"3I with {len(k)} other MPC rows:", [(str(n), round(float(s),3)) for n,s in matches])
        counts.append(len(matches))
    labels, n = np.unique(counts[1:],return_counts=True)
    print("Control hit histogram (number of hits, number of rotations):", list(zip(labels.tolist(),n.tolist())))
    return counts


if __name__ == "__main__":
    run()
