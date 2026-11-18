"""Ellipsoidal ray geometry for synthetic/calibration use, units km.

Offset-axis shells approximate, but are not identical to, constant geodetic height.
No 2026 pointing or quality calibration is implied by these functions.
"""
import math
import numpy as np
A=6378.137
B=A*(1-1/298.257223563)

def surface(lon,lat,a,b):
    lon,lat=np.radians([lon,lat]); n=a*a/np.hypot(a*np.cos(lat),b*np.sin(lat))
    return np.array([n*np.cos(lat)*np.cos(lon),n*np.cos(lat)*np.sin(lon),n*(b/a)**2*np.sin(lat)])

def intersect(origin,point,a,b):
    origin=np.asarray(origin,dtype=float); point=np.asarray(point,dtype=float)
    if origin.shape!=(3,) or point.shape!=(3,) or not np.isfinite(origin).all() or not np.isfinite(point).all() or a<=0 or b<=0: raise ValueError('bad geometry')
    d=point-origin; norm=np.linalg.norm(d)
    if norm==0: raise ValueError('degenerate ray')
    d=d/norm; s=np.array([a,a,b]); x=origin/s; v=d/s
    aa=np.dot(v,v); bb=2*np.dot(x,v); cc=np.dot(x,x)-1; disc=bb*bb-4*aa*cc
    if disc<0: raise ValueError('ray misses shell')
    roots=[t for t in ((-bb-np.sqrt(disc))/(2*aa),(-bb+np.sqrt(disc))/(2*aa)) if t>0]
    if not roots: raise ValueError('shell behind observer')
    return origin+min(roots)*d

def triangulate(station1,coord1,station2,coord2):
    v1=surface(*coord1,A+14,B+6)-station1; v2=surface(*coord2,A+14,B+6)-station2
    matrix=np.column_stack((v1,-v2)); condition=np.linalg.cond(matrix)
    if not np.isfinite(condition) or condition>1e4: raise ValueError('ill-conditioned rays')
    t=np.linalg.lstsq(matrix,station2-station1,rcond=None)[0]
    if np.any(t<=0): raise ValueError('backward ray')
    p=station1+t[0]*v1; q=station2+t[1]*v2
    return (p+q)/2,float(np.linalg.norm(p-q))

def altitude_family(station,lon,lat,heights):
    if not heights or any(not np.isfinite(h) or not 40<=h<=120 for h in heights): raise ValueError('invalid height family')
    point=surface(lon,lat,A+14,B+6)
    return [intersect(station,point,A+h,B+h) for h in heights]
