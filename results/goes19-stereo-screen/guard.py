"""No target reader exists in this metadata-stage packet."""
import datetime as dt

def gate(phase, satellite, timestamps, qualities):
    if phase != 'synthetic': raise ValueError('target arrays sealed pending reviewed protocol')
    if satellite not in (18,19): raise ValueError('unknown spacecraft')
    if len(timestamps)!=len(qualities) or not timestamps: raise ValueError('invalid cohort')
    if any(not isinstance(t,(int,float)) or not 0<=t<120 for t in timestamps): raise ValueError('outside frozen interval')
    if any(b<a for a,b in zip(timestamps,timestamps[1:])): raise ValueError('time order')
    if any(q not in (0,1) for q in qualities): raise ValueError('unknown quality encoding')
    return [t for t,q in zip(timestamps,qualities) if q==0]

def predict_ray(observer,point,height):
    """Synthetic spherical Earth test only; NOT operational GLM ellipsoid navigation."""
    import math
    if height not in (60.,80.,100.): raise ValueError('height outside frozen family')
    if len(observer)!=3 or len(point)!=3 or not all(math.isfinite(v) for v in (*observer,*point)): raise ValueError('bad coordinates')
    norm=math.sqrt(sum((a-b)**2 for a,b in zip(point,observer)))
    if norm==0: raise ValueError('degenerate ray')
    d=[(b-a)/norm for a,b in zip(observer,point)]; dot=sum(a*b for a,b in zip(observer,d)); disc=dot*dot-(sum(a*a for a in observer)-(6371.+height)**2)
    if disc<0: raise ValueError('ray misses assumed shell')
    t=-dot-math.sqrt(disc)
    if t<=0: raise ValueError('shell behind observer')
    return tuple(a+t*b for a,b in zip(observer,d))
