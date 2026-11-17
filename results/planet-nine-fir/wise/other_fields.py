"""Geometric sample of additional underdetermined historic candidates for WISE."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'other'))
from candidates import circular_family,other_sky

def fields():
    return [f'{name}-h{heading:03d}' for name in ['Chen1','Chen2','RR'] for heading in range(0,360,45)]

def state(field):
    name,heading=field.split('-h');r,v,t=circular_family(name)
    radius_steps=7 if name=='RR' else 11
    idx=int(heading)//5*radius_steps+radius_steps//2
    return r[idx:idx+1],v[idx:idx+1],t

def position(field,t):
    return other_sky(state(field),t)
