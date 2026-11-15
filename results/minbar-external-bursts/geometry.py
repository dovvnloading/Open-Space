"""Check Wind-to-Earth delays against the pointing and independent localizations."""
import csv
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def direction(ra, dec):
    a, d = math.radians(ra), math.radians(dec)
    return [math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)]


def projection(xyz, ra, dec):
    return sum(x*y for x,y in zip(xyz,direction(ra,dec)))/299792.458


def sep(ra,dec,a,d):
    dot=sum(x*y for x,y in zip(direction(ra,dec),direction(a,d)))
    return math.degrees(math.acos(min(1,max(-1,dot))))


def geometry():
    rows = list(csv.DictReader((ROOT/'data/wind_vectors.csv').open()))
    events = {r['seq']:r for r in csv.DictReader((ROOT/'data/events.csv').open())}
    out = {}
    for row in rows:
        seq = row['seq']
        v = [float(row[k]) for k in ('x_km','y_km','z_km')]
        event = events[seq]
        ra, dec = float(event['rxte_ra']),float(event['rxte_dec'])
        m = {'wind_light_time_s':round(math.sqrt(sum(x*x for x in v))/299792.458,3),
             'pointing_delay_s':round(projection(v,ra,dec),3)}
        if seq == '8279':
            m['ipn_center_delay_s']=round(projection(v,134.486,-31.972),3)
            m['ipn_offset_deg']=round(sep(ra,dec,134.486,-31.972),3)
        if seq == '3165':
            m['ipn_endpoint_delays_s']=[round(projection(v,a,d),3) for a,d in [(246,-42),(253,-33)]]
            m['distance_to_full_annulus_deg']=round(abs(sep(ra,dec,222.855,-17.387)-31.387),3)
        if seq == '3234':
            m['solar_direction_delay_s']=round(projection(v,218.56,-15.12),3)
        out[seq]=m
    return out
