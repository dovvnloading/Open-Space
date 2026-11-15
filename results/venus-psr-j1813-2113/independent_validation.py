#!/usr/bin/env python3
"""Independent Moshier validation of the 6 February 2027 Venus pulsar event.

NO JPL ephemeris files: SWIEPH is forbidden; this pins the semi-analytic
Moshier code (fit historically to DE404), not another JPL integration.
Reads frozen DE441 samples in PR #21 only after computing an independent
30-minute prediction for new radio sites. Install pyswisseph==2.10.3.2.
"""
import csv
import json
import math
import sys
from datetime import datetime,timedelta,timezone
from pathlib import Path

import swisseph as swe

ROOT=Path(__file__).resolve().parent
START=datetime(2027,2,6,16,40,tzinfo=timezone.utc)
FLAGS=swe.FLG_MOSEPH|swe.FLG_EQUATORIAL|swe.FLG_J2000|swe.FLG_NOABERR|swe.FLG_NOGDEFL
AU=149597870.7
RADIUS_VENUS_KM=6051.8
# WGS84 approximate sites, positive east, m above ellipsoid. Reuse original
# eight coordinates plus four sites not tested in PR #21.
SITES={
 "Green_Bank":(-79.8398,38.4331,824),
 "Hobart":(147.43917,-42.805,43),
 "IAR":(-58.140289,-34.86593,30),
 "McMurdo":(166.67,-77.85,20),
 "South_Pole":(0,-90,2835),
 "VLA":(-107.6177,34.0784,2124),
 "Warkworth":(174.663,-36.434,100),
 "geocenter":None,
 "MeerKAT":(21.443803,-30.712925,1050),
 "Parkes":(148.26351,-32.99840,415),
 "ATCA":(149.55012,-30.31289,237),
 "ASKAP":(116.5,-26.7,350),
}
POSITIONS={
 "MeerKAT_2024":(15*(18+13/60+39.9057931/3600),-(21+13/60+10.22721/3600)),
 "Thesis_2023":(15*(18+13/60+39.87/3600),-(21+13/60+0.4/3600)),
 "Parkes_2002":(15*(18+13/60+39.87/3600),-(21+13/60)),
}

def jd(seconds):
    t=START+timedelta(seconds=float(seconds))
    return swe.julday(t.year,t.month,t.day,t.hour+t.minute/60+t.second/3600+t.microsecond/3.6e9)

def venus(seconds,site):
    if site is not None:swe.set_topo(*site)
    options=FLAGS|(swe.FLG_TOPOCTR if site is not None else 0)
    (ra,dec,dist,*unused),ret=swe.calc_ut(jd(seconds),swe.VENUS,options)
    assert ret & swe.FLG_MOSEPH
    return ra,dec,math.degrees(math.asin(RADIUS_VENUS_KM/(AU*dist)))*3600

def sep(seconds,site,position):
    ra,dec,r=venus(seconds,site)
    ra0,dec0=position
    d=math.radians(dec);d0=math.radians(dec0);dl=math.radians(ra-ra0)
    x=math.cos(d)*math.sin(dl)
    y=math.cos(d0)*math.sin(d)-math.sin(d0)*math.cos(d)*math.cos(dl)
    z=math.sin(d0)*math.sin(d)+math.cos(d0)*math.cos(d)*math.cos(dl)
    return math.hypot(x,y)/z*206264.806247,r

def root(f,a,b):
    sign=f(a)>=0
    for _ in range(36):
        c=(a+b)/2
        if (f(c)>=0)==sign:a=c
        else:b=c
    return (a+b)/2

def contacts(site,position):
    # Golden-section closest approach; then two bisection limb contacts.
    a,b=0.,1800.
    g=(math.sqrt(5)-1)/2
    c=b-g*(b-a);d=a+g*(b-a)
    while b-a>.001:
        if sep(c,site,position)[0]<sep(d,site,position)[0]:
            b,d=d,c
            c=b-g*(b-a)
        else:
            a,c=c,d
            d=a+g*(b-a)
    t=(a+b)/2
    impact,radius=sep(t,site,position)
    if impact>=radius:return t,impact,None,None
    f=lambda sec:sep(sec,site,position)[0]-sep(sec,site,position)[1]
    assert f(0)>0 and f(1800)>0
    return t,impact,root(f,0,t),root(f,t,1800)

def altitude(seconds,site,body):
    if site is None:return None
    swe.set_topo(*site)
    equatorial,ret=swe.calc_ut(jd(seconds),body,
       swe.FLG_MOSEPH|swe.FLG_EQUATORIAL|swe.FLG_TOPOCTR)
    assert ret & swe.FLG_MOSEPH
    return swe.azalt(jd(seconds),swe.EQU2HOR,site,0.,10.,equatorial)[1]

def report():
    results={}
    for name,site in SITES.items():
        entry={}
        for key,pos in POSITIONS.items():
            t,impact,begin,end=contacts(site,pos)
            entry[key]={"closest_s":round(t,2),"impact_arcsec":round(impact,3),
                "radius_arcsec":round(sep(t,site,pos)[1],3),
                "ingress_s":None if begin is None else round(begin,2),
                "egress_s":None if end is None else round(end,2)}
        if site is not None:
            t=entry["MeerKAT_2024"]["closest_s"]
            entry["venus_alt_deg"]=round(altitude(t,site,swe.VENUS),2)
            entry["sun_alt_deg"]=round(altitude(t,site,swe.SUN),2)
        results[name]=entry
    delta=[]
    # Frozen JPL DE441 track is held out until after the new model calculation.
    for row in csv.DictReader((ROOT/"horizons.csv").open(newline="")):
        ra,dec,r=venus(float(row["seconds"]),SITES[row["site"]])
        dr=(ra-float(row["ra_deg"]))*math.cos(math.radians(float(row["dec_deg"])))*3600
        dd=(dec-float(row["dec_deg"]))*3600
        delta.append(math.hypot(dr,dd))
    assert len(delta)==248 and max(delta)<.12,(len(delta),max(delta))
    return {"model":"Moshier DE404-fit, ICRS J2000 astrometric airless",
            "version":swe.version,"horizons_frozen_rows":len(delta),
            "max_difference_arcsec":round(max(delta),5),
            "mean_difference_arcsec":round(sum(delta)/len(delta),5),
            "sites":results}

def frozen_checks(o):
    v=o["sites"]
    assert v["Warkworth"]["MeerKAT_2024"]["ingress_s"] is not None
    assert v["Warkworth"]["Thesis_2023"]["ingress_s"] is None
    assert v["South_Pole"]["Thesis_2023"]["ingress_s"] is not None
    assert v["MeerKAT"]["venus_alt_deg"] < -20
    assert v["ASKAP"]["venus_alt_deg"] < -15
    assert v["Parkes"]["venus_alt_deg"] < 30.25
    assert v["ATCA"]["venus_alt_deg"] < 12
    assert v["Warkworth"]["venus_alt_deg"] > 25
    # Model independent to <few seconds; these are frozen broad
    # contact windows, not fictitious sub-second prediction precision.
    for site,first,last in (("Warkworth",608,613),("Parkes",590,597),
                             ("South_Pole",764,770),("MeerKAT",868,875)):
        x=v[site]["MeerKAT_2024"]
        assert first < x["ingress_s"] < last
    # 2024 Warkworth chord margin and 2023 Warkworth miss must be >2 arcseconds.
    x=v["Warkworth"]
    assert x["MeerKAT_2024"]["radius_arcsec"]-x["MeerKAT_2024"]["impact_arcsec"]>7
    assert x["Thesis_2023"]["impact_arcsec"]-x["Thesis_2023"]["radius_arcsec"]>2.5

if __name__=="__main__":
    o=report()
    frozen_checks(o)
    print(json.dumps(o,sort_keys=True,indent=2))
