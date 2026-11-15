"""Reproduce the July 13 source audit from the shipped extracts (Python 3.10+)."""
from collections import defaultdict
import csv
from datetime import datetime, timedelta
import json
import math
from pathlib import Path
import re
from statistics import mean, median

ROOT = Path(__file__).resolve().parent
IN = ROOT / 'inputs'
D = datetime.fromisoformat

def csv_rows(name):
    with (IN / name).open() as handle:
        return list(csv.DictReader(handle))

def utc(value):
    return value.isoformat(timespec='seconds')

def angular_distance(a, b):
    return abs((a - b + 180) % 360 - 180)

def numeric_class(text):
    return float(text[1:]) * {'A':1e-8,'B':1e-7,'C':1e-6,'M':1e-5,'X':1e-4}[text[0]]

def coordinates(longitude, latitude, observer_lon, observer_lat):
    difference = math.radians(longitude - observer_lon)
    b, b0 = map(math.radians, (latitude,observer_lat))
    x = math.cos(b) * math.sin(difference)
    y = math.sin(b) * math.cos(b0) - math.cos(b) * math.cos(difference) * math.sin(b0)
    z = math.sin(b) * math.sin(b0) + math.cos(b) * math.cos(difference) * math.cos(b0)
    return x, y, z

catalog = csv_rows('stix_extract.csv')
radio = [s.split() for s in (IN/'backside_radio.txt').read_text().splitlines() if s.strip()]
source = next(r for r in catalog if r['flare_id']=='2107131835')
july = [r for r in catalog if r['peak_UTC'].startswith('2021-07')]
header = json.loads((IN/'stix_headers.json').read_text())['SCI']
# The science file header is centered only ~25 min after the flare peak.
# Use the flare-epoch spacecraft Carrington longitude from its catalog coordinates.
observer_lon = (float(source['hgc_lon']) - float(source['hgs_lon']) + float(source['solo_position_lon'])) % 360
observer_lat = float(source['solo_position_lat'])
x, y = (float(source[k]) for k in ['hpc_x_solo','hpc_y_solo'])
observed_pa = math.degrees(math.atan2(-x,y)) % 360
rho = math.hypot(x,y) / header['RSUN_ARC']
alternatives = []
for line in (IN/'noaa_srs_20210722.txt').read_text().splitlines():
    match = re.match(r'(\d{4})\s+([NS])(\d+)[EW]\d+\s+(\d+)\b',line)
    if not match:
        continue
    region, hemisphere, lat, lon = match.groups()
    lat = int(lat) * (1 if hemisphere=='N' else -1)
    lon = int(lon)
    px,py,pz = coordinates(lon, lat, observer_lon, observer_lat)
    pa = math.degrees(math.atan2(-px,py)) % 360
    slip = []
    # Longitudinal drift over the eight-day interval and latitude-of-emission offsets.
    for dl in range(-5,6):
        for db in range(-3,4):
            sx,sy,sz = coordinates(lon+dl,lat+db,observer_lon,observer_lat)
            slip.append(angular_distance(observed_pa, math.degrees(math.atan2(-sx,sy))))
    alternatives.append({'noaa':10000+int(region),'longitude':lon,'latitude':lat,
        'predicted_pa_deg':round(pa,2),'pa_residual_deg':round(angular_distance(observed_pa,pa),2),
        'implied_height_Rsun':round(rho/math.hypot(px,py)-1,3),
        'pa_residual_range_drift_deg':[round(min(slip),2),round(max(slip),2)]})
alternatives.sort(key=lambda r:r['pa_residual_deg'])

constrained_longitudes = {}
# Enumerate radial heights and both allowed branches along the same image ray.
b0 = math.radians(observer_lat)
for h in [0.04,0.06,0.10,0.15]:
    allowed = []
    for i in range(1001):
        rr = rho + (1+h-rho)*i/1000
        X = (x/header['RSUN_ARC'])/rr
        Y = (y/header['RSUN_ARC'])/rr
        for sign in [-1,1]:
            Z = sign*math.sqrt(max(0,1-X*X-Y*Y))
            dl = math.degrees(math.atan2(X,Z*math.cos(b0)-Y*math.sin(b0)))
            allowed.append((observer_lon+dl)%360)
    constrained_longitudes[str(h)] = [round(min(allowed),2),round(max(allowed),2)]

rates = defaultdict(list)
for r in csv_rows('stix_rates.csv'):
    rates[r['product'],r['band_keV']].append((D(r['spacecraft_utc']),float(r['count_rate_s'])))
raw_peaks = {}
for (product, band), observations in rates.items():
    background = median(rate for tm,rate in observations if D('2021-07-13T17:40')<=tm<D('2021-07-13T18:00'))
    tm,rate = max(((tm,rate) for tm,rate in observations if D('2021-07-13T18:00')<=tm<D('2021-07-13T19:05')),key=lambda p:p[1])
    raw_peaks[product+' '+band] = {'spacecraft_utc':utc(tm),'earth_equivalent_utc':utc(tm+timedelta(seconds=header['EAR_TDEL'])),
        'rate_counts_s':round(rate,3),'baseline_counts_s':round(background,3),'peak_to_baseline':round(rate/background,2)}

heights = []
for row in (IN/'lasco_height_time.txt').read_text().splitlines():
    if not row.startswith('#'):
        rs,dt,tm,*extra = row.split()
        heights.append((D(dt.replace('/','-')+'T'+tm),float(rs)))
t0 = heights[0][0]
times = [(tm-t0).total_seconds() for tm,rs in heights]
rs = [r for tm,r in heights]
slope = sum((t-mean(times))*(r-mean(rs)) for t,r in zip(times,rs))/sum((t-mean(times))**2 for t in times)
intercept = mean(rs)-slope*mean(times)

# Test the same window about all July 13 clock phases (one-minute centers) and
# repeat over July. This is a descriptive occupancy, not an experiment-wide p-value.
jpeaks = [D(r['peak_UTC'])+timedelta(seconds=header['EAR_TDEL']) for r in july if r['hgs_lon'] and numeric_class(r['goes_estimated_mean_class'])>=1e-5]
def hits(center):
    return any(abs((t-center).total_seconds()) <= 30*60 for t in jpeaks)
control = {}
for label,first,days in [('same_day','2021-07-13',1),('month','2021-07-01',31)]:
    occupied = sum(hits(D(first)+timedelta(minutes=i)) for i in range(days*1440))
    control[label] = [occupied,days*1440]

# Report screen rows so misses, ambiguous regions and duplicate flare triggers are inspectable.
# One-AU approximation for selection only, with the exact FITS time correction used above for July 13.
screen = []
for row in radio:
    rt = D(row[0].replace('/','-')+'T'+row[1])
    candidates = []
    for r in catalog:
        if not r['hgs_lon']:
            continue
        peak = D(r['peak_UTC'])+timedelta(seconds=499.005*(1-float(r['solo_position_AU_distance'])))
        lag = (peak-rt).total_seconds()/60.
        if -35<=lag<=25 and numeric_class(r['goes_estimated_mean_class'])>=1e-5:
            candidates.append([r['flare_id'],round(lag,1),r['goes_estimated_mean_class'],round(float(r['hgs_lon']),1),round(float(r['hgs_lat']),1)])
    screen.append({'radio_utc':utc(rt),'catalog_location':row[6],'bright_flare_candidates':candidates})

out = {'radio_screen':screen,'stix_catalog_july_rows':len(july),'stix_catalog_localized_M_or_better_july':len(jpeaks),
    'july13_catalog_M_or_better':sum(t.day==13 for t in jpeaks),'flare_pa_deg':round(observed_pa,2),'flare_radius_Rsun':round(rho,4),
    'screen_spacecraft_carrington_lon':round(observer_lon,3),'raw_fits_peaks':raw_peaks,
    'lasco_linear_fit':{'number_of_points':len(heights),'km_s':round(slope*695700,1),
      'crosses_one_solar_radius_utc':utc(t0+timedelta(seconds=(1-intercept)/slope))},
    'noaa_region_controls':alternatives,'allowed_carrington_lon_deg_by_max_height_Rsun':constrained_longitudes,
    'catalog_timing_occupancy':control}
(ROOT/'measurements.json').write_text(json.dumps(out,indent=2)+'\n')
