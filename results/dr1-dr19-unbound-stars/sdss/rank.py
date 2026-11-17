"""Rank the loose heliocentric screen using point estimates; subsequent checks are independent."""
import csv
import math
from pathlib import Path
from fetch_survey import HERE, WORK, TABLES
# Standard ICRS -> Galactic rotation, rows l=0, l=90, north.
ROT = [(-.0548755604,-.8734370902,-.4838350155),
       (.4941094279,-.4448296300,.7469822445),
       (-.8676661490,-.1980763734,.4559837762)]
SOLAR = [11.1,248.5,7.25] # screening approximation in Galactic XYZ, km/s
def speed(r, distance):
    a,d=map(math.radians,(r['ra'],r['dec']))
    rv=r['v_rad']
    er=[math.cos(a)*math.cos(d),math.sin(a)*math.cos(d),math.sin(d)]
    ea=[-math.sin(a),math.cos(a),0.]
    ed=[-math.cos(a)*math.sin(d),-math.sin(a)*math.sin(d),math.cos(d)]
    v=[rv*x+4.74047*distance*(r['pmra']*y+r['pmde']*z) for x,y,z in zip(er,ea,ed)]
    return math.sqrt(sum((sum(q*t for q,t in zip(row,v))+sun)**2 for row,sun in zip(ROT,SOLAR)))
def read_raw(table,label):
    with (WORK/f'{table}_{label}.csv').open() as f:
        next(f)
        for row in csv.DictReader(f): yield row
def rank(row):
    r={k:float(row[k]) for k in ['ra','dec','pmra','pmde','v_rad','plx','e_plx','r_lo_photogeo','r_med_photogeo','r_hi_photogeo']}
    l,b=map(math.radians,map(float,(row['l'],row['b'])))
    rvgsr=r['v_rad']+SOLAR[0]*math.cos(l)*math.cos(b)+SOLAR[1]*math.sin(l)*math.cos(b)+SOLAR[2]*math.sin(b)
    inv=pg=pg_lo=None
    if all(abs(r[k])<9998 for k in ['pmra','pmde']):
        if r['e_plx']>0 and r['plx']>4*r['e_plx']: inv=speed(r,1/r['plx'])
        if 0<r['r_lo_photogeo']<=r['r_med_photogeo']<=r['r_hi_photogeo']:
            pg=speed(r,r['r_med_photogeo']/1000)
            pg_lo=speed(r,r['r_lo_photogeo']/1000)
    return {'rv_gsr':round(rvgsr,2),'v_inv':round(inv,2) if inv is not None else '',
            'v_pg':round(pg,2) if pg is not None else '',
            'v_pg_lo':round(pg_lo,2) if pg_lo is not None else ''}
def selected(r):
    return (abs(r['rv_gsr'])>=525 or (r['v_inv'] or 0)>=550 or (r['v_pg_lo'] or 0)>=550)
def main():
    for label in ['prefilter','control']:
        rows=[]
        for table in TABLES:
            for row in read_raw(table,label):
                row.update(rank(row))
                if selected(row):
                    rows.append({'survey':table,**row})
        if not rows: continue
        keys=list(dict.fromkeys(k for r in rows for k in r))
        with (HERE/f'{label}_tail.csv').open('w') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
        print(label,len(rows),len(set(r['sdss_id'] for r in rows)),sum(abs(r['rv_gsr'])>=525 for r in rows))
if __name__=='__main__': main()
