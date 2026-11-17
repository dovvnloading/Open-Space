"""Cross-check a frozen sample against SIMBAD, KYSO and a published disk census."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode

import numpy as np
import requests
from astropy.io.votable import parse_single_table
from astropy.coordinates import SkyCoord
import astropy.units as u

OUT=Path(__file__).resolve().parent/'classifications'
RAW=OUT/'raw'
SIM=['https://simbad.cds.unistra.fr/simbad/sim-tap/sync',
     'https://simbad.cfa.harvard.edu/simbad/sim-tap/sync','https://simbad.idia.ac.za/simbad/sim-tap/sync']
VIZ=['https://vizier.idia.ac.za/viz-bin/votable','https://vizier.cfa.harvard.edu/viz-bin/votable',
     'https://vizier.cds.unistra.fr/viz-bin/votable']
GAIA=['https://gea.esac.esa.int/tap-server/tap/sync']
SECURE_YOUNG={'Y*O','TT*','Ae*','Or*'}
SECURE_GAL={'G','AGN','SyG','Sy1','Sy2','QSO','Bla','BLL','LIN','EmG','rG','SBG','LSB','bCG','H2G','GiP','GiG','GiC','BiC'}
SECURE_MULTIPLE={'**','SB*','EB*','El*','RS*','Sy*','XB*','LXB','HXB','CV*'}
POSSIBLE={'Y*?','TT?','Ae?','G?','AGN?','AG?','Q?','Bla?','BLL?','**?','EB?','SB?','Sy?'}
receipts={}
errors={}


def fetch(name, hosts, params, parser):
    path=RAW/(name+'.gz')
    if path.exists():
        data=gzip.decompress(path.read_bytes())
        if hashlib.sha256(data).hexdigest()!=receipts[name]['sha256']:
            raise ValueError(f'Changed replay input: {name}')
        return parser(data)
    failures=[]
    for host in hosts:
        try:
            url=host+'?'+urlencode(params)
            response=requests.get(url,timeout=(12,80),headers={'User-Agent':'Open-Space-astronomical-research/1.0'})
            response.raise_for_status()
            parsed=parser(response.content)
            receipts[name]=dict(url=response.url,utc=datetime.now(timezone.utc).isoformat(),
                                 http=response.status_code,bytes=len(response.content),
                                 sha256=hashlib.sha256(response.content).hexdigest())
            path.write_bytes(gzip.compress(response.content,mtime=0))
            print(name,len(parsed),'records',flush=True)
            return parsed
        except Exception as exc:
            failures.append(f'{host}: {str(exc)[:260]}')
    raise RuntimeError('; '.join(failures))


def adql(name,hosts,query,cols,maxrows=100000):
    def parse(data):
        reader=csv.DictReader(io.StringIO(data.decode('utf-8-sig')))
        if not set(cols)<=set(reader.fieldnames or []):
            raise ValueError('Missing expected fields: '+str((reader.fieldnames or [])[:9]))
        vals=list(reader)
        if len(vals)>=maxrows: raise ValueError('Possibly truncated TAP response')
        return vals
    return fetch(name,hosts,dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC=maxrows,QUERY=query),parse)


def viz(name,catalog,minrows):
    def parse(data):
        tab=parse_single_table(io.BytesIO(data)).to_table(use_names_over_ids=True)
        if not minrows<=len(tab)<30000: raise ValueError('Incomplete published table: '+str(len(tab)))
        return tab
    return fetch(name,VIZ,{'-source':catalog,'-out.max':'30000','-out':'**','-out.add':'_RAJ,_DEJ'},parse)


def cells(table,key):
    return [str(x).strip() if not np.ma.is_masked(x) else '' for x in table[key]]


def coords(table):
    names=table.colnames
    ra=next((n for n in ['_RAJ2000','_RAJ','RAJ2000','RAdeg','RA_ICRS','RA_ICRS_','_RA'] if n in names),None)
    dec=next((n for n in ['_DEJ2000','_DEJ','DEJ2000','DEdeg','DE_ICRS','DE_ICRS_','_DE'] if n in names),None)
    if not ra or not dec: raise ValueError('No coordinates: '+str(names))
    def as_deg(x,y):
        if not x or not y: return [np.nan,np.nan]
        try: return [float(x),float(y)]
        except ValueError:
            c=SkyCoord(x,y,unit=(u.hourangle,u.deg))
            return [c.ra.deg,c.dec.deg]
    return np.array([as_deg(x,y) for x,y in zip(cells(table,ra),cells(table,dec))])


def sky_distance(ra,dec,xy):
    lon=np.deg2rad(xy[:,0]-ra);lat=np.deg2rad(xy[:,1]-dec)
    h=np.sin(lat/2)**2+math.cos(math.radians(dec))*np.cos(np.deg2rad(xy[:,1]))*np.sin(lon/2)**2
    return np.rad2deg(2*np.arcsin(np.sqrt(np.clip(h,0,1))))*3600


def typed(types):
    key=set(types.split('|'))
    found=set()
    if key & SECURE_YOUNG: found.add('young star')
    if key & SECURE_GAL: found.add('extragalactic')
    if key & SECURE_MULTIPLE: found.add('multiple')
    return found,bool(key & POSSIBLE)


def main():
    p=argparse.ArgumentParser();p.add_argument('inputs',type=Path);args=p.parse_args()
    RAW.mkdir(parents=True,exist_ok=True)
    f=OUT/'receipts.json'
    if f.exists(): receipts.update(json.loads(f.read_text())['primary'])
    def read(name):
        return list(csv.DictReader((args.inputs/name).open()))
    allstar={r['source_id']:r for r in read('stars.csv')}
    ids=sorted(r['source_id'] for r in read('screen.csv') if r['status'] in ['screened','photometric excess'])
    dr1={r['source_id'] for r in read('dr1_screen.csv') if r['extraction']=='strict' and r['status'] in ['screened','trigger']}
    assert len(ids)==len(set(ids))==404 and len(dr1)==55 and dr1<=set(ids)
    def jobsafe(name,function):
        try: return function()
        except Exception as exc:
            errors[name]=str(exc)
            print('INCOMPLETE',name,str(exc),flush=True)
            return None
    def sim_field():
        where="1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',180,30,2.05))"
        q='SELECT b.oid,b.main_id,b.ra,b.dec,b.otype,a.otypes FROM basic AS b LEFT OUTER JOIN alltypes AS a ON a.oidref=b.oid WHERE '+where
        data=adql('simbad_field',SIM,q,['oid','main_id','ra','dec','otype','otypes'])
        count=adql('simbad_count',SIM,'SELECT COUNT(*) AS n FROM basic AS b WHERE '+where,['n'])
        if not data or len(data)!=int(count[0]['n']) or len({r['oid'] for r in data})!=len(data):raise ValueError('SIMBAD field/count mismatch')
        return data
    def exact():
        out=[]
        for i in range(0,len(ids),42):
            labels=[f"'{prefix} {sid}'" for sid in ids[i:i+42] for prefix in ['Gaia DR3','Gaia EDR3']]
            q='SELECT i.id,i.oidref FROM ident AS i WHERE i.id IN ('+','.join(labels)+')'
            out+=adql(f'simbad_ids_{i:03d}',SIM,q,['id','oidref'])
        return out
    def variable():
        q=("SELECT g.source_id,g.non_single_star,c.best_class_name FROM gaiadr3.gaia_source AS g "
           "LEFT OUTER JOIN gaiadr3.vari_classifier_result AS c ON g.source_id=c.source_id WHERE g.source_id IN (")
        out=[]
        for i in range(0,len(ids),80):
            out+=adql(f'gaia_class_{i:03d}',GAIA,q+','.join(ids[i:i+80])+')',['source_id','non_single_star','best_class_name'])
        if {r['source_id'] for r in out}!=set(ids): raise ValueError('Missing Gaia IDs')
        return out
    with ThreadPoolExecutor(5) as pool:
        tasks={k:pool.submit(jobsafe,k,fn) for k,fn in {
            'simbad_field':sim_field,'simbad_ident':exact,'gaia_variability':variable,
            'published_kyso':lambda:viz('kyso','J/A+A/674/A21/kyso',11671),
            'published_ir_prime':lambda:viz('cotten_song_table3','J/ApJS/225/15/table3',400),
            'published_ir_reserved':lambda:viz('cotten_song_table4','J/ApJS/225/15/table4',1200)}.items()}
        data={k:v.result() for k,v in tasks.items()}
    sim=data['simbad_field'] or []
    sim_by_id={s['oid']:s for s in sim}
    sim_xy=np.array([[float(r['ra']),float(r['dec'])] for r in sim]).reshape(-1,2)
    ident={sid:set() for sid in ids}
    for r in data['simbad_ident'] or []: ident[r['id'].split()[-1]].add(r['oidref'])
    gaia={r['source_id']:r for r in data['gaia_variability'] or []}
    kyso=data['published_kyso']
    young=set()
    young_xy=np.empty((0,2));young_ids=[]
    if kyso is not None:
        young_ids=cells(kyso,'GaiaDR3')
        young=set(young_ids)
        young_xy=coords(kyso)
    disk_catalogs={key:coords(data[key]) for key in ['published_ir_prime','published_ir_reserved'] if data[key] is not None}
    summary=[];nearby=[]
    for sid in ids:
        star=allstar[sid]
        ra,dec,pmra,pmde=map(float,[star[k] for k in ('ra','dec','pmra','pmdec')])
        ra0=ra-16*pmra/3.6e6/math.cos(math.radians(dec));dec0=dec-16*pmde/3.6e6
        ra10=ra-6*pmra/3.6e6/math.cos(math.radians(dec));dec10=dec-6*pmde/3.6e6
        known=set();uncertain=set();names=[];types=[]
        if sid in young: known.add('KYSO exact published Gaia ID')
        if len(young_xy) and any((not young_ids[k]) for k in np.flatnonzero(sky_distance(ra0,dec0,young_xy)<=2)):
            uncertain.add('KYSO position only')
        # Literature coordinates may trace J2000 or WISE-epoch positions.
        for catalog,disk_xy in disk_catalogs.items():
            sep=np.minimum(sky_distance(ra0,dec0,disk_xy),sky_distance(ra10,dec10,disk_xy))
            for i in np.flatnonzero(sep<=3):
                uncertain.add('published IR-excess coordinate match: '+catalog)
                nearby.append(dict(source_id=sid,catalog=catalog,object_id=str(i+1),
                                  distance_arcsec=round(float(sep[i]),3),types='IR-excess list'))
        gr=gaia.get(sid,{})
        if int(gr.get('non_single_star') or 0)>0:known.add('Gaia non-single-star solution')
        if gr.get('best_class_name') in ['YSO','AGN','QSO','GALAXY']:
            uncertain.add('Gaia photometric '+gr['best_class_name'])
        listed=ident[sid]
        if len(listed)>1: uncertain.add('multiple SIMBAD identities for ID')
        for oid in listed:
            if oid not in sim_by_id:
                uncertain.add('SIMBAD ID outside retrieved footprint or missing position')
                continue
            s=sim_by_id[oid];names.append(s['main_id'])
            typ=s['otypes']+'|'+s['otype'];types.append(typ)
            ph,potential=typed(typ)
            if len(listed)==1: known.update('SIMBAD '+p for p in ph)
            if potential: uncertain.add('SIMBAD candidate classification')
        d=sky_distance(ra0,dec0,sim_xy) if len(sim_xy) else []
        close=0
        for ii in np.flatnonzero(np.asarray(d)<=12):
            s=sim[ii];t=s['otypes']+'|'+s['otype']
            nearby.append(dict(source_id=sid,catalog='SIMBAD',object_id=s['main_id'],
                               distance_arcsec=round(float(d[ii]),3),types=t))
            if s['oid'] in listed: continue
            close+=1
            if listed or d[ii]>2: uncertain.add('separate or unresolved SIMBAD neighbor within 12 arcsec')
            else:
                ph,potential=typed(t)
                if ph or potential: uncertain.add('position-only SIMBAD natural-class match')
        status='known association' if known else 'ambiguous association' if uncertain else 'no listed association'
        summary.append(dict(source_id=sid,in_strict_dr1=int(sid in dr1),status=status,
             secure=' | '.join(sorted(known)),ambiguous=' | '.join(sorted(uncertain)),
             simbad_id_count=len(listed),simbad_names=' | '.join(names),simbad_types=' | '.join(types),
             other_simbad_within_12as=close,gaia_variability=gr.get('best_class_name',''),
             gaia_non_single_star=gr.get('non_single_star',''),known_catalog_ID=int(sid in young)))
    def write(name,records,columns):
        with (OUT/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(records)
    write('stars.csv',summary,list(summary[0]))
    write('neighbors.csv',nearby,['source_id','catalog','object_id','distance_arcsec','types'])
    def stats(ss):
        return dict(total=len(ss),secure=sum(bool(s['secure']) for s in ss),
            ambiguous_only=sum(s['status']=='ambiguous association' for s in ss),
            both_secure_and_ambiguous=sum(bool(s['secure']) and bool(s['ambiguous']) for s in ss),
            retain_after_secure=sum(not s['secure'] for s in ss),
            retain_after_secure_and_ambiguous=sum(not s['secure'] and not s['ambiguous'] for s in ss),
            simbad_exact_ID=sum(s['simbad_id_count']>0 for s in ss),
            no_simbad_exact_ID=sum(s['simbad_id_count']==0 for s in ss),
            secure_reasons=dict(Counter(why for s in ss for why in s['secure'].split(' | ') if why)),
            ambiguous_reasons=dict(Counter(why for s in ss for why in s['ambiguous'].split(' | ') if why)))
    final=dict(complete=not errors,errors=errors,wise=stats(summary),spherex=stats([s for s in summary if s['in_strict_dr1']]))
    inputs={n:hashlib.sha256((args.inputs/n).read_bytes()).hexdigest() for n in ['screen.csv','stars.csv','dr1_screen.csv']}
    (OUT/'counts.json').write_text(json.dumps(final,indent=2)+'\n')
    (OUT/'receipts.json').write_text(json.dumps(dict(primary=receipts,inputs=inputs,errors=errors),indent=2)+'\n')
    print(json.dumps(final,indent=2),flush=True)


if __name__=='__main__':
    main()
