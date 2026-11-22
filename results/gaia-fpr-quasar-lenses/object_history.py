"""Independent publication checks on every unlisted two-filter imaging pass."""
import csv
import gzip
import hashlib
import io
import json
import math
import re
import sys
import subprocess
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests
from astropy.io.votable import parse

# Catalogue file populations published in the authors' CDS data descriptions.
TABLES = {
    'J/A+A/692/A154/golden': ('2024A&A...692A.154W', 1867, 'RAdeg', 'DEdeg'),
    'J/A+A/692/A154/superset': ('2024A&A...692A.154W', 5286, 'RAdeg', 'DEdeg'),
    'J/ApJS/264/4/table1': ('2023ApJS..264....4M', 3140, 'RA_ICRS', 'DE_ICRS'),
}
TOPIC = re.compile(r'lens|dual.{0,12}(agn|quasar)|binary.{0,12}(agn|quasar)|'
                   r'(agn|quasar).{0,12}pair|double.{0,8}(nuclei|quasar)|merg(e|ing)', re.I)


def get(url, params, dest):
    res = requests.get(url, params=params, timeout=(6, 30))
    res.raise_for_status()
    packed=gzip.compress(res.content,mtime=0)
    dest.write_bytes(packed)
    return res.content, dict(url=res.url,upstream_bytes=len(res.content),file=dest.name,
                             upstream_sha256=hashlib.sha256(res.content).hexdigest(),
                             sha256=hashlib.sha256(packed).hexdigest())


def votable(data):
    doc = ET.fromstring(data)
    status = [e for e in doc.iter() if e.tag.endswith('INFO') and e.get('name','').upper()=='QUERY_STATUS']
    for item in status:
        if item.get('value','').upper() in {'ERROR','OVERFLOW'}:
            raise ValueError('VO service status '+str(item.get('value'))+': '+str(item.text)[:160])
    output = []
    for table in parse(io.BytesIO(data), verify='ignore').iter_tables():
        cols = [f.name or f.ID for f in table.fields]
        arr = table.array
        def string(v):
            if getattr(v,'mask',False):
                return ''
            return v.decode('utf-8','replace').strip() if isinstance(v,bytes) else str(v).strip()
        rows = [dict(zip(cols,(string(x) for x in row))) for row in arr]
        output.append(dict(columns=cols, rows=rows))
    return output


def distance(ra,dec,other_ra,other_dec):
    dr = math.radians((ra-other_ra+180)%360-180)
    dd = math.radians(dec-other_dec)
    a = math.sin(dd/2)**2 + math.cos(math.radians(dec))*math.cos(math.radians(other_dec))*math.sin(dr/2)**2
    return 7200*math.degrees(math.asin(min(1,math.sqrt(a))))


def match_table(table, packets, out):
    bib, count, ra_key, dec_key = TABLES[table]
    result = dict(table=table, bibcode=bib, expected_rows=count, attempts=[], matches=[])
    endpoints=[('static','https://vizier.cfa.harvard.edu/ftp/cats/'),('vo','https://vizier.cfa.harvard.edu'),
               ('vo','https://vizier.cds.unistra.fr'),('vo','https://vizier.idia.ac.za')]
    for kind,host in endpoints:
        try:
            if kind=='static':
                prefix,file=table.rsplit('/',1)
                base=host+prefix+'/'
                metadata,receipt=get(base+'ReadMe',{},out/(table.replace('/','_')+'_ReadMe.gz'))
                result['attempts'].append(receipt)
                definition=metadata.decode('utf-8','replace').split('Description of file: '+file+'.dat',1)[1]
                definition=definition.split('Description of file:',1)[0]
                spans={}
                for label in ('RAdeg','DEdeg'):
                    span=re.search(r'^\s*(\d+)-\s*(\d+)\s+F\S+\s+deg\s+'+label+r'\s',definition,re.M)
                    if span is None:
                        raise ValueError('CDS dictionary missing coordinate '+label)
                    spans[label]=slice(int(span[1])-1,int(span[2]))
                for ext in ('.dat','.dat.gz'):
                    try:
                        raw,receipt=get(base+file+ext,{},out/(table.replace('/','_')+ext+'.gz'))
                        result['attempts'].append(receipt)
                        data=gzip.decompress(raw) if raw.startswith(bytes((31,139))) else raw
                        lines=data.decode('ascii').splitlines()
                        if len(lines)!=count:
                            raise ValueError('CDS file row count '+str(len(lines)))
                        break
                    except Exception as exc:
                        result['attempts'].append(dict(url=base+file+ext,error=str(exc)[:150]))
                else:
                    raise ValueError('Both CDS byte-table encodings unavailable')
                valid=[[{ra_key:line[spans['RAdeg']].strip(),dec_key:line[spans['DEdeg']].strip(),
                         'published_row':line,'row_number':i} for i,line in enumerate(lines,1)]]
            else:
                data, receipt = get(host+'/viz-bin/votable',
                                    {'-source':table, '-out.all':1, '-out.max':'unlimited'},
                                    out/(table.replace('/','_')+'.xml.gz'))
                result['attempts'].append(receipt)
                tables = votable(data)
                valid = [t['rows'] for t in tables if {ra_key,dec_key} <= set(t['columns'])]
                if len(valid)!=1 or len(valid[0])!=count:
                    raise ValueError('unexpected rows by returned table '+str([len(t['rows']) for t in tables]))
            matches = []
            for row in valid[0]:
                ra, dec = float(row[ra_key]),float(row[dec_key])
                if not (math.isfinite(ra) and 0<=ra<360 and -90<=dec<=90):
                    raise ValueError('invalid published coordinate')
                for packet in packets:
                    near = min(distance(ra,dec,c['ra'],c['dec']) for c in packet['components'])
                    if near<=8:
                        matches.append(dict(system_id=packet['system_id'],nearest_image_arcsec=near,row=row))
            result.update(status='complete', matches=matches, rows=count)
            break
        except Exception as exc:
            result['attempts'].append(dict(host=host,error=str(exc)[:200]))
    result.setdefault('status','unavailable')
    return result


def ned(packet,out):
    cen = packet['imaging']['lens_model']['center_radec_deg']
    components = packet['components']
    records = []
    # Every object in an 8" lens-centred cone, including alternate IDs or neighbours.
    sql = ("SELECT * FROM NEDTAP.objdir WHERE 1=CONTAINS(POINT('J2000',ra,dec),"
           f"CIRCLE('J2000',{cen[0]:.10f},{cen[1]:.10f},{8/3600:.12f}))")
    params = dict(LANG='ADQL',REQUEST='doQuery',FORMAT='csv',MAXREC=1000,QUERY=sql)
    audit = []
    for host in ('ned.ipac.caltech.edu','vo.ned.ipac.caltech.edu'):
        try:
            raw,receipt = get('https://'+host+'/tap/sync',params,out/'ned_cone.csv.gz')
            audit.append(receipt)
            reader = csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
            names = {str(k).lower():k for k in reader.fieldnames or []}
            nm = next(names[k] for k in ('prefname','objname','obj_name','object_name','name') if k in names)
            if 'ra' not in names or 'dec' not in names:
                raise ValueError('Missing object coordinates')
            rows = list(reader)
            if len(rows)>=1000:
                raise ValueError('Response may be truncated')
            for row in rows:
                ra,dec=float(row[names['ra']]),float(row[names['dec']])
                records.append(dict(name=row[nm],catalogue=row,
                    nearest_image_arcsec=min(distance(ra,dec,c['ra'],c['dec']) for c in components)))
            break
        except Exception as exc:
            records.clear()
            audit.append(dict(host=host,error=repr(exc)[:250]))
    else:
        return dict(status='unavailable',attempts=audit,objects=[])
    def details(service,name,idx):
        r = dict(service=service,object=name,status='unavailable',attempts=[])
        for host in ('ned.ipac.caltech.edu','vo.ned.ipac.caltech.edu'):
            try:
                data,rec = get('https://'+host+'/NED::API/'+service,
                               {'TARGET':name,'MAXREC':5000},out/f'ned_{idx}_{service}.xml.gz')
                r['attempts'].append(rec)
                tables=votable(data)
                if not tables or any(len(t['rows'])>=5000 for t in tables):
                    raise ValueError('missing tables or incomplete response')
                r.update(status='response_received', tables=tables)
                break
            except Exception as exc:
                r['attempts'].append(dict(host=host,error=str(exc)[:150]))
        return r
    with ThreadPoolExecutor(4) as pool:
        jobs = [(r,pool.submit(details,s,r['name'],i)) for i,r in enumerate(records) for s in
                ['CrossidsOfObject','ReferencesOfObject','NotesOfObject','ClassificationsOfObject','RedshiftsOfObject']]
        for r,f in jobs:
            r.setdefault('object_services',[]).append(f.result())
    return dict(status='complete_cone',objects=records,attempts=audit)


def arxiv(packet,out,known_aliases):
    # Compact coordinate aliases catch common system designations in survey papers.
    # These are metadata searches; association requires inspection of the returned paper.
    def aliases(ra,dec):
        hours = ra/15
        hm = int(hours)*100+int(hours*60)%60
        dm = int(abs(dec))*100+int(abs(dec)*60)%60
        sign = '+' if dec>=0 else '-'
        return {f'J{hm:04d}{sign}{dm:04d}',
                f'J{hm:04d}{int(hours*3600)%60:02d}{sign}{dm:04d}{int(abs(dec)*3600)%60:02d}'}
    terms=set()
    for c in packet['components']:
        terms.update(aliases(c['ra'],c['dec']))
        terms.add(c['source_id'])
    terms.update(aliases(*packet['imaging']['lens_model']['center_radec_deg']))
    terms.update(a.strip() for a in known_aliases if a.strip())
    params={'search_query':' OR '.join('all:"'+t+'"' for t in sorted(terms)), 'max_results':50}
    try:
        data,rec=get('https://export.arxiv.org/api/query',params,out/'arxiv_atom.xml.gz')
        doc=ET.fromstring(data)
        total=int(doc.findtext('{http://a9.com/-/spec/opensearch/1.1/}totalResults'))
        papers=[]
        for e in doc.findall('{http://www.w3.org/2005/Atom}entry'):
            title=' '.join(e.findtext('{http://www.w3.org/2005/Atom}title','').split())
            abstract=e.findtext('{http://www.w3.org/2005/Atom}summary','')
            if title.lower()=='error':
                raise ValueError(abstract[:160])
            papers.append(dict(title=title,url=e.findtext('{http://www.w3.org/2005/Atom}id'),
                               matched_alias=[t for t in terms if t.lower() in (title+' '+abstract).lower()],
                               topic=bool(TOPIC.search(title+' '+abstract))))
        return dict(status='complete' if total==len(papers) else 'truncated',total=total,queries=sorted(terms),papers=papers,receipt=rec)
    except Exception as exc:
        return dict(status='unavailable',queries=sorted(terms),error=str(exc)[:180])


def host_publications(packet):
    # Three directly germane imaging/spectroscopy samples and a SIMBAD-linked
    # optical-polarimetry sample serve as independently auditable full-text checks.
    papers = [('Hutchings et al. low-z host galaxy imaging','astro-ph/0304075'),
              ('Marble et al. HST WFPC2 red QSO host imaging','astro-ph/0303184'),
              ('Rose et al. low-z red AGN spectroscopy','1303.7112'),
              ('Smith et al. near-IR-selected QSO polarimetry','astro-ph/0112334')]
    ra,dec=packet['imaging']['lens_model']['center_radec_deg']
    h=ra/15
    hm=f'{int(h):02d}{int(h*60)%60:02d}'
    d=f'{int(abs(dec)):02d}'
    dm=f'{int(abs(dec)*60)%60:02d}'
    sign = r'[+＋]' if dec>=0 else r'[-−–]'
    # Full names and common coordinate abbreviations. The shortest names
    # are review leads; a similar coordinate prefix need not refer to this object.
    pattern=re.compile(r'(?<![0-9])'+hm+r'(?:[ :.]*[0-9]{2}(?:[.][0-9]*)?)?\s*'+sign+
                       r'\s*'+d+r'(?:[ :]*'+dm+r'[0-9.]*)?(?![0-9])',re.I)
    def inspect(item):
        title,code=item
        errors=[]
        for host in ('arxiv.org','export.arxiv.org'):
            url='https://'+host+'/pdf/'+code
            try:
                r=requests.get(url,timeout=(6,32))
                r.raise_for_status()
                if not r.content.startswith(b'%PDF'):
                    raise ValueError('non-PDF response')
                try:
                    proc=subprocess.run(['pdftotext','-layout','-','-'],input=r.content,
                                        capture_output=True,check=True,timeout=35)
                    pages=proc.stdout.decode('utf-8','replace').split(chr(12))
                except (OSError,subprocess.SubprocessError):
                    from pypdf import PdfReader
                    pages=[p.extract_text() or '' for p in PdfReader(io.BytesIO(r.content)).pages]
                if len(' '.join(pages).strip())<900:
                    raise ValueError('no searchable full text')
                hits=[]
                for page_no,page in enumerate(pages,1):
                    for line in page.splitlines():
                        if pattern.search(line):
                            hits.append(dict(page=page_no,coordinate_fragment=pattern.search(line).group(0),
                                             row=' '.join(line.split())[:230]))
                return dict(article=title,arxiv=code,url=r.url,sha256=hashlib.sha256(r.content).hexdigest(),
                            bytes=len(r.content),status='searched',pages=len([p for p in pages if p.strip()]),
                            coordinate_matches=hits,prior_errors=errors)
            except Exception as exc:
                errors.append(dict(url=url,error=str(exc)[:160]))
        return dict(article=title,arxiv=code,status='unavailable',errors=errors)
    with ThreadPoolExecutor(2) as pool:
        return list(pool.map(inspect,papers))


def mast(packet,out):
    ra,dec=packet['imaging']['lens_model']['center_radec_deg']
    request=dict(service='Mast.Caom.Cone',params=dict(ra=ra,dec=dec,radius=2/3600),format='json',
                 pagesize=2000,page=1)
    try:
        raw,rec=get('https://mast.stsci.edu/api/v0/invoke',{'request':json.dumps(request)},out/'mast.json.gz')
        doc=json.loads(raw)
        if doc.get('status')!='COMPLETE' or not isinstance(doc.get('data'),list):
            raise ValueError('MAST status '+str(doc.get('status')))
        selected=[{k:r.get(k) for k in ('obs_collection','obs_id','obsid','instrument_name',
                  'filters','target_name','proposal_id','s_region','dataproduct_type','dataURL')}
                  for r in doc['data'] if r.get('obs_collection') in ('HST','HLA','JWST')]
        partial=len(doc['data'])>=2000 or int(doc.get('paging',{}).get('pagesFiltered',1))>1
        return dict(status='truncated' if partial else 'complete',receipt=rec,total_returned=len(doc['data']),
                    high_resolution_observations=selected)
    except Exception as exc:
        return dict(status='unavailable',error=str(exc)[:220])


def run(evidence,out):
    out.mkdir(parents=True,exist_ok=True)
    packets = [json.loads(x) for x in (evidence/'review_packets.jsonl').read_text().splitlines() if x.strip()]
    selected=[p for p in packets if not p['prior'] and p['imaging']['imaging_model_pass']]
    with ThreadPoolExecutor(3) as pool:
        tables=list(pool.map(lambda t: match_table(t,selected,out),TABLES)) if selected else []
    objects=[]
    sim_path=out.parent/'literature'/'individual_literature.json'
    if not sim_path.exists():
        sim_path=evidence/'secondary'/'literature'/'individual_literature.json'
    sim=json.loads(sim_path.read_text()) if sim_path.exists() else []
    for p in selected:
        folder=out/p['system_id']
        folder.mkdir(exist_ok=True)
        aliases=[a['id'] for system in sim if system['system_id']==p['system_id'] for a in system['aliases']]
        objects.append(dict(system_id=p['system_id'],components=[c['source_id'] for c in p['components']],
                            ned=ned(p,folder), arxiv=arxiv(p,folder,aliases),mast=mast(p,folder),
                            previous_red_quasar_samples=host_publications(p) if any('2MASS' in a for a in aliases) else []))
    (out/'object_history.json').write_text(json.dumps(dict(tables=tables,objects=objects,
        scope='Eight arcseconds from any image for the published pair tables; eight arcseconds from the mass centre for NED. NED bibliographic tags and arXiv title/abstract matching have incomplete coverage. A catalogue row establishes previous selection; per-object classifications and papers require human inspection.'),indent=2)+'\n')


if __name__=='__main__':
    run(Path(sys.argv[1]),Path(sys.argv[2]))
