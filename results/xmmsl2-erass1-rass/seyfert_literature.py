"""Targeted source-name and object-bibliography audit for the RASS-bright Seyfert."""
import argparse
import concurrent.futures as futures
import hashlib
import io
import json
import re
import shutil
import subprocess
import threading
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

import requests
from astropy.io.votable import parse
from bs4 import BeautifulSoup

SIMBAD = 'https://simbad.cfa.harvard.edu/simbad/sim-tap/sync'
OA = 'https://api.openalex.org/works'
CR = 'https://api.crossref.org/works'
PRIMARY = '10.1046/j.1365-8711.2001.04763.x'
CONTROL = '10.1051/0004-6361/201628737'
SOURCES = {
    'seyfert': {'pos': (194.71454, 23.92411), 'id': '2MASX J12585149+2355268',
                'aliases': ['KUG 1256+241', 'RX J1258.8+2355', '2MASX J12585149+2355268',
                            '1RXS J125851.4+235532', '2RXS J125851.4+235533', 'XMMSL2 J125851.2+235522',
                            '2SXPS J125851.4+235529', '3eRASS J125851.6+235526', '3eRASS J125851.6+235527']},
    'control': {'pos': (97.69125, -60.51944), 'id': 'XMMSL1 J063045.9-603110',
                'aliases': ['XMMSL1 J063045.9-603110', 'XMMSL2 J063045.7-603110']}}
TERMS = {'change': r'\b(?:variab\w*|chang\w*|declin\w*|turn[ -]?off|switch\w*|fad\w*|dimm\w*|brighten\w*|transient\w*|flare\w*|outburst\w*|tidal disruption)\b',
         'rosat': r'\b(?:ROSAT|2RXS|1RXS|RASS)\b', 'erosita': r'\b(?:eROSITA|eRASS\d?|[123]eRASS)\b',
         'slew': r'\b(?:XMMSL\d?|slew)\b'}
receipts = []
lock = threading.Lock()
services = {OA: threading.BoundedSemaphore(2), CR: threading.BoundedSemaphore(2)}


def clean(s):
    return ' '.join(str(s or '').translate(str.maketrans({'−':'-', '–':'-', '＋':'+'})).split())


def norm(s):
    return re.sub(r'[^a-z0-9+-]', '', clean(s).lower())


def fetch(url, params=None, cap=8_000_000):
    item = {'url':url, 'parameters':params, 'utc':datetime.now(timezone.utc).isoformat()}
    try:
        with requests.get(url, params=params, timeout=(5,12), stream=True,
                          headers={'User-Agent':'Open-Space astronomical archive study/1.0 (https://github.com/dovvnloading/Open-Space)'}) as r:
            item.update(status=r.status_code, final=r.url)
            r.raise_for_status()
            data=bytearray()
            for chunk in r.iter_content(65536):
                if len(data)+len(chunk)>cap:
                    raise ValueError('response exceeds size cap')
                data.extend(chunk)
        item.update(sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
        return bytes(data)
    except Exception as exc:
        item['error']=str(exc)[:430]
        raise
    finally:
        with lock:
            receipts.append(item)


def tap(sql):
    raw=fetch(SIMBAD, {'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'votable','MAXREC':5000,'QUERY':sql})
    v=parse(io.BytesIO(raw)); flags=list(v.infos)+[i for r in v.resources for i in r.infos]
    for i in flags:
        if i.name=='QUERY_STATUS' and str(i.value).upper() in ('ERROR','OVERFLOW'):
            raise ValueError(str(i.content or i.value))
    return v.get_first_table().to_table(use_names_over_ids=True).to_pandas().fillna('').to_dict('records')


def mentions(text, aliases):
    text=clean(text)
    matches=[]
    for alias in aliases:
        # Retain coordinate signs and digit boundaries, allow catalogue spacing/punctuation differences.
        chars=norm(alias)
        if len(chars)<7: continue
        rx=r'(?<![A-Za-z0-9])'+r'[\s.:_]*'.join(re.escape(x) for x in chars)+r'(?![0-9])'
        for m in list(re.finditer(rx, text, re.I))[:3]:
            context=text[max(0,m.start()-160):m.end()+190]
            near=[k for k,pat in TERMS.items() if re.search(pat,context,re.I)]
            matches.append({'alias':alias,'local_terms':near,'context':context})
    return matches


def record(title,doi,year,abstract,aliases,**more):
    title,abstract=clean(title),clean(abstract)
    body=title+' '+abstract
    return {'title':title,'doi':doi or '', 'year':year, **more, 'abstract_available':bool(abstract),
            'alias_contexts':mentions(body,aliases), 'terms_in_title_or_abstract':[k for k,p in TERMS.items() if re.search(p,body,re.I)]}


def simbad_source(key):
    source=SOURCES[key]; ra,dec=source['pos']; name=source['id']; aliases=source['aliases']
    literal="'"+name.replace("'","''")+"'"
    sql=('SELECT DISTINCT b.oid,b.main_id,b.ra,b.dec,b.otype,i2.ids FROM ident AS i '
         'JOIN basic AS b ON i.oidref=b.oid JOIN ids AS i2 ON i2.oidref=b.oid WHERE i.id='+literal)
    fields={}; errors={}; linked=[]
    for tag,query in [('exact_identifier',sql),('positional_check',"SELECT b.oid,b.main_id,b.ra,b.dec,b.otype,i2.ids FROM basic AS b JOIN ids AS i2 ON i2.oidref=b.oid WHERE 1=CONTAINS(POINT('ICRS',b.ra,b.dec),CIRCLE('ICRS',%.7f,%.7f,0.0008333333))"%(ra,dec))]:
        try: fields[tag]=tap(query)
        except Exception as e: errors[tag]=str(e)[:350]
    exact=fields.get('exact_identifier',[])
    if len(exact)!=1:
        exact=[x for x in fields.get('positional_check',[]) if norm(x['main_id'])==norm(name)]
    if len(exact)==1:
        obj=exact[0]; linked=sorted({clean(x) for x in obj['ids'].split('|') if clean(x)} | {name})
        try:
            rows=tap('SELECT r.bibcode,r.title,r.doi,r.abstract FROM has_ref AS h JOIN ref AS r ON h.oidbibref=r.oidbib WHERE h.oidref='+str(int(obj['oid'])))
            fields['all_bibliography']=[record(r['title'],r.get('doi'),str(r['bibcode'])[:4],r.get('abstract'),set(aliases+linked),bibcode=r['bibcode']) for r in rows]
        except Exception as e: errors['all_bibliography']=str(e)[:350]
    else: errors['all_bibliography']='could not verify unique source identity'
    return {**fields,'errors':errors,'resolved_aliases':linked}


def abstract_oa(work):
    positions=work.get('abstract_inverted_index') or {}
    return ' '.join(t for i,t in sorted((i,t) for t,p in positions.items() for i in p))


def search(provider, key, query, doi='', author_refs=False):
    aliases=list(SOURCES[key]['aliases'])
    if query and query not in aliases and 'ROSAT eROSITA variability' not in query:
        aliases.append(query)
    base=OA if provider=='OpenAlex' else CR
    url=base; params={}; cap=10; works=[]; details={}
    if doi:
        url += '/'+('https://doi.org/' if base==OA else '')+quote(doi,safe='/')
    else:
        params=({'search':query,'per_page':cap} if base==OA else {'query.bibliographic':query,'rows':cap})
    with services[base]:
        payload=json.loads(fetch(url,params))
    if base==OA:
        items=[payload] if doi else payload['results']; total=1 if doi else payload['meta'].get('count')
        for w in items:
            work=record(w.get('display_name'),w.get('doi'),w.get('publication_year'),abstract_oa(w),aliases,
                        work_id=w.get('id'),url=(w.get('primary_location') or {}).get('landing_page_url'))
            if author_refs:
                work['referenced_work_ids']=w.get('referenced_works',[])
                work['cited_by_count']=w.get('cited_by_count')
                details['seed_id']=str(w.get('id','')).rsplit('/',1)[-1]
            works.append(work)
    else:
        msg=payload['message']; items=[msg] if doi else msg['items']; total=1 if doi else msg.get('total-results')
        for w in items:
            a=BeautifulSoup(str(w.get('abstract','')),'html.parser').get_text(' ',strip=True)
            date=(w.get('issued',{}).get('date-parts') or [[None]])[0]
            work=record(' '.join(w.get('title') or []),w.get('DOI'),date[0] if date else None,a,aliases,url=w.get('URL'))
            if author_refs:
                work['reference_dois']=[x['DOI'] for x in w.get('reference',[]) if x.get('DOI')]
            works.append(work)
    return {'provider':provider,'subject':key,'query':query,'exact_doi':doi,'total':total,'returned':len(works),'works':works,**details}


def citing(seed, key):
    # Search the descendants of the source-selected 2001 parent article. Citing it alone never proves object identity.
    p={'filter':'cites:'+seed,'search':'ROSAT eROSITA variability','per_page':10}
    with services[OA]: data=json.loads(fetch(OA,p))
    return {'query':p,'total':data['meta'].get('count'),'works':[record(w.get('display_name'),w.get('doi'),w.get('publication_year'),abstract_oa(w),SOURCES[key]['aliases'],work_id=w.get('id')) for w in data['results']]}


def primary_text(key, arxiv_id, pdf=False):
    path=('pdf/' if pdf else 'abs/')+arxiv_id
    raw=fetch('https://arxiv.org/'+path)
    if pdf:
        if not raw.startswith(b'%PDF') or not shutil.which('pdftotext'):
            raise ValueError('original PDF or pdftotext unavailable')
        p=subprocess.run(['pdftotext','-layout','-','-'],input=raw,capture_output=True,timeout=12,check=True)
        body=p.stdout.decode('utf8','replace')
    else:
        doc=BeautifulSoup(raw,'html.parser'); abstract=doc.select_one('blockquote.abstract'); title=doc.select_one('h1.title')
        if not abstract or not title: raise ValueError('structured arXiv abstract/title not found')
        body=title.get_text(' ',strip=True)+' '+abstract.get_text(' ',strip=True)
    return {'subject':key,'arxiv':arxiv_id,'section':'full PDF' if pdf else 'title and abstract',
            'characters':len(body),'alias_contexts':mentions(body,SOURCES[key]['aliases']),
            'terms_in_text':[k for k,p in TERMS.items() if re.search(p,body,re.I)]}


def main(out):
    result={'query_started_utc':datetime.now(timezone.utc).isoformat(),'target_identifiers':SOURCES,
            'coverage':'Searches cover the returned OpenAlex and Crossref records, all Harvard SIMBAD references attached to each uniquely resolved object, and the specified arXiv articles. Search strings have at most ten ranked works each. Per-record exact-name contexts use titles and available abstracts unless tagged full PDF. Coordinate rounding, unindexed work, abstract omissions and missing archive metadata limit recall. A citation or a mission name alone does not report variability. No absent or failed query establishes novelty.'}
    results={}
    def collect(jobs):
        for f in futures.as_completed(jobs):
            tag=jobs[f]
            try: results[tag]=f.result()
            except Exception as e: results[tag]={'error':str(e)[:440]}
    with futures.ThreadPoolExecutor(max_workers=10) as pool:
        jobs={pool.submit(simbad_source,k):'simbad_'+k for k in SOURCES}
        for provider in ('OpenAlex','Crossref'):
            for key,terms in [('seyfert',['KUG 1256+241','RX J1258.8+2355','2MASX J12585149+2355268','XMMSL2 J125851.2+235522','2RXS J125851.4+235533','KUG 1256+241 ROSAT eROSITA variability']),('control',['XMMSL1J063045.9-603110'])]:
                for i,q in enumerate(terms): jobs[pool.submit(search,provider,key,q)]=f'{provider}_{key}_search_{i}'
            for key,doi in [('seyfert',PRIMARY),('control',CONTROL)]:
                jobs[pool.submit(search,provider,key,'',doi,True)]=f'{provider}_{key}_primary'
        for key,aid,pdf in [('control','1605.06133',False),('seyfert','astro-ph/0106311',False),('seyfert','astro-ph/0106311',True)]:
            jobs[pool.submit(primary_text,key,aid,pdf)]=f'arxiv_{key}_{"pdf" if pdf else "abstract"}'
        collect(jobs)
        more={}
        seed=results['OpenAlex_seyfert_primary'].get('seed_id')
        if re.fullmatch(r'W\d+',seed or ''): more[pool.submit(citing,seed,'seyfert')]='OpenAlex_targeted_citing_works'
        # Prefer other historical identifying names, not positional neighbours or additional searches of modern source integers.
        linked=results['simbad_seyfert'].get('resolved_aliases',[])
        for a in linked:
            if a not in SOURCES['seyfert']['aliases'] and re.match(r'^([MP]G[C ]|Mrk |[CH]S |[HB89]+ |IRAS |SDSS J|RBS )',a,re.I):
                for provider in ('OpenAlex','Crossref'):
                    more[pool.submit(search,provider,'seyfert',a)]=f'{provider}_other_historical_alias_{a}'
                if len(more)>=7: break
        collect(more)
    checks={}
    for provider in ('OpenAlex','Crossref'):
        named=results[provider+'_control_search_0']; primary=results[provider+'_control_primary']
        checks[provider]={'name_search_retrieves_known_control_doi':any(CONTROL in str(x['doi']).lower() for x in named.get('works',[])) if 'works' in named else None,
                          'doi_lookup_recovers_control':any(CONTROL in str(x['doi']).lower() for x in primary.get('works',[])) if 'works' in primary else None}
    for key,expected in [('seyfert',PRIMARY),('control',CONTROL)]:
        s=results['simbad_'+key]; refs=s.get('all_bibliography')
        checks['SIMBAD_'+key]={'references_returned':len(refs) if refs is not None else None,
                              'known_primary_doi_recovered':any(expected in str(x['doi']).lower() for x in refs) if refs is not None else None}
    result.update(checks=checks,results=results,receipts=sorted(receipts,key=lambda x:(x['url'],str(x['parameters']))),query_completed_utc=datetime.now(timezone.utc).isoformat())
    (out/'seyfert_literature.json').write_text(json.dumps(result,indent=2,default=str,allow_nan=False)+'\n')
    print('Seyfert publication audit:',json.dumps(checks),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('out',type=Path); main(p.parse_args().out)
