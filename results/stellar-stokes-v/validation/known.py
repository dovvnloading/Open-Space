"""Identify published radio counterparts and count independent, dated V observations."""
import argparse
import hashlib
import json
import re
from pathlib import Path
from urllib.request import urlopen

import pandas as pd
from lxml import html

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
def sha(b):
    return hashlib.sha256(b).hexdigest()

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('cache',type=Path)
    args=ap.parse_args()
    def read(name,url):
        p=args.cache/name
        p.parent.mkdir(parents=True,exist_ok=True)
        if not p.exists():
            with urlopen(url,timeout=180) as f:
                p.write_bytes(f.read())
        data=p.read_bytes()
        receipts[name]=dict(url=url,sha256=sha(data),bytes=len(data))
        return data
    receipts={}
    srsc='https://radiostars.org/data/download/'
    stars=read('srsc/stars.dat',srsc+'stars.dat/').decode().splitlines()
    radio=read('srsc/radio.dat',srsc+'radio.dat/').decode().splitlines()
    read('srsc/ReadMe',srsc+'ReadMe/')
    hosts=pd.read_csv(ROOT/'samples/hosts.csv.gz',dtype=str,keep_default_na=False)
    host_map=hosts[hosts.gaia_dr3_id.ne('')].groupby('gaia_dr3_id').hostname.agg(list)
    ids={}
    for row in stars:
        gid=row[40:68].strip().removeprefix('Gaia DR3 ')
        for name in host_map.get(gid,[]):
            ids.setdefault(row[:10],[]).append((name,gid,row[11:39].strip(),row[306:309].strip()))
    fields={'radio_id':(17,53),'field':(55,75),'start_utc':(77,95),'exposure_s':(97,104),
            'radio_ra_deg':(106,115),'radio_dec_deg':(117,125),'survey':(215,230),'mhz':(139,144),
            'i_peak_mjy_beam':(146,153),'v_peak_catalog_mjy_beam':(232,237),'v_peak_err_mjy_beam':(239,243),
            'v_int_catalog_mjy':(245,250),'v_int_err_mjy':(252,256),'telescope':(258,264),'ref':(278,304)}
    rows=[]
    for r in radio:
        for name,gid,simbad,multiple in ids.get(r[:10],[]):
            rec={k:r[a-1:b].strip() for k,(a,b) in fields.items()}
            rec={k:('' if v=='None' else v) for k,v in rec.items()}
            rows.append(dict(hostname=name,gaia_dr3_id=gid,srsc=r[:10],published_alias=simbad,
                radio_multiple=multiple,has_v=bool(rec['v_peak_catalog_mjy_beam'] or rec['v_int_catalog_mjy']),**rec))
    m=pd.DataFrame(rows).sort_values(['hostname','start_utc','radio_id'])
    m.to_csv(HERE/'srsc_hosts.csv',index=False)
    # Pull published names from all three tables. Aliases here are limited to the
    # V-LoTSS names encountered locally, with a separate column for the source form.
    yiu=html.fromstring(read('yiu.html','https://arxiv.org/html/2312.07162'))
    pub=[]
    for fig in yiu.xpath('//figure[.//table]'):
        caps=fig.xpath('.//figcaption')
        desc=' '.join(caps[0].text_content().split()) if caps else ''
        tab=re.match(r'Table (A\.[123])',desc)
        if tab:
            for row in fig.xpath('.//table//tr[td]'):
                name=' '.join(row.xpath('./td')[0].text_content().split())
                pub.append(dict(table=tab[1],name=name))
    if not any(x['table']=='A.2' for x in pub):
        raise ValueError('Expected Yiu Stokes V table A.2')
    def key(s):
        s=s.casefold().replace('σ2\\sigma^{2} crb (tz crb)','sig crb').replace('σ\\sigma gem','sig gem')
        return re.sub(r'\W','',s.split(' (')[0])
    published=pd.DataFrame(pub)
    yiu_names={key(a.name):a.name for a in published[published.table=='A.2'].itertuples()}
    v=pd.read_csv(ROOT/'lofar/data/vlotss_matches.csv',dtype=str).query("target_position=='on'")
    v=v[['vlotss_name','published_name','date','field']].drop_duplicates().sort_values('vlotss_name')
    v['yiu_table_A2_name']=v.published_name.map(lambda s:yiu_names.get(key(s),''))
    v['status']='published V-LoTSS counterpart'
    v.to_csv(HERE/'lofar_known.csv',index=False)
    racs=pd.read_csv(ROOT/'askap/racs_matches.csv',dtype=str)
    old=racs[['radio_id','source_association','source_type','sbid','mjd']].drop_duplicates()
    old.sort_values(['source_association','mjd']).to_csv(HERE/'racs_known.csv',index=False)
    data=dict(input_files=receipts,srsc_stars=len(stars),srsc_measurements=len(radio),
              host_matches=m.hostname.nunique(),v_hosts=m[m.has_v].hostname.nunique(),
              all_srsc_host_rows=len(m),srsc_v_rows=int(m.has_v.sum()),
              vlotss_unique=len(v),yiu_A2_name_overlaps=int(v.yiu_table_A2_name.ne('').sum()),
              racs_unique_ids=old.radio_id.nunique(),racs_source_types=old.source_type.value_counts().to_dict())
    (HERE/'known_receipt.json').write_text(json.dumps(data,indent=2)+'\n')
if __name__=='__main__':
    main()
