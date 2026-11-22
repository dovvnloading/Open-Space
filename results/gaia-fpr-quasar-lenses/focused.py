"""Follow all previously unlisted two-filter morphology passes, including unanchored fields."""
import csv
import gzip
import hashlib
import io
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

import deflector
import epochs

FPR_INPUTS = {
    1:'Ducourant et al. Gaia DR3 morphology',2:'Gravitational Lens database',
    3:'Gravitational Lens candidates',4:'Milliquas 7.5b',5:'Milliquas 7.5',
    6:'Milliquas 7.4d',7:'Milliquas 7.4c',8:'Milliquas 7.1b',9:'Milliquas 7.0',
    10:'Milliquas 6.4',11:'Assef R90',12:'Shu et al. 2019',13:'Assef C75',
    14:'Klioner input list'}
TABLES = {
    'gaiafpr.lens_catalogue_name':['source_id','catalogue_id'],
    'gaiadr3.gaia_source':['source_id','in_qso_candidates','classprob_dsc_combmod_quasar',
                           'classprob_dsc_combmod_galaxy','classprob_dsc_combmod_star'],
    'gaiadr3.qso_candidates':['source_id','classlabel_dsc','classlabel_dsc_joint','classlabel_oa',
        'redshift_qsoc','redshift_qsoc_lower','redshift_qsoc_upper','flags_qsoc','zscore_qsoc',
        'ccfratio_qsoc','astrometric_selection_flag','gaia_crf_source','vari_best_class_name',
        'vari_agn_membership_score'],
    'gaiadr3.astrophysical_parameters':['source_id'] + [
        'classprob_dsc_'+method+'_'+c for method in ('combmod','specmod')
        for c in ('quasar','galaxy','star','whitedwarf','binarystar')]}


def lookup(table,ids,out):
    """Keep source IDs as exact decimal strings, and distinguish a failed query from no row."""
    if not ids:
        return [], dict(status='no_inputs',table=table)
    if not all(s.isdecimal() for s in ids):
        raise ValueError('Noninteger Gaia ID in a publication packet')
    cols = TABLES[table]
    query = f"SELECT {','.join(cols)} FROM {table} WHERE source_id IN ({','.join(sorted(ids))})"
    receipt = dict(table=table,query=query,source_ids=sorted(ids),errors=[])
    for endpoint in ('https://gaia.aip.de/tap/sync','https://gea.esac.esa.int/tap-server/tap/sync'):
        try:
            resp = requests.get(endpoint, params=dict(QUERY=query,LANG='ADQL',REQUEST='doQuery',
                                                       FORMAT='csv',MAXREC=10000),timeout=(10,55))
            resp.raise_for_status()
            reader=csv.DictReader(io.StringIO(resp.content.decode('utf-8-sig')))
            if not reader.fieldnames or set(cols)-{c.lower() for c in reader.fieldnames}:
                raise ValueError('Invalid query response: '+resp.text[:180])
            rows=[{k.lower():v.strip() if v else None for k,v in r.items()} for r in reader]
            returned={r['source_id'] for r in rows}
            if len(rows)>=10000 or not returned<=ids:
                raise ValueError('Unexpected or truncated query result')
            if table in ('gaiadr3.gaia_source','gaiafpr.lens_catalogue_name') and returned!=ids:
                raise ValueError('Missing mandatory parent/source IDs')
            path=out/(table+'.csv.gz')
            path.write_bytes(gzip.compress(resp.content,mtime=0))
            receipt.update(status='complete',endpoint=endpoint,sha256=hashlib.sha256(resp.content).hexdigest(),
                           raw=path.name,bytes=len(resp.content),rows=len(rows),no_catalogue_row=sorted(ids-returned))
            return rows,receipt
        except Exception as exc:
            receipt['errors'].append(dict(endpoint=endpoint,error=str(exc)[:450]))
    return [],dict(receipt,status='unavailable')


def run(evidence,out):
    evidence,out=Path(evidence),Path(out)
    out.mkdir(parents=True,exist_ok=True)
    packet_file=evidence/'review_packets.jsonl'
    packets=[json.loads(line) for line in packet_file.read_text().splitlines() if line.strip()]
    packets=[p for p in packets if str(p['prior']).lower()=='false' and p['imaging']['imaging_model_pass']]
    systems=[]
    parents=set()
    targets=set()
    for packet in packets:
        c=packet['components']
        parent=str(c[0]['field_id'])
        if not all(str(row['field_id'])==parent for row in c):
            raise ValueError('Inconsistent FPR field identifiers')
        parents.add(parent)
        targets.update(str(row['source_id']) for row in c)
        systems.append(dict(id=packet['system_id'],field_id=parent,components=c,
                           center_radec_deg=packet['imaging']['lens_model']['center_radec_deg']))
    # Retrieve component and parent classifications separately when the FPR seed is outside the selected subset.
    classified=targets|parents
    queries={}
    receipts={}
    with ThreadPoolExecutor(3) as pool:
        epoch_job=pool.submit(epochs.assess,systems,out/'epochs')
        deblend_job=pool.submit(deflector.assess,systems,out/'deflectors')
        for table in TABLES:
            ids=parents if table.startswith('gaiafpr.') else classified
            queries[table],receipts[table]=lookup(table,ids,out)
        variability={r['id']:r['status'] for r in epoch_job.result()}
        deblends=deblend_job.result()
    result=[]
    bibliography=[]
    for packet,sysobj in zip(packets,systems):
        ident,parent=sysobj['id'],sysobj['field_id']
        ids=[str(c['source_id']) for c in sysobj['components']]
        members=set(ids)|{parent}
        rows={t:[r for r in raw if r['source_id'] in (members if not t.startswith('gaiafpr.') else {parent})]
              for t,raw in queries.items()}
        for r in rows['gaiafpr.lens_catalogue_name']:
            r['catalogue_name']=FPR_INPUTS.get(int(r['catalogue_id']),'unmapped')
        result.append(dict(system_id=ident,field_id=parent,component_source_ids=ids,
            fpr_component_ids=[int(c['image_id']) for c in sysobj['components']],
            external_quasar_anchor=packet['catalogue_flags']['qso_evidence'],prior_disposition=packet['disposition'],
            measured_model_centre_radec_deg=sysobj['center_radec_deg'],catalogues=rows,
            epoch_status=variability.get(ident),deblend_results=[r for r in deblends if r['system_id']==ident]))
        ra,dec=sysobj['center_radec_deg']
        bibliography.append(dict(system_id=ident,ra=ra,dec=dec))
    (out/'focused_summary.json').write_text(json.dumps(dict(
        systems=result,field_count=len(result),unique_gaia_image_sets=len({tuple(sorted(r['component_source_ids'])) for r in result}),
        input_review_packets_sha256=hashlib.sha256(packet_file.read_bytes()).hexdigest(),query_receipts=receipts,
        interpretation='FPR seed provenance and Gaia DSC/QSOC redshifts or probabilities are supporting catalogue information. No spectroscopic confirmation or promotion is assigned by this module.'),indent=2,allow_nan=False)+'\n')
    with (out/'bibliography_targets.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['system_id','ra','dec'])
        writer.writeheader()
        writer.writerows(bibliography)


if __name__=='__main__':
    run(*sys.argv[1:])
