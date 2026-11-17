"""Explicitly development-only VAST reader; do not execute before review.
PyArrow may decode entire row-group column chunks before applying predicates.
The returned table is rechecked against the allowed source and image IDs.
"""
from pathlib import Path
import csv,hashlib
import pyarrow.parquet as pq
from packet_guard import guard
from conditional_pipeline import phase_rows

COLUMNS=['source','image_id','flux_peak','flux_peak_err','local_rms','forced','flag_c4','has_siblings']

def digest(path):
 h=hashlib.sha256()
 with path.open('rb') as f:
  for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
 return h.hexdigest()

def read_development(source_id,parquet_path,image_csv):
 packet=guard();source_id=int(source_id);p=Path(parquet_path)
 if digest(Path(image_csv))!='70ce692b3c8d067d53e248b6265420a26880147bd7bc15327930880bd8395195':raise ValueError('image metadata digest mismatch')
 if digest(p)!=packet['measurement_parquet_sha256']:raise ValueError('parquet digest mismatch')
 item=packet['sources'][str(source_id)];images={i for group in item['development'] for i in group}
 hold={i for group in item['holdout'] for i in group}
 if images&hold:raise ValueError('overlapping phase')
 controls={int(r['control_source']) for r in csv.DictReader((Path(__file__).parent/'control-metadata.csv').open()) if int(r['target_source'])==source_id};sources=controls|{source_id}
 tab=pq.read_table(p,columns=COLUMNS,filters=[('source','in',list(sources)),('image_id','in',list(images))])
 # Check the predicate before converting or returning any values.
 if not set(tab.column('source').to_pylist())<=sources or not set(tab.column('image_id').to_pylist())<=images:raise ValueError('reader returned forbidden row')
 meta={}
 with Path(image_csv).open(newline='',encoding='utf-8') as f:
  for r in csv.DictReader(f):
   if int(r['image_id']) in images:meta[int(r['image_id'])]=(r['beam_bmaj'],r['beam_bmin'])
 if set(meta)!=images:raise ValueError('incomplete image metadata')
 rows=tab.to_pylist()
 for r in rows:r['beam_bmaj'],r['beam_bmin']=meta[int(r['image_id'])]
 return phase_rows(rows,images,sources)
