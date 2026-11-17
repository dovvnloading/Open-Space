"""Replay the bounded IRSA image eligibility check from saved responses."""
import csv
import io
import json
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parent
ns = {'v': 'http://www.ivoa.net/xml/VOTable/v1.3'}
root = ET.parse(ROOT / 'inputs/sia-vis.xml').getroot()
status = root.find('.//v:INFO[@name="QUERY_STATUS"]', ns).attrib['value']
if status != 'OK':
    raise SystemExit('SIA response is incomplete: ' + status)
table = root.find('.//v:TABLE', ns)
names = [field.attrib['name'] for field in table.findall('v:FIELD', ns)]
rows = [dict(zip(names, [(cell.text or '') for cell in tr.findall('v:TD', ns)])) for tr in table.findall('.//v:TR', ns)]
science = sorted((r for r in rows if r['dataproduct_subtype'] == 'science' and r['energy_bandpassname'] == 'VIS'), key=lambda r: r['obs_id'])
with (ROOT / 'inventory.csv').open('w', newline='', encoding='utf-8') as f:
    w = csv.writer(f, lineterminator='\n')
    w.writerow(['obs_id', 'tile', 't_min', 't_max', 'access_url'])
    for r in science:
        w.writerow([r['obs_id'], r['obs_id'][:9], r['t_min'], r['t_max'], r['access_url']])
header_results = []
request = json.loads((ROOT / 'inputs/requests.json').read_text())
for item in request['headers']:
    data = (ROOT / 'inputs' / item['file']).read_bytes()
    cards = [data[i:i+80].decode('ascii', errors='replace') for i in range(0, len(data), 80)]
    names_present = [k for k in ('DATE-OBS','MJD-OBS','MJD-BEG','MJD-END') if any(c[:8].strip() == k for c in cards)]
    header_results.append({'obs_id': item['observation_id'], 'date_keys': names_present, 'http_status': item['http_status']})
summary = {'sia_status': status, 'all_rows': len(rows), 'vis_science_rows': len(science), 'tiles': sorted(set(r['obs_id'][:9] for r in science)), 'science_rows_with_both_times': sum(bool(r['t_min'] and r['t_max']) for r in science), 'examined_headers': header_results, 'verified_disjoint_visits': None, 'eligible_three_visit_sets': None}
(ROOT / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n', encoding='utf-8')
print(json.dumps(summary, indent=2))
