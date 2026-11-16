"""Summarize reported MPC magnitudes; do not infer calibrated comet photometry."""
import argparse,json,hashlib,collections,statistics
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--input',required=True);p.add_argument('--output',required=True);a=p.parse_args();src=Path(a.input);rows=[]
for line in src.read_text().splitlines():
 if len(line)!=80:continue
 date=line[15:32].strip()
 if not ('2025 12 20'<=date<='2026 01 05.999999'):continue
 try:mag=float(line[65:70])
 except ValueError:continue
 rows.append(dict(date=date,mag=mag,band=line[70],site=line[77:80],line=line))
groups=collections.defaultdict(list)
for r in rows:groups[(r['date'][:10],r['site'],r['band'])].append(r['mag'])
summary=[dict(date=k[0],site=k[1],band=k[2],n=len(v),median=statistics.median(v),minimum=min(v),maximum=max(v)) for k,v in sorted(groups.items())]
Path(a.output).write_text(json.dumps(dict(url='https://neoexchange.lco.global/target/156360/measurements/mpc/download/',sha256=hashlib.sha256(src.read_bytes()).hexdigest(),all_lines=len(src.read_text().splitlines()),selected_magnitude_records=len(rows),caveat='Heterogeneous reported MPC magnitudes; not calibrated dust photometry or IR confirmation; multi-observatory aggregation.',groups=summary,records=rows),indent=2)+'\n')
