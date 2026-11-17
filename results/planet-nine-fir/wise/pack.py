"""Compact the per-field IRSA ledgers and W2 completeness records for Git."""
import argparse,csv,gzip,json
from pathlib import Path
HERE=Path(__file__).resolve().parent

def pack(source,dest):
    dest.mkdir(parents=True,exist_ok=True)
    receipts=[]; results=[]
    for p in sorted((source/'receipts').glob('*.json')):receipts.append(json.loads(p.read_text()))
    for p in sorted((source/'search_results').glob('*.json')):results.append(json.loads(p.read_text()))
    for name,records in [('acquisition',receipts),('field_results',results)]:
        compressed=gzip.compress(('\n'.join(json.dumps(x,separators=(',',':')) for x in records)+'\n').encode(),mtime=0)
        (dest/(name+'.jsonl.gz')).write_bytes(compressed)
    result={r['target']:r for r in results}
    fields=['target','science_frames_found','selected_frames','usable_frames','searched_tracks','eligible_tracks','crossings','eligible_controls','control_crossings','injection_draws','rec_W2_13','rec_W2_14','rec_W2_14.5','rec_W2_15','rec_W2_16','max_track_score','status']
    with (dest/'field_summary.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(fields)
        for a in receipts:
            name=a['target'];r=result.get(name)
            start=[name,a['science_count'],len(a['products']),sum(all(k+'_sha256' in x for k in ['int','msk']) for x in a['products'])]
            if r:
                n={row[0]:row[1] for row in r['injections']};highest=r['highest_mean_scores']
                w.writerow(start+[r[k] for k in ['templates','eligible','threshold_crossings','eligible_controls','control_crossings']]+[len(r['injections']) and r['injections'][0][2],*[n[m] for m in [13,14,14.5,15,16]],highest[0]['mean'] if highest else '', 'searched' if r['eligible'] else 'no eligible test'])
            else:w.writerow(start+['']*12+['insufficient valid images or controls'])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--dest',type=Path,default=HERE/'compact');args=p.parse_args();pack(HERE,args.dest)
