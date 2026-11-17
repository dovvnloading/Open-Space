"""Four-header eligibility audit, no access to primary science arrays."""
import argparse
import json
from pathlib import Path
from geometry import read_header
from metadata_audit import audit
from footprints import footprints

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('files',nargs=4,type=Path);args=parser.parse_args()
    h=[read_header(p) for p in args.files];a=[audit(p) for p in args.files]
    for i in range(4):
        for j in range(i):
            if a[j]['date_end']>=a[i]['date_beg'] or set(a[j]['timestamps'])&set(a[i]['timestamps']) or set(a[j]['operand_ids'])&set(a[i]['operand_ids']):raise ValueError('acquisition overlap')
    patches=footprints(h[1],h)
    print(json.dumps({'total':len(patches),'eligible':sum(x['eligible'] for x in patches),'nodes_per_patch':len(patches[0]['maps'][0]['x']),'patches':[{'patch_id':x['patch_id'],'eligible':x['eligible'],'reasons':x['reasons']} for x in patches]},indent=2))
