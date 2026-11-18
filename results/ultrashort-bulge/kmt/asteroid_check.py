"""Topocentric known-object checks across the CTIO peak of KMT-2024-BLG-3326."""
import json
from urllib.parse import urlencode
from pathlib import Path
from freeze import HERE,get
from neighbours import angle


def main():
    import argparse
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--cache',type=Path,required=True)
    args=p.parse_args()
    r=next(r for r in map(json.loads,(args.cache/'index.jsonl').read_text().splitlines()) if r['event']=='KMT-2024-BLG-3326')
    recs=[];targets=(2460570.48,2460570.55,2460570.62)
    for hjd in targets:
        # HJD-UTC displacement is < 8.4 minutes. Keep explicitly approximate dates and a 5' radius.
        url='https://ssp.imcce.fr/webservices/skybot/api/conesearch.php?'+urlencode(dict(RA=f"{angle(r['ra'],True):.9f}",
            DEC=f"{angle(r['dec']):.9f}",SR='0.083333',TIME=hjd,**{'-loc':'807','-mime':'text','-output':'basic'}))
        try:
            raw,rec=get(url);rec.update(approx_epoch_jd=hjd,reply=raw.decode())
        except OSError as exc:
            rec=dict(url=url,approx_epoch_jd=hjd,status=getattr(exc,'code',None),error=str(exc))
        recs.append(rec)
    (HERE/'asteroid-3326.json').write_text(json.dumps(recs,indent=2)+'\n')


if __name__=='__main__':main()
