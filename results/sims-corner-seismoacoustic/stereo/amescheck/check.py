"""Cross-check Washington groups against the NASA Ames published event selection."""
import argparse
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
EVENT = '2026-08-14 07:48:36'
EPOCH_MS = 1000 * datetime.fromisoformat(EVENT).replace(tzinfo=timezone.utc).timestamp()
URL = 'https://neo-bolide.ndc.nasa.gov/service/event/6a80a3b3f102fba8e7ba7a19'


def download(output):
    raw = urlopen(URL, timeout=45).read()
    row = json.loads(raw)['data'][0]
    extract = {k:row[k] for k in ('_id','publishedDate','lastModifiedDate','humanClassification',
                                   'nearbyLightningActivity','detectedBy','beginningAltitude','csv')}
    extract['source_url'] = URL
    extract['source_sha256'] = hashlib.sha256(raw).hexdigest()
    extract['groups'] = []
    for item in row['attachments']:
        for group in item['geoData']:
            lon,lat = group['location']['coordinates']
            extract['groups'].append([int(item['platformId'][1:]), group['time'], lon, lat, group['energy']])
    extract['group_columns']=['sat','epoch_ms','longitude_deg','latitude_deg','joules']
    output.write_text(json.dumps(extract, separators=(',',':'))+'\n')


def match(groups_path, selected_path):
    selected = json.loads(selected_path.read_text())
    with groups_path.open() as f:
        ours=[r for r in csv.DictReader(f) if r['event']==EVENT]
    accepted=set()
    differences=[]
    for sat,millis,lon,lat,energy in selected['groups']:
        # Both products report the same native L2 coordinates/energy. Restrict
        # candidates before matching the independently formatted millisecond time.
        candidates=[(i,r) for i,r in enumerate(ours) if i not in accepted and int(r['sat'])==sat
                    and abs(float(r['lon'])-lon)<1e-6 and abs(float(r['lat'])-lat)<1e-6
                    and abs(float(r['joules'])/energy-1)<1e-6]
        if candidates:
            i,r=min(candidates,key=lambda item:abs(float(item[1]['seconds'])*1000-(millis-EPOCH_MS)))
            diff=float(r['seconds'])*1000-(millis-EPOCH_MS)
            if abs(diff)<2:
                accepted.add(i)
                differences.append(diff)
    bins=stereo_start_times_ms(ours)
    output=[]
    for start_ms in bins:
        row={'start_offset_ms':start_ms}
        for sat in (18,19):
            subset=[i for i,r in enumerate(ours) if int(r['sat'])==sat and
                    start_ms<=float(r['seconds'])*1000<start_ms+40]
            row[f'noaa{sat}']=len(subset)
            row[f'ames_matched{sat}']=sum(i in accepted for i in subset)
        output.append(row)
    return {'beginningAltitude':selected['beginningAltitude'], 'ames_total':len(selected['groups']),
            'noaa_selection_total':len(ours),'matched_total':len(accepted),
            'matched_time_difference_ms':[min(differences),max(differences)],'bins':output}


def stereo_start_times_ms(rows):
    for k in range(-75,75):
        t=k*40
        sats={int(r['sat']) for r in rows if t<=float(r['seconds'])*1000<t+40}
        if len(sats)==2:
            yield t


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--groups',type=Path,default=ROOT.parent/'groups.csv')
    parser.add_argument('--selection',type=Path,default=ROOT/'ames-groups.json')
    parser.add_argument('--download',action='store_true')
    parser.add_argument('--output',type=Path,default=ROOT/'matches.json')
    args=parser.parse_args()
    if args.download:
        download(args.selection)
    args.output.write_text(json.dumps(match(args.groups,args.selection),indent=2)+'\n')
