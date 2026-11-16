"""Download station metadata and the public numeric plot series from DGFI-TUM."""
import csv
import datetime as dt
import hashlib
import io
import json
import math
import re
import urllib.request
from pathlib import Path

OUT=Path(__file__).resolve().parent
URL="https://dahiti.dgfi.tum.de/en/virtual-stations/africa/mali/"
SITES={"mopti":(-4.2014431,14.53268014), "dire":(-3.3563,16.2666)}


def fetch(url):
    with urllib.request.urlopen(url,timeout=50) as resp:
        body=resp.read()
        return body,{"url":url,"http_date":resp.headers.get("date"),
                     "sha256_body":hashlib.sha256(body).hexdigest()}


def run():
    raw,meta=fetch(URL)
    country=json.loads(re.search(r'var targets = (\[.*?\]);',raw.decode(),re.S)[1])
    meta["chart_object"]="var targets"
    (OUT/"station_index.json").write_text(json.dumps(country,indent=2)+"\n")
    out={"index":meta,"targets":[]}
    for region,(lon,lat) in SITES.items():
        for record in country:
            if record["target_name"]!="Niger, River" or record["permission_water_level_altimetry"]!="public":
                continue
            x,y=record["longitude"],record["latitude"]
            dist=111*math.hypot((x-lon)*math.cos(math.radians(lat)),y-lat)
            if dist>60:
                continue
            station=record["id"]
            body,m=fetch(f"https://dahiti.dgfi.tum.de/en/{station}/water-level-altimetry/")
            html=body.decode()
            mission=re.findall(r'<td style="width: 250px;">([^<]+)</td>\s*<td>([^<]+)</td>',html)
            series=[json.loads(s) for s in re.findall(r'var data_water_level(?:_error)? = JSON.parse\("(\[\[.*?\]\])"\)',html)]
            height,err=series
            assert [v[0] for v in height]==[v[0] for v in err]
            r=[]
            for (stamp,z),(_,e) in zip(height,err):
                t=dt.datetime.fromtimestamp(stamp/1000,dt.UTC)
                r.append([t.isoformat(),z,e])
            name=f"station_{station}.csv"
            with (OUT/name).open("w",newline="") as f:
                writer=csv.writer(f,lineterminator="\n")
                writer.writerow(["time_utc","water_level_m","uncertainty_m"])
                writer.writerows(r)
            m.update(region=region,id=station,longitude=x,latitude=y,distance_km=dist,
                     name=record["target_name"],file=name,sha256_csv=hashlib.sha256((OUT/name).read_bytes()).hexdigest(),
                     count=len(r),first=r[0][0],last=r[-1][0],satellite_track=mission)
            print(region,station,round(dist,1),r[0][0][:10],r[-1][0][:10],len(r),flush=True)
            out["targets"].append(m)
    (OUT/"manifest.json").write_text(json.dumps(out,indent=2)+"\n")


if __name__=="__main__":
    run()
