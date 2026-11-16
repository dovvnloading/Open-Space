"""Describe historical-only station pairs separately and the fixed-date downstream series."""
from datetime import date
import itertools
import json
from pathlib import Path
from compare import nearest,read,save

OUT=Path(__file__).resolve().parent


def endpoints(series,y):
    return [nearest(series,date(y,7,22),0.30),nearest(series,date(y+1,1,26),0.30)]


def main():
    rows=json.loads((OUT/"manifest.json").read_text())["targets"]
    hist={name:[t for t in rows if t["region"]==name and t["last"]<"2018"] for name in ("mopti","dire")}
    pairs=[]
    for m,d in itertools.product(hist['mopti'],hist['dire']):
        for y in range(2002,2016):
            e=endpoints(read(m),y)+endpoints(read(d),y)
            mf=e[1][1]-e[0][1] if all(e[:2]) else ""
            df=e[3][1]-e[2][1] if all(e[2:]) else ""
            pairs.append(dict(mopti_id=m["id"],dire_id=d["id"],season=y,m_july=str(e[0][0]) if e[0] else "",
                              m_jan=str(e[1][0]) if e[1] else "",d_july=str(e[2][0]) if e[2] else "",
                              d_jan=str(e[3][0]) if e[3] else "",mopti_delta_m=mf,dire_delta_m=df,
                              four_dates_present=all(e),
                              mopti_within_20cm=all(e) and abs(mf)<=0.20))
    save("historical_only_pairs.csv",pairs)
    solo=[]
    for target in rows:
        starts=int(target['first'][:4]); ends=int(target['last'][:4])
        for y in range(starts,ends):
            r=endpoints(read(target),y)
            solo.append(dict(site=target['region'],id=target['id'],season=y,
                             july_date=str(r[0][0]) if r[0] else "",july_wse_m=r[0][1] if r[0] else "",
                             january_date=str(r[1][0]) if r[1] else "",january_wse_m=r[1][1] if r[1] else "",
                             change_m=r[1][1]-r[0][1] if all(r) else ""))
    save("single_station_seasons.csv",solo)


if __name__=='__main__':
    main()
