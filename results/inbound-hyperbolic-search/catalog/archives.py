"""Search pre-discovery exposure footprints for the catalogued comet C/2025 V1."""
import csv
import hashlib
import io
import json
import math
import sys
from bisect import bisect_right
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

OUT = Path(sys.argv[1])
OUT.mkdir(exist_ok=True, parents=True)
AGENT = "Open-Space-image-search/1.0 (+https://github.com/dovvnloading/Open-Space)"
receipts = []


def fetch(name, url, params):
    uri = url + "?" + urlencode(params)
    rec = {"file": name, "url": uri}
    receipts.append(rec)
    try:
        with urlopen(Request(uri, headers={"User-Agent": AGENT}), timeout=150) as resp:
            data = resp.read()
            rec.update(status=resp.status, http_date=resp.headers.get("Date"),
                       bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
            (OUT / name).write_bytes(data)
            return data.decode("utf-8")
    except Exception as exc:
        rec["error"] = str(exc)
        raise
    finally:
        (OUT / "archive-receipts.json").write_text(json.dumps(receipts, indent=2) + "\n")


def ztf(weeks=None):
    # Six-hourly topocentric astrometric ICRF coordinates at Palomar.
    params = dict(format="json", COMMAND="'DES=C/2025 V1;'", EPHEM_TYPE="'OBSERVER'",
                  CENTER="'I41'", START_TIME="'2025-07-01'", STOP_TIME="'2025-11-02'",
                  STEP_SIZE="'6 h'", QUANTITIES="'1,9,23,36'", ANG_FORMAT="'DEG'",
                  CSV_FORMAT="'YES'", CAL_FORMAT="'JD'", EXTRA_PREC="'YES'")
    raw = fetch("horizons.json", "https://ssd.jpl.nasa.gov/api/horizons.api", params)
    result = json.loads(raw)["result"]
    before, sep, rest = result.partition("$$SOE")
    if not sep:
        raise ValueError("No Horizons ephemeris: " + result[:700])
    header = next(csv.reader([next(x for x in reversed(before.splitlines()) if "R.A." in x)]))
    ri = next(i for i,x in enumerate(header) if "R.A." in x)
    di = next(i for i,x in enumerate(header) if "DEC" in x.upper())
    pts = []
    for row in csv.reader(rest.partition("$$EOE")[0].strip().splitlines()):
        pts.append((float(row[0]), float(row[ri]), float(row[di])))
    if len(pts) < 490:
        raise ValueError("Incomplete ephemeris")
    dates = [p[0] for p in pts]
    candidates, used_ids, gaps = [], set(), []
    completed = 0
    starts = list(range(0, len(pts)-1, 28))
    if weeks is not None:
        starts = starts[-weeks:]
    for start in reversed(starts):
        group = pts[start:start+29]
        mid = group[len(group)//2]
        # Deliberately generous axis-aligned envelope for each week's track.
        dx = lambda a,b: (a-b+180)%360-180
        width = 2*max(abs(dx(p[1],mid[1])) for p in group) + 1
        height = 2*max(abs(p[2]-mid[2]) for p in group) + 1
        params = {"POS": f"{mid[1]:.7f},{mid[2]:.7f}", "SIZE": f"{width:.6f},{height:.6f}",
                  "INTERSECT": "OVERLAPS", "WHERE": f"obsjd>={group[0][0]} AND obsjd<{group[-1][0]}",
                  "COLUMNS": "filefracday,field,ccdid,qid,filtercode,obsjd,ra,dec,ra1,dec1,ra2,dec2,ra3,dec3,ra4,dec4",
                  "ct":"csv"}
        try:
            raw = fetch(f"ztf-{start:03d}.csv","https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci",params)
        except Exception as exc:
            gaps.append([start,str(exc)])
            continue
        table = csv.DictReader(io.StringIO(raw))
        if "ra1" not in (table.fieldnames or []):
            gaps.append([start,"Malformed ZTF response: " + raw[:250]])
            continue
        completed += 1
        print("ZTF window",start,group[0][0],group[-1][0],flush=True)
        for row in table:
            jd = float(row["obsjd"])
            i = min(max(0,bisect_right(dates,jd)-1),len(pts)-2)
            j0,r0,d0=pts[i]
            j1,r1,d1=pts[i+1]
            f=(jd-j0)/(j1-j0)
            ra=(r0+f*dx(r1,r0))%360
            dec=d0+f*(d1-d0)
            # Tangent-plane polygon (ZTF quadrant spans about one degree).
            poly=[(dx(float(row[f"ra{k}"]),ra)*math.cos(math.radians(dec)),
                   float(row[f"dec{k}"])-dec) for k in range(1,5)]
            sides=[a*y-b*x for (a,b),(x,y) in zip(poly,poly[1:]+poly[:1])]
            if not (min(sides)>=0 or max(sides)<=0):
                continue
            keys=(row['filefracday'],row['ccdid'],row['qid'],row['filtercode'])
            if keys in used_ids:
                continue
            used_ids.add(keys)
            candidates.append({**row,"pred_ra_deg":ra,"pred_dec_deg":dec})
    return {"completed_weekly_windows": completed, "windows_requested":len(starts), "errors":gaps, "grid_points":len(pts),
            "ephemeris_columns":header,"ztf_nominal_footprints":candidates}


def cadc():
    params=dict(search="bynameCADC",object="C/2025 V1",epoch1="2025-07-01",epoch2="2025-11-02",
                xyres="no",format="tsv",lang="en")
    raw=fetch("cadc.tsv","https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/cadcbin/ssos/ssosclf.pl",params)
    return {"response_excerpt":raw[:1600], "line_count":len(raw.splitlines())}


if __name__ == "__main__":
    results={}
    # --recent covers the final two weekly windows, including published precoveries.
    recent = "--recent" in sys.argv[2:]
    for name,job in (("ztf",lambda: ztf(2 if recent else None)),("cadc",cadc)):
        try:
            results[name]=job()
        except Exception as exc:
            results[name]={"error":str(exc)}
        print(name,json.dumps(results[name])[:2800],flush=True)
    (OUT/"archive-results.json").write_text(json.dumps(results,indent=2)+"\n")
