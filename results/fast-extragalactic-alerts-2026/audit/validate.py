"""Bounded, independently receipted 32-position Gaia, PS1 and TNS check.

Uses the frozen pair tables; never reruns discovery or expands its time window.
Unavailable services are marked UNKNOWN, never treated as empty cones.
"""
import csv
import gzip
import hashlib
import io
import json
import math
import re
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(__file__).resolve().parent
GAIA = "https://gea.esac.esa.int/tap-server/tap/sync"
PS1 = "https://catalogs.mast.stsci.edu/api/v0.1/panstarrs/dr2/mean.csv"
TNS = "https://www.wis-tns.org/sites/default/files/tns_public_objects/tns_public_objects.csv.zip"
ZTF_FIXED = {
    "ZTF19acymkzw": ("stellar flare/variable", "Gaia parallax 2.716 ± 0.014 mas; bright reference star"),
    "ZTF26acetuvm": ("stellar flare/variable", "SIMBAD eclipsing binary TYC 5341-65-1; Gaia G=11.53"),
    "ZTF26acfmvvf": ("stellar flare/variable", "Probable Galactic outburst: PS1 point source; CV/nova classifier 0.975; type unconfirmed"),
    "ZTF26acfplkt": ("bogus", "Dipole/subtraction at two parallax stars; older residuals"),
    "ZTF26acfwmzf": ("stellar flare/variable", "Rapid flare on a red PS1 point source; subtype unconfirmed"),
    "ZTF26acgcmxl": ("asteroid", "ZTF MPC 18150 and six Horizons-resolved positions"),
    "ZTF26acgrjsd": ("asteroid", "Second broker ID of asteroid (18150)"),
    "ZTF26acgwnuz": ("stellar flare/variable", "Gaia red star G=17.86, parallax 2.646 ± 0.139 mas; rapid flare"),
    "ZTF26acheiuo": ("bogus", "Bright-star dipole, DRB=0.001; older subtraction residuals"),
}
NOW = datetime.now(timezone.utc).isoformat()
receipts = {"retrieved_utc": NOW, "method": "Frozen pairs; 10 arcsec Gaia/PS1 cones, TNS full public CSV; classify matches <=1.0 arcsec; catalogues may be incomplete", "inputs": {}, "http": []}

def sha(data):
    return hashlib.sha256(data).hexdigest()

def query(label, url, timeout=40):
    entry = {"label": label, "url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat()}
    try:
        req = Request(url, headers={"User-Agent": "Open-Space-public-archive-audit/1.0 (astronomy research)", "Accept": "*/*"})
        with urlopen(req, timeout=timeout) as response:
            data = response.read()
            entry.update(status=response.status, bytes=len(data), sha256=sha(data), content_type=response.headers.get("Content-Type"))
        return data
    except Exception as err:
        entry.update(status=getattr(err, "code", None), error=str(err)[:300])
        return None
    finally:
        receipts["http"].append(entry)

def rows(path):
    p = ROOT / path
    data = p.read_bytes()
    receipts["inputs"][path] = {"sha256": sha(data), "bytes": len(data)}
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))

def sep(a, b, c, d):
    a,b,c,d = map(math.radians,(a,b,c,d))
    return 2 * math.degrees(math.asin(min(1, math.sqrt(math.sin((b-d)/2)**2 +
                            math.cos(b)*math.cos(d)*math.sin((a-c)/2)**2)))) * 3600

def f(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None

def nearest(data, a, b, xkey="ra", ykey="dec"):
    cand = [(sep(a,b,f(r.get(xkey)),f(r.get(ykey))),r)
            for r in data if f(r.get(xkey)) is not None and f(r.get(ykey)) is not None]
    return min(cand, key=lambda x: x[0]) if cand else (None, None)

items = {}
for r in rows("rubin/pairs.csv"):
    oid = r["oid"]
    if oid not in items:
        items[oid] = dict(oid=oid, survey="Rubin", ra=f(r["ra"]), dec=f(r["dec"]),
                          flags=set(), first_t=f(r["t1"])-37/86400)
    if r["flags"].strip():
        items[oid]["flags"].add(r["flags"].strip())
for r in rows("ztf/derived/ztf_pairs.csv"):
    oid = r["oid"]
    if oid not in items:
        items[oid] = dict(oid=oid, survey="ZTF", ra=f(r["meanra"]), dec=f(r["meandec"]),
                          flags=set(), first_t=f(r["mjd1"]))
assert len(items) == 32 and sum(x["survey"]=="Rubin" for x in items.values())==23
assert len(ZTF_FIXED) == 9

def get_catalog(item):
    oid = item["oid"]
    ra,dec = item["ra"],item["dec"]
    sql = ("SELECT source_id,ra,dec,parallax,parallax_error,phot_g_mean_mag,bp_rp,ruwe "
           "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),"
           "CIRCLE('ICRS',%.9f,%.9f,0.002777778))" % (ra,dec))
    gaia_url = GAIA + "?" + urlencode(dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=sql))
    body = query("gaia_" + oid, gaia_url)
    gaia = None
    if body is not None:
        try:
            gaia = list(csv.DictReader(io.StringIO(body.decode("utf-8-sig"))))
            if gaia and "ra" not in gaia[0]:
                gaia = None
        except Exception:
            gaia = None
    ps1 = None
    if dec > -30:
        params = dict(ra=ra, dec=dec, radius=10/3600,
                      columns="[objID,raMean,decMean,nDetections,gMeanPSFMag,rMeanPSFMag,iMeanPSFMag,zMeanPSFMag,iMeanKronMag]")
        body = query("ps1_"+oid, PS1+"?"+urlencode(params))
        if body is not None:
            try:
                ps1 = list(csv.DictReader(io.StringIO(body.decode("utf-8-sig"))))
                if ps1 and "raMean" not in ps1[0]:
                    ps1 = None
            except Exception:
                ps1 = None
    return oid, gaia, ps1

catalogs = {}
with ThreadPoolExecutor(max_workers=3) as pool:
    futures = [pool.submit(get_catalog, item) for item in items.values()]
    for fut in as_completed(futures):
        oid, gaia, ps1 = fut.result()
        catalogs[oid] = (gaia, ps1)

tns_data = query("tns_public_full", TNS, timeout=80)
tns = None
if tns_data is not None:
    try:
        archive = zipfile.ZipFile(io.BytesIO(tns_data))
        names = [n for n in archive.namelist() if n.lower().endswith(".csv")]
        with archive.open(names[0]) as z:
            tns = list(csv.DictReader(io.TextIOWrapper(z, encoding="utf-8-sig", errors="replace")))
        if len(tns) < 50000:
            raise ValueError("Unexpectedly short full TNS CSV")
        receipts["tns_rows"] = len(tns)
        receipts["tns_columns"] = list(tns[0])
    except Exception as err:
        receipts["tns_parse_error"] = str(err)[:250]
        tns = None

def convert(v, is_ra=False):
    v = str(v or "").strip()
    x = f(v)
    if x is not None:
        return x
    nums = [f(t) for t in re.split(r"[:\s]+", v) if t]
    if len(nums) < 3 or any(n is None for n in nums[:3]):
        return None
    sign = -1 if v.startswith("-") else 1
    return sign * (abs(nums[0])+nums[1]/60+nums[2]/3600)*(15 if is_ra else 1)

tns_sky = []
if tns is not None:
    def pick(r, *names):
        d = {re.sub("[^a-z0-9]", "", k.lower()):v for k,v in r.items() if k}
        return next((d.get(re.sub("[^a-z0-9]", "", n.lower())) for n in names
                     if d.get(re.sub("[^a-z0-9]", "", n.lower()))), None)
    for r in tns:
        ra,dec = convert(pick(r,"ra","radeg"),True),convert(pick(r,"declination","dec","decl","decdeg"))
        if ra is None or dec is None:
            continue
        name = " ".join(str(pick(r,"nameprefix","prefix") or "AT").split())+" "+str(pick(r,"name","objname") or "")
        tns_sky.append((ra,dec,name.strip()))
    receipts["tns_positions_parsed"] = len(tns_sky)
    if len(tns_sky) < 50000:
        receipts["tns_parse_error"] = "Fewer than 50,000 usable TNS coordinates"
        tns_sky = []
else:
    receipts["tns_positions_parsed"] = 0

out = []
for oid,item in sorted(items.items(),key=lambda kv: (kv[1]["survey"],kv[0])):
    gaia,ps1 = catalogs[oid]
    gsep,g = nearest(gaia or [],item["ra"],item["dec"])
    psep,p = nearest(ps1 or [],item["ra"],item["dec"],"raMean","decMean")
    sig = (f(g.get("parallax"))/f(g.get("parallax_error")) if g and f(g.get("parallax")) is not None
           and f(g.get("parallax_error")) is not None and f(g.get("parallax_error")) > 0 else None)
    point = (p and f(p.get("iMeanPSFMag")) is not None and f(p.get("iMeanKronMag")) is not None
             and abs(f(p["iMeanPSFMag"])-f(p["iMeanKronMag"])) < 0.1
             and (f(p.get("nDetections")) or 0)>=3)
    tmatch = min(((sep(item["ra"],item["dec"],ra,dec),n) for ra,dec,n in tns_sky),
                 default=(None,""), key=lambda x:x[0])
    match = tmatch[1] if tmatch[0] is not None and tmatch[0]<=10 else ""
    if item["survey"]=="ZTF":
        disposition, evidence = ZTF_FIXED[oid]
    elif "isdipole" in item["flags"]:
        disposition,evidence="bogus","At least one frozen passing pair has isdipole quality flag; source nature not established"
    elif oid=="171198361324486700":
        disposition,evidence="pre-existing point source","Legacy Surveys source 0.11 arcsec away; direct-image fade only 0.67 sigma in independent later check"
    elif gsep is not None and gsep <= 1 and sig is not None and sig >= 5:
        disposition,evidence="stellar flare/variable","Gaia DR3 <=1 arcsec with parallax S/N >=5; physical variability type not established"
    elif (gsep is not None and gsep<=1) or (psep is not None and psep<=1 and point):
        disposition,evidence="pre-existing point source","Prior Gaia/PS1 point-like counterpart <=1 arcsec; no extragalactic distance established"
    else:
        disposition,evidence="unresolved","No secure class from frozen pairs and available independent positional checks; not an extragalactic validation"
    row = dict(oid=oid,survey=item["survey"],ra=round(item["ra"],7),dec=round(item["dec"],7),
       disposition=disposition,evidence=evidence, flags=";".join(sorted(item["flags"])),
       gaia_query="ok" if gaia is not None else "unavailable",gaia_sep_arcsec=round(gsep,3) if gsep is not None else "",
       gaia_parallax_snr=round(sig,2) if sig is not None else "",
       gaia_G=g.get("phot_g_mean_mag","") if g else "",
       ps1_query="not_in_footprint" if item["dec"] <= -30 else "ok" if ps1 is not None else "unavailable",
       ps1_sep_arcsec=round(psep,3) if psep is not None else "",
       ps1_pointlike="yes" if point else "unknown" if p is None else "no",
       tns_query="full_public_csv" if tns_sky else "unavailable",
       tns_match_within_10arcsec=match,tns_sep_arcsec=round(tmatch[0],3) if match else "",
       gcn="No supported match in 234 scoped GCN comparisons (coincidences/sky_time_matches.csv)")
    out.append(row)

with (OUT/"dispositions.csv").open("w",newline="") as fh:
    writer=csv.DictWriter(fh,fieldnames=list(out[0]))
    writer.writeheader()
    writer.writerows(out)
receipts["output"]={"dispositions.csv":{"sha256":sha((OUT/"dispositions.csv").read_bytes()),"rows":len(out)}}
(OUT/"query_receipts.json").write_text(json.dumps(receipts,indent=2,sort_keys=True)+"\n")
from collections import Counter
print("32-position audit:",dict(Counter(x["disposition"] for x in out)),
      "TNS rows",len(tns_sky),"catalog successes",sum(x["gaia_query"]=="ok" for x in out),
      sum(x["ps1_query"]=="ok" for x in out))
