"""Check external catalogues and fold available optical epochs for each trigger host."""
import argparse
import csv
import gzip
import io
import hashlib
import json
from time import sleep
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
from datetime import datetime,timezone
import numpy as np
from astropy import units as u
from astropy.coordinates import EarthLocation,SkyCoord
from astropy.time import Time
from astropy.utils import iers
iers.conf.auto_download=False
VIZ="https://vizier.cfa.harvard.edu/viz-bin/asu-tsv?"
PLANETS="https://exoplanetarchive.ipac.caltech.edu/TAP/sync?"
PAL=EarthLocation.from_geodetic(-116.863*u.deg,33.356*u.deg,1706*u.m)
def tap(url,q):return url+urlencode(dict(query=q,format="csv",lang="ADQL",request="doQuery"))
def tsv(raw):
    lines=[s for s in raw.decode(errors="replace").splitlines() if s.strip() and not s.startswith("#")]
    if not lines:return []
    return [{k:v.strip() for k,v in r.items()} for r in csv.DictReader([lines[0],*lines[3:]],delimiter="\t")]
def fetch(url,path):
    receipt=path.with_suffix(".receipt")
    if path.exists() and receipt.exists():
        meta=json.loads(receipt.read_text())
        raw=gzip.decompress(path.read_bytes())
        if meta.get("url")==url and hashlib.sha256(raw).hexdigest()==meta.get("sha256"):
            return raw,meta
    meta=dict(url=url,retrieved_utc=datetime.now(timezone.utc).isoformat(),failed_attempts=[])
    raw=b""
    attempts=6 if url.startswith(VIZ) else 3
    for attempt in range(attempts):
        # A failed response can remain cached at this VizieR mirror. Sorting
        # requests the same rows with an explicit order on a retry.
        route=url
        if attempt and url.startswith(VIZ):
            sort_key="TimeG" if "epphot" in url else "_r"
            route=url+"&"+urlencode({"-sort":sort_key})
            if attempt==2:
                route+="&-oc.form=d"
            if attempt==3:
                route=url+"&-c.eq=J2000"
            if attempt==4:
                route=url+"&-out.all=1"
            if attempt==5:
                route=url+"&-oc.form=dec"
        try:
            with urlopen(route,timeout=65) as response:
                raw=response.read()
                meta["http"]=response.status
            if not raw or b'QUERY_STATUS" value="ERROR"' in raw or b'#INFO\tError' in raw or raw.lstrip().startswith(b'<'):
                raise ValueError("service error or non-tabular response")
            meta.update(response_url=route,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw))
            path.write_bytes(gzip.compress(raw,mtime=0))
            break
        except (ValueError,OSError) as exc:
            meta["failed_attempts"].append(dict(url=route,message=str(exc),bytes=len(raw),sha256=hashlib.sha256(raw).hexdigest()))
            raw=b""
            if attempt+1<attempts:
                sleep(1+attempt)
            else:
                meta["error"]=str(exc)
    receipt.write_text(json.dumps(meta)+"\n")
    return raw,meta


def records(path):
    paths=sorted(path.glob("*.json")) if path.is_dir() else [path]
    for file in paths:
        opener=gzip.open if file.suffix==".gz" else open
        with opener(file,"rt") as stream:
            data=json.load(stream)
        yield from data if isinstance(data,list) else [data]


def historical_events(path):
    for row in records(path):
        for e in row['image_checks']:
            if e['diagnostic'] in ('individual trigger','periodic trigger'):
                yield str(row['tic']), dict(product=f"tic{row['tic']}-s{row['sector']}.npz",
                    kind='periodic' if e.get('period_days') else 'isolated',survey='historical_ffi',
                    epoch_btjd=e['epoch_btjd'],period_days=e.get('period_days'),duration_min=e['duration_min'])


def get_events(short,extended,historical=None,recent=None):
    by={}
    for path in short:
        for a in records(path):
            for kind in ("periodic","isolated"):
                for ev in a.get("search",a).get(kind,[]):
                    by.setdefault(str(a["tic"]),[]).append(dict(product=a["filename"],kind=kind,survey="spoc",**ev))
    for a in records(extended):
        for kind in ("periodic","isolated"):
            for ev in a.get("search",{}).get(kind,[]):
                by.setdefault(str(a["tic"]),[]).append(dict(product=a["lc_file"],kind=kind,survey="ffi",**ev))
    if historical:
        for tic,event in historical_events(historical):
            by.setdefault(tic,[]).append(event)
    if recent:
        for a in records(recent):
            for kind in ("periodic","isolated"):
                for e in a['search'][kind]:
                    by.setdefault(a['tic'],[]).append(dict(e,kind=kind,survey='raw_ffi',product=a['frames_sha256']))
    return by
def visit_epochs(row,cache,queries):
    tic=row["tic"]
    valid={}
    for tag,url in queries.items():
        path=cache/(tic+"-"+tag+".gz")
        native=cache/(tic+"-ztf-native.gz")
        native_receipt=native.with_suffix(".receipt")
        if tag=="ztf" and native.exists() and native_receipt.exists():
            raw=gzip.decompress(native.read_bytes())
            rec=json.loads(native_receipt.read_text())
            if hashlib.sha256(raw).hexdigest()!=rec["sha256"] or not raw:
                raise ValueError("Native ZTF photometry digest mismatch for "+tic)
            rec=dict(rec,cache_name=native.name)
        else:
            raw,rec=fetch(url,path)
        valid[tag]=(raw,rec)
    return parse_epochs(row,valid)


def parse_epochs(row,valid):
    ra,dec,ra0,dec0=map(float,(row["ra_2016"],row["dec_2016"],row["ra_2000"],row["dec_2000"]))
    zraw=valid["ztf"][0]
    if zraw:
        reader=csv.DictReader(io.StringIO(zraw.decode()))
        if not {"oid","mjd","mag","magerr","ra","dec","catflags","exptime"}.issubset(reader.fieldnames or []):
            raise ValueError("Unexpected ZTF response schema for "+row["tic"])
        zrows=list(reader)
    else:
        zrows=[]
    zsel=[]
    for r in zrows:
        year=2000+(float(r['mjd'])-51544.5)/365.25
        dra=(ra-ra0+180)%360-180
        dist=3600*np.hypot((((float(r['ra'])-ra-dra*(year-2016)/16)+180)%360-180)*np.cos(np.radians(dec)),float(r['dec'])-dec-(dec-dec0)*(year-2016)/16)
        if int(r["catflags"])==0 and dist<=2 and float(r["magerr"])>0:zsel.append(r)
    epochs={}
    if zsel:
        zt=Time([float(r["mjd"])+.5*float(r["exptime"])/86400 for r in zsel],format="mjd",scale="utc",location=PAL)
        sc=SkyCoord(ra*u.deg,dec*u.deg);t=(zt.tdb+zt.light_travel_time(sc)).jd-2457000
        flux=10**(-.4*(np.array([float(r["mag"]) for r in zsel])-13.))
        epochs["ztf"]=(t,flux,.4*np.log(10)*flux*np.array([float(r["magerr"]) for r in zsel]),np.array([r["oid"] for r in zsel]))
    gr=tsv(valid["gaiaepoch"][0])
    for band in ("G","BP","RP"):
        g=[r for r in gr if r.get(band+"rVFlag")=="0" and r.get("Time"+band) and r.get("F"+band) and r.get("e_F"+band) and (band=="G" or r.get(band+"rFlag")=="0")]
        if g:
            t=Time([float(r["Time"+band])+2455197.5 for r in g],format="jd",scale="tcb").tdb.jd-2457000
            epochs["gaia_"+band]=(t,np.array([float(r["F"+band]) for r in g]),np.array([float(r["e_F"+band]) for r in g]),np.repeat(band,len(g)))
    fields={}
    for tag in ("simbad","ps","toi","vsx","eb","atlas"):
        raw=valid[tag][0]
        fields[tag]=tsv(raw) if tag in ("vsx","eb","atlas") else list(csv.DictReader(io.StringIO(raw.decode())))
    return valid,fields,epochs
def fold(ep,event):
    t,f,s,g=ep;p=event.get("period_days");width=event["duration_min"]/1440
    phase=(t-event["epoch_btjd"]+.5*p)%p-.5*p if p else t-event["epoch_btjd"]
    inside=abs(phase)<=width/2
    ans=dict(matched_points=int(inside.sum()),points=len(t))
    if not inside.any():return ans
    # Compare within the same survey oid/band; estimate one common relative deficit.
    numerator=denom=0.
    for name in np.unique(g[inside]):
        sel=(g==name)&np.isfinite(f+s)&(s>0)
        if (sel&~inside).sum()<5:continue
        norm=np.median(f[sel&~inside])
        if norm<=0:continue
        sig=s[sel]/norm
        floor=1.4826*np.median(abs(f[sel&~inside]/norm-1));w=1/np.maximum(sig,max(floor,1e-8))**2
        a=inside[sel]; y=f[sel]/norm
        if a.any() and (~a).any():
            wa=w[a].sum();wb=w[~a].sum(); iv=1/(1/wa+1/wb)
            numerator+=iv*(np.average(y[~a],weights=w[~a])-np.average(y[a],weights=w[a]));denom+=iv
    if denom:ans.update(depth=float(numerator/denom),formal_sigma=float(denom**-.5))
    return ans
def check(item,cache,out,known):
    row,ev=item;ra,dec=float(row["ra_2016"]),float(row["dec_2016"]);ra0,dec0=float(row["ra_2000"]),float(row["dec_2000"]);tic=row["tic"]
    cone={"-c":f"{ra0:.9f} {dec0:+.9f}","-c.rs":"63","-out.max":"500"}
    qs={key:VIZ+urlencode({"-source":tab,**cone}) for key,tab in (("vsx","B/vsx/vsx"),("eb","J/ApJS/258/16/tess-ebs"),("atlas","J/AJ/156/241"))}
    qs.update(simbad=tap("https://simbad.harvard.edu/simbad/sim-tap/sync?",f"select top 150 main_id,otype,sp_type,ra,dec from basic where 1=contains(point('ICRS',ra,dec),circle('ICRS',{ra0},{dec0},0.0175))"),
              ps=tap(PLANETS,f"select pl_name,hostname,ra,dec,pl_orbper from ps where default_flag=1 and ra between {ra0-.018/max(.01,np.cos(np.radians(dec)))} and {ra0+.018/max(.01,np.cos(np.radians(dec)))} and dec between {dec0-.018} and {dec0+.018}"),
              toi=tap(PLANETS,"select tid,toi,tfopwg_disp,pl_orbper from toi where tid="+tic),
              gaiaepoch=VIZ+urlencode({"-source":"I/355/epphot","Source":row["gaia"],"-out.all":"1","-out.max":"unlimited"}),
              ztf="https://irsa.ipac.caltech.edu/cgi-bin/ZTF/nph_light_curves?"+urlencode(dict(POS=f"CIRCLE {ra} {dec} 0.006",COLLECTION="ztf_dr24",FORMAT="csv")))
    v,fields,epochs=visit_epochs(row,cache,qs)
    for r in fields["simbad"]:
        r["sep_epoch2000_as"]=3600*np.hypot((float(r["ra"])-ra0)*np.cos(np.radians(dec)),float(r["dec"])-dec0)
    if fields["simbad"]:fields["simbad"].sort(key=lambda r:r["sep_epoch2000_as"])
    wd=SkyCoord(ra*u.deg,dec*u.deg);seps=wd.separation(known[1]).arcsec;k=seps.argmin()
    result=dict(tic=tic,gaia=row["gaia"],wd=row["name"],queries={k:m for k,(_,m) in v.items()},fields=fields,
                closest_debris_monitor=dict(name=known[0][k],sep_arcsec=float(seps[k])),epochs_available={k:len(row[0]) for k,row in epochs.items()},
                events=[dict(product=e["product"],kind=e["kind"],survey=e["survey"],epoch_btjd=e["epoch_btjd"],period_days=e.get("period_days"),duration_min=e["duration_min"],optical={k:fold(data,e) for k,data in epochs.items()}) for e in ev])
    (out/(tic+".json")).write_text(json.dumps(result,allow_nan=False)+"\n")
if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    for name in ("matches","short","extension","cache","out"):p.add_argument(name,type=Path)
    p.add_argument("--extra-short",action="append",type=Path,default=[])
    p.add_argument("--historical",type=Path)
    p.add_argument("--recent",type=Path)
    p.add_argument("--tic",action="append",help="Refresh only these TICs in an existing output directory")
    a=p.parse_args();a.cache.mkdir(exist_ok=True,parents=True);a.out.mkdir(exist_ok=True,parents=True)
    url="https://janvanroestel.github.io/WD-transitingdebris-alerts/"
    path=a.cache/"debris.html"
    if not path.exists():path.write_bytes(urlopen(url,timeout=60).read())
    entries=json.JSONDecoder().raw_decode(path.read_text().split("const DATA = ")[1])[0]
    centers=SkyCoord([v["ra_sex"] for v in entries],[v["dec_sex"] for v in entries],unit=(u.hourangle,u.deg))
    known=([v["name"] for v in entries],centers)
    (a.out/"debris_list.json").write_text(json.dumps(dict(url=url,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),sources=[{k:v[k] for k in ("name","ra_sex","dec_sex","group")} for v in entries]),indent=1)+"\n")
    by=get_events([a.short,*a.extra_short],a.extension,a.historical,a.recent)
    jobs=[(r,by[r["tic"]]) for r in csv.DictReader(a.matches.open()) if r["tic"] in by and (not a.tic or r["tic"] in a.tic)]
    with ThreadPoolExecutor(6) as ex:list(ex.map(lambda item:check(item,a.cache,a.out,known),jobs))
