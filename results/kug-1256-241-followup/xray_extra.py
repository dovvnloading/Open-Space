"""Independent archive checks for dataset timing, individual observations and limits."""
import csv
import gzip
import inspect
import io
import json
import math
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.parse import urljoin
from xml.etree import ElementTree as ET

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io.votable import parse
from astropy.time import Time
from astropy.wcs import WCS

import xray as xr
from common import HOST_RA, HOST_DEC, OUT, event, fetch, status

results=[]
dates=[]
NASA=xr.NASA
XSA="https://nxsa.esac.esa.int/tap-server/tap/sync"


def dump(name, value):
    (OUT/("xray-extra-"+name+".json")).write_text(json.dumps(value,indent=2,default=str)+"\n")


def norm(v):
    return re.sub(r"[^a-z0-9]","",str(v).lower())


def report(arc, ident, start, end, f, note, kind="detection", error=None):
    results.append(dict(archive=arc,identifier=ident,mjd=start,end_mjd=end,
                        flux=f,error=error if error is not None else "",kind=kind,notes=note))
    event(arc,xr.ENERGY,start,f,error if error is not None else "",xr.UNIT,kind,str(ident),note,end_mjd=end)


def dates_from_table(tag, rows):
    for original in rows:
        r={norm(k):v for k,v in original.items()}
        name=xr.first(r,"iauname","name","srcname","erass1")
        if not name or "125851" not in str(name).replace(" ",""):
            continue
        lo=xr.mjd(xr.first(r,"mjdmin","mjdstart","mjdbeg","mjdbegin","tstart"))
        hi=xr.mjd(xr.first(r,"mjdmax","mjdend","mjdstop","tstop"))
        mid=xr.mjd(xr.first(r,"mjd","mjdobs","mjdmean"))
        dates.append(dict(table=tag,name=name,mjd_min=lo,mjd=mid,mjd_max=hi,
                          flux_fields={k:v for k,v in r.items() if k in ("mlflux0","emlflux0","mlfluxerr0","mlflux1","emlflux1")}))
        if lo!="" or mid!="":
            event("eROSITA timing "+tag,"X-ray",lo if lo!="" else mid,kind="catalogue date cross-check",
                  end_mjd=hi,source=str(name),note="Independent cone query; duplicate timing for the same survey detection.")
            status("eROSITA timing "+tag,"measured",f"{name} MJD min/mean/max {lo}/{mid}/{hi}")


def erass():
    failures=[]
    def vizier(cat):
        for host in ("vizier.cfa.harvard.edu","vizier.cds.unistra.fr"):
            try:
                raw=fetch("xr-extra-viz-"+cat.replace("/","-")+"-"+host,
                          "https://"+host+"/viz-bin/votable",params={"-source":cat,
                          "-c":f"{HOST_RA} {HOST_DEC}","-c.rs":9,
                          "-out.all":"","-out.max":50},timeout=24)
                tables=list(parse(io.BytesIO(raw)).iter_tables())
                for v in tables:
                    t=v.to_table(use_names_over_ids=True)
                    rows=[{k:xr.clean(r[k]) for k in t.colnames} for r in t]
                    if rows:
                        dump(cat.replace("/","-")+"-"+re.sub(r"[^a-zA-Z0-9]","-",str(v.name)), rows)
                        dates_from_table(str(v.name),rows)
                if any(len(t.array) for t in tables):
                    return
            except Exception as exc:
                failures.append((host,cat,str(exc)[:180]))
    vizier("J/A+A/682/A34")
    vizier("J/A+A/712/A171")
    try:
        names=xr.tap("extra-erass-tables",NASA,
                   "SELECT table_name FROM TAP_SCHEMA.tables WHERE table_name LIKE '%erass1%'")
        for item in names:
            name=item["table_name"]
            if "main" not in name.lower() or not re.fullmatch(r"[A-Za-z0-9_.]+",name):
                continue
            rows=xr.tap("extra-"+name,NASA,f"SELECT * FROM {name} WHERE "+xr.cone("ra","dec",10))
            dates_from_table(name,rows)
            dump("erass-heasarc",rows)
    except Exception as exc:
        failures.append(("HEASARC","eRASS1",str(exc)[:180]))
    dump("erass-dates",dict(rows=dates,unavailable=failures))
    if not any(d["mjd"]!="" or d["mjd_min"]!="" for d in dates):
        status("Independent eROSITA dates","unavailable",f"No time fields recovered. {failures}")


def rass_epoch():
    # 2RXS source-level entry, not the sequence nearest the sky position.
    path=xr.PREV/"rass_hardness/rass2rxs_native.csv"
    from common import input_file
    with input_file(path).open() as handle:
        row=next(r for r in csv.DictReader(handle) if r["name"]=="2RXS J125851.4+235533")
    prefix="rosat/data/pspc/processed_data/rass/release/rs"+row["skyfield_number"]
    base="https://nasa-heasarc.s3.amazonaws.com/"
    index=ET.fromstring(fetch("xr-rass-field-inventory",base,params={"list-type":"2","prefix":prefix},timeout=25))
    files=[a.text for a in index.findall(".//{*}Key") if a.text]
    names=[p for p in files if re.search(r"_(bas|events|evt|evt2)\.fits?(\.gz|\.Z)?$",p,re.I)]
    if not names:
        raise ValueError(f"No verified archive event list for field {prefix}; {len(files)} objects")
    packed=fetch("xr-ROSAT-field-photons",base+names[0],timeout=42)
    if packed[:2]==b"\x1f\x9d":
        packed=subprocess.run(["gzip","-dc"],input=packed,capture_output=True,check=True,timeout=45).stdout
    if packed[:2]==b"\x1f\x8b":
        packed=gzip.decompress(packed)
    with fits.open(io.BytesIO(packed),memmap=False) as hdus:
        h=next(h for h in hdus if "TIME" in (getattr(getattr(h,"columns",None),"names",[]) or [])
               and ("X" in h.columns.names or "RA" in h.columns.names))
        cols=list(h.columns.names)
        if {"RA","DEC"}<=set(cols):
            pos=SkyCoord(h.data["RA"]*u.deg,h.data["DEC"]*u.deg)
            inside=pos.separation(xr.target).arcsec<=90
        else:
            xi,yi=cols.index("X")+1,cols.index("Y")+1
            wh=fits.Header()
            for i,k in enumerate((xi,yi),1):
                for std,table in (("CTYPE","TCTYP"),("CRVAL","TCRVL"),("CRPIX","TCRPX"),("CDELT","TCDLT"),("CUNIT","TCUNI")):
                    if table+str(k) in h.header:
                        wh[std+str(i)]=h.header[table+str(k)]
            if not all(x in wh for x in ("CTYPE1","CTYPE2","CRVAL1","CRVAL2","CDELT1","CDELT2")):
                # US ROSAT RDF uses these REFX* / REFY* image reference keywords.
                for i,axis in enumerate(("X","Y"),1):
                    for std,src in (("CTYPE","CTYP"),("CRVAL","CRVL"),("CRPIX","CRPX"),("CDELT","CDLT")):
                        k="REF"+axis+src
                        if k in h.header:
                            wh[std+str(i)]=h.header[k]
            w=WCS(wh)
            if not w.has_celestial:
                raise ValueError("No reliable celestial coordinates for PSPC aperture")
            x,y=w.all_world2pix(HOST_RA,HOST_DEC,1)
            scale=np.abs(np.asarray(w.wcs.cdelt))*3600
            near=(np.abs(h.data["X"]-x)*scale[0]<92)&(np.abs(h.data["Y"]-y)*scale[1]<92)
            ids=np.flatnonzero(near)
            a,b=w.all_pix2world(h.data["X"][ids],h.data["Y"][ids],1)
            inside=ids[SkyCoord(a,b,unit="deg").separation(xr.target).arcsec<=90]
        ts=np.asarray(h.data["TIME"][inside],float)
        headers={**dict(hdus[0].header),**dict(h.header)}
        ref=xr.number(headers.get("MJDREF"))
        if ref is None and "MJDREFI" in headers:
            ref=float(headers["MJDREFI"])+float(headers.get("MJDREFF",0))
        if ref is None:
            raise ValueError("No documented event time reference")
        k=cols.index("TIME")+1
        unit=h.header.get("TUNIT"+str(k),headers.get("TIMEUNIT","s"))
        ts=ref+(ts+float(headers.get("TIMEZERO",0)))*u.Unit(unit).to(u.day)
        ts=ts[np.isfinite(ts)]
        if len(ts)<10 or not 47900<min(ts)<=max(ts)<49000:
            raise ValueError("Unverified source position/timestamps in ROSAT aperture")
        payload=dict(field=names[0],source=row["name"],n_aperture=len(ts),radius_arcsec=90,
                     first_mjd=float(min(ts)),last_mjd=float(max(ts)),mjdref=ref,unit=unit,
                     time_system=headers.get("TIMESYS","unspecified"),median_mjd=float(np.median(ts)))
        dump("rass-timing",payload)
        event("ROSAT native photon times","X-ray",payload["first_mjd"],kind="aperture photon dates",
              source=row["name"],end_mjd=payload["last_mjd"],
              note=f"{len(ts)} detector events within 90 arcsec, including background. ROSAT time system {payload['time_system']}; timing precision one day.")
        status("ROSAT native timing","measured",f"{len(ts)} aperture event times in RASS field "+row["skyfield_number"])


def uk_upper_limits():
    """UKSSDC archive service performs retrospective XRT photometry."""
    import pandas as pd
    from swifttools.ukssdc.data import SXPS
    from swifttools.ukssdc import main as uk
    import requests

    class Client:
        def __init__(self):
            self.n=0

        def __getattr__(self,k):
            return getattr(requests,k)

        def post(self,url,**kwargs):
            kwargs.setdefault("timeout",35)
            r=requests.post(url,**kwargs)
            self.n+=1
            xr.archived("extra-swift-ul-response-"+str(self.n),r.content)
            xr.archived("extra-swift-ul-request-"+str(self.n),dict(url=url,parameters=kwargs.get("json",kwargs.get("data",{}))))
            return r

    if hasattr(uk,"requests"):
        uk.requests=Client()
    func=SXPS.getUpperLimits
    available=inspect.signature(func).parameters
    options=dict(ra=HOST_RA,dec=HOST_DEC,bands=["total"],whichData="all",timeFormat="MJD",
                 sigma=3,doDet=True,timeout=150,waitTime=10,pandas=True,silent=True)
    product=func(**{k:v for k,v in options.items() if k in available})
    if not isinstance(product,dict):
        raise ValueError("Unexpected upper-limit return type "+str(type(product)))
    def walk(v, path):
        if isinstance(v,pd.DataFrame):
            yield path,v
        elif isinstance(v,list) and v and all(isinstance(i,dict) for i in v):
            yield path,pd.DataFrame(v)
        elif isinstance(v,dict):
            for k,item in v.items():
                yield from walk(item,path+"-"+str(k))
    n=0
    prior=json.loads((xr.PREV/"rass_hardness/swift2sxps_review.json").read_text())
    ecf=prior["measurement_fields"]["powcf_flux_a"]["value"]
    factor=xr.scale(.3,10)
    for path,t in walk(product,"xrt"):
        t.to_csv(OUT/("xray-extra-"+re.sub(r"[^a-zA-Z0-9_-]","-",path)+".csv"),index=False)
        for original in t.to_dict("records"):
            r={norm(k):v for k,v in original.items()}
            def val(*k):
                return next((r[a] for a in k if a in r and pd.notna(r[a])),None)
            rate=xr.number(val("totalupperlimit","upperlimittotal","totalratelimit","totalul","upperlimit"))
            if rate is None or rate<=0:
                continue
            # Restrict to a named total band; an unlabelled table can include multiple bands.
            band=val("band","energyband")
            if band is not None and str(band).lower() not in ("total","t","full","0.3-10"):
                continue
            star=xr.mjd(val("starttime","tstart","mjdstart","startmjd","mjd","t0","start"))
            end=xr.mjd(val("stoptime","endtime","stopmjd","endmjd","mjdstop","t1","end"))
            ident=val("obsid","datasetid","dataset","targetid","imageid") or path
            report("Swift LSXPS position limit",ident,star,end,rate*ecf*factor,
                   f"UKSSDC 3-sigma total-band upper bound {rate:.5g} count/s; "
                   f"projected with source 2SXPS ECF {ecf:.4g} and band factor {factor:.4g}. "
                   "A limit is an upper bound even if source emission is detected.",kind="upper limit")
            n+=1
    if not n:
        raise ValueError("Upper-limit endpoint returned no recognized total-band rows. Product keys: "+repr(list(product))[:250])
    status("UKSSDC position limits","measured",f"{n} calibrated XRT flux upper limits. See raw response receipts.")


def bare_heasarc(table,radius):
    rows=xr.tap("extra-bare-"+table,NASA,f"SELECT * FROM {table} WHERE "+xr.cone("ra","dec",radius))
    dump(table,rows)
    ids=[]
    for r in rows:
        start=xr.mjd(xr.first(r,"time","mjd_start","start_time","mjd"))
        end=xr.mjd(xr.first(r,"end_time","mjd_stop","stop_time"))
        ident=xr.first(r,"obsid","obs_id","sequence_number","name")
        ids.append(str(ident))
        if table=="xmmssc":
            if xr.sep(r)>6:
                continue
            flux=[xr.number(xr.first(r,f"ep_{i}_flux",f"ep{i}_flux",f"flux_band_{i}",
                           f"flux_b{i}",f"flux_{i}")) for i in (1,2,3)]
            err=[xr.number(xr.first(r,f"ep_{i}_flux_err",f"ep{i}_flux_err",f"flux_band_{i}_error",
                          f"flux_b{i}_error",f"flux_{i}_error")) for i in (1,2,3)]
            if None not in flux:
                report("Pointed XMM EPIC",ident,start,end,sum(flux),"EPIC catalogue observed bands 1+2+3.",
                       error=math.hypot(*err) if None not in err else None)
        else:
            event("HEASARC "+table,"X-ray coverage",start,kind="pointing candidate",
                  end_mjd=end,source=str(ident),note=f"{xr.sep(r)/60:.2f} arcmin off centre. Footprint unverified.")
    status("HEASARC "+table+" direct","coverage_only" if rows else "no_match",f"{len(rows)} rows within {radius/60:.2f} arcmin; {ids[:30]}")


def rapidxmm():
    tables=xr.tap("extra-rapidxmm-list",XSA,
          "SELECT table_name,description FROM TAP_SCHEMA.tables WHERE table_name LIKE '%rapid%'")
    notes=[]
    tables=sorted(tables,key=lambda r:("pixel" not in str(r["table_name"]).lower(),len(str(r["table_name"]))))
    for item in tables:
        table=item["table_name"]
        if not re.fullmatch(r"[a-zA-Z0-9_.]+",table):
            continue
        cols=xr.tap("extra-rapid-"+table+"-cols",XSA,
            f"SELECT column_name,unit,description FROM TAP_SCHEMA.columns WHERE table_name='{table}'")
        definitions={r["column_name"].lower():r for r in cols}
        dump("rapid-schema-"+table,dict(table=item,columns=cols))
        x=next((k for k in ("ra","ra_deg","ra_centre","raj2000") if k in definitions),None)
        y=next((k for k in ("dec","dec_deg","dec_centre","dej2000") if k in definitions),None)
        if not x or not y:
            continue
        nside=None
        nesting=None
        descriptions=" ".join([str(item.get("description",""))]+[str(r.get("description","")) for r in cols])
        # A geometrical cone around the optical nucleus catches the centres
        # of both pointed and slew pixels. It does not assume an undocumented order.
        rows=xr.tap("extra-rapid-"+table,XSA,f"SELECT * FROM {table} WHERE "+xr.cone(x,y,12))
        dump("rapid-values-"+table,rows)
        notes.append((table,len(rows)))
        if not rows:
            continue
        obskey=next((k for k in ("obs_id","obsid","observation_id") if k in definitions),None)
        # A catalogue rate limit cannot be converted by a universal count-rate
        # factor when the camera, filter, energy band or sky pixel is unknown.
        # Store every native bound at its angular offset, including all cells.
        for r in rows:
            offset=xr.sep(r,x,y)
            candidates=[k for k in definitions if re.search(r"ul.*3|3.*ul|upper.*3|3.*upper",k)]
            for k in candidates:
                value=xr.number(r[k])
                unit=definitions[k].get("unit") or ""
                if value is None or value<=0 or not unit:
                    continue
                ident=str(r.get(obskey,""))
                band=xr.first(r,"band","energy_band","eband","band_id")
                event("RapidXMM nearby pixel",f"native XMM band {band}",xr.mjd(xr.first(r,"mjd","mjd_start","time")),
                      value,unit=unit,kind="map limit at pixel centre",source=ident,
                      note=f"{offset:.2f} arcsec from nucleus. {k}. No assigned source-specific common flux without matched camera/filter ECF.")
    status("RapidXMM exposure and sensitivity","coverage_only" if any(n for t,n in notes) else (
           "no_match" if notes else "unavailable"),f"Tables with usable coordinates and row counts at 12 arcsec: {notes}")


def direct_csc():
    # Native CSC per-observation rows join through the unique stack/region key.
    # CXC TAP deployments expose different schemas. Try the published query
    # directly on both release endpoints, independent of service-wide metadata.
    schemas=(("https://cda.cfa.harvard.edu/csc21tap/sync","csc21"),
             ("https://cda.cfa.harvard.edu/csc2tap/sync","csc2"),
             ("https://cda.harvard.edu/csc2tap/sync","csc2"))
    failures=[]
    for index,(url,schema) in enumerate(schemas):
        try:
            try:
                master=xr.tap("extra-csc-master-by-name-"+str(index),url,
                    f"SELECT name,ra,dec,flux_aper_s,flux_aper_m,flux_aper_b FROM {schema}.master_source WHERE name='2CXO J125851.4+235526'")
                dump("csc-master-exact-"+str(index),master)
                if master:
                    add_csc_master(master[0],schema)
            except Exception as exc:
                failures.append((schema+" master",str(exc)[:160]))
            stacks=xr.tap("extra-csc-stack-"+str(index),url,
                f"SELECT * FROM {schema}.master_stack_assoc WHERE name='2CXO J125851.4+235526'")
            if not stacks:
                continue
            dump("csc-exact-associations",stacks)
            obsrows=[]
            for stack in stacks:
                sid=str(stack.get("detect_stack_id",""))
                rid=xr.number(stack.get("region_id"))
                if not re.fullmatch(r"[\w.+:-]+",sid) or rid is None:
                    continue
                for name in ("observation_source","per_observation","obs_source"):
                    try:
                        rows=xr.tap("extra-csc-"+str(index)+"-"+name+"-"+sid,url,
                          f"SELECT * FROM {schema}.{name} WHERE detect_stack_id='{sid}' AND region_id={int(rid)}")
                        obsrows+=rows
                        if rows:
                            break
                    except Exception as exc:
                        failures.append((name,str(exc)[:200]))
            dump("csc-observation-rows",obsrows)
            n=0
            for r in obsrows:
                ident=xr.first(r,"obsid","obs_id","detect_stack_id")
                start=xr.mjd(xr.first(r,"gti_mjd_obs","mjd_obs","gti_start","date_obs","obs_start"))
                end=xr.mjd(xr.first(r,"gti_mjd_end","mjd_end","gti_end","date_end","obs_end"))
                for tag,bounds in (("soft+medium",(.5,2)),("b",(.5,7)),("s",(.5,1.2))):
                    names=("s","m") if tag=="soft+medium" else (tag,)
                    v=[xr.number(r.get("flux_aper_"+b)) for b in names]
                    if any(x is None or x<=0 for x in v):
                        continue
                    factor=xr.scale(*bounds)
                    lo=[xr.number(r.get("flux_aper_lolim_"+b)) for b in names]
                    hi=[xr.number(r.get("flux_aper_hilim_"+b)) for b in names]
                    err=max(sum(v)-sum(lo),sum(hi)-sum(v))*factor if None not in lo+hi else None
                    report("Chandra individual exposure",ident,start,end,sum(v)*factor,
                           f"Exact CSC associated stack/region; aperture {tag} {bounds} keV × {factor:.4g}.",error=err)
                    n+=1
                    break
            if n:
                status("CSC native per-observation rows","measured",f"{n} individual fluxes")
                return
        except Exception as exc:
            failures.append((url,str(exc)[:200]))
    status("CSC native per-observation rows","unavailable",str(failures)[:420])


def add_csc_master(row,table):
    r={norm(k):v for k,v in row.items()}
    # HEASARC and VizieR projections can use different field spellings.
    f=[xr.number(xr.first(r,"fluxaper"+b,"aperflux"+b,"flux"+b,"apflux"+b)) for b in ("s","m","b")]
    if all(x is not None and x>0 for x in f[:2]):
        bands=(.5,2)
        fl=f[0]+f[1]
    elif f[2] is not None and f[2]>0:
        bands=(.5,7)
        fl=f[2]
    else:
        return False
    if not 1e-20<fl<1e-7:
        return False
    fac=xr.scale(*bands)
    report("Chandra master-source catalogue","2CXO J125851.4+235526","","",fl*fac,
           f"Undated master from {table}. Observed native {bands} keV flux {fl:.5g} times {fac:.5g}.")
    status("CSC "+table,"measured","Exact master source energy flux. Undated combination.")
    return True


def csc_mirrors():
    names=xr.tap("extra-csc-heasarc-schema",NASA,"SELECT table_name,description FROM TAP_SCHEMA.tables "
                 "WHERE table_name LIKE '%csc%' OR table_name LIKE '%chandra%'")
    for r in names:
        name=r["table_name"]
        if not re.fullmatch(r"[A-Za-z0-9_.]+",name):
            continue
        if not re.search(r"csc.*(2|master)|chandrasrccat",name.lower()):
            continue
        def check():
            rows=xr.tap("extra-csc-nasa-source-"+name,NASA,f"SELECT * FROM {name} WHERE "+xr.cone("ra","dec",4))
            dump("csc-heasarc-"+name,rows)
            for rec in rows:
                if xr.sep(rec)<2.5:
                    add_csc_master(rec,name)
        xr.guard("CSC HEASARC "+name,check)


def corrected_curves():
    # Repeat the SDK requests using its documented array of bands and the
    # run-one query selections. Keep all response bodies independently.
    from swifttools.ukssdc.data import SXPS
    import pandas as pd
    prior=json.loads((xr.PREV/"rass_hardness/swift2sxps_review.json").read_text())
    ecf=prior["measurement_fields"]["powcf_flux_a"]["value"]
    path=OUT/"xray-inputs/swift-LSXPS-cone.gz"
    selections=[("2SXPS",{"sourceName":prior["returned_name"]})]
    if path.is_file():
        with gzip.open(path,"rt") as src:
            obs=list(csv.DictReader(src))
        if len(obs)==1:
            record={norm(k):v for k,v in obs[0].items()}
            lid=xr.first(record,"lsxpsid","sourceid")
            if lid:
                selections.append(("LSXPS",{"sourceID":int(lid)}))
    def one(cat,kw):
        err=""
        for bands in (["total"],["Total"],"all"):
            try:
                p=SXPS.getLightCurves(**kw,cat=cat,bands=bands,binning="observation",timeFormat="MJD",
                                      returnData=True,saveData=False,groupRates=False,groupULs=False)
                break
            except Exception as exc:
                err=str(exc)
        else:
            raise RuntimeError(err)
        def frames(v,path):
            if isinstance(v,pd.DataFrame):
                yield path,v
            elif isinstance(v,dict):
                for k,z in v.items():
                    yield from frames(z,path+"-"+str(k))
        n=0
        for name,df in frames(p,cat):
            if bands=="all" and "total" not in name.lower():
                continue
            safe=re.sub(r"[^a-zA-Z0-9_-]","-",name)
            df.to_csv(OUT/("xray-extra-lc-"+safe+".csv"),index=False)
            for rec in df.to_dict("records"):
                r={norm(k):v for k,v in rec.items() if pd.notna(v)}
                t=xr.mjd(xr.first(r,"mjd","time","t"))
                if t=="":
                    continue
                tneg=xr.number(xr.first(r,"mjdneg","timeneg","tneg","mjdminus","timeerrneg")) or 0
                tpos=xr.number(xr.first(r,"mjdpos","timepos","tpos","mjdplus","timeerrpos")) or 0
                rate=xr.number(xr.first(r,"rate","countrate"))
                lim=xr.number(xr.first(r,"upperlimit","ul","ratelimit","upperlimitrate"))
                if lim is None and "ul" in name.lower():
                    lim=rate
                upper=lim is not None and lim>0
                val=lim if upper else rate
                if val is None or val<=0:
                    continue
                er=[xr.number(xr.first(r,k)) for k in ("ratepos","rateneg","rateerr","error")]
                err=max((abs(x) for x in er if x is not None),default=None)
                f=ecf*xr.scale(.3,10)
                ident=xr.first(r,"obsid","obsids","observationid") or str(kw)
                report("Swift XRT "+cat+" visit",ident,t-abs(tneg),t+abs(tpos),val*f,
                       f"{name}; total band {val:.4g} count/s, 2SXPS observed ECF {ecf:.4g}, common-band factor {f/ecf:.5g}.",
                       kind="upper limit" if upper else "detection",error=err*f if not upper and err is not None else None)
                n+=1
        if not n:
            raise ValueError("No measured or upper-limit total-band bins; returned keys "+str(list(p))[:220])
        status("Swift XRT "+cat+" epochs","measured",f"{n} measured or upper-limit bins")
    for cat,kw in selections:
        xr.guard(cat+" Swift corrected curve",lambda:one(cat,kw))


def run():
    import importlib.util
    spec=importlib.util.spec_from_file_location("slew_calibration",xr.BASE/"catalogs.py")
    xr.model=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(xr.model)
    # A separate script invocation has empty module buffers. No baseline rows
    # are re-emitted or included in the supplement ledger.
    tasks={"Independent eROSITA epochs":erass,"Source-specific ROSAT epochs":rass_epoch,
           "CSC exact association":direct_csc,"CSC public mirrors":csc_mirrors,"CSC VizieR":csc_vizier,
           "RapidXMM":rapidxmm,"Einstein Probe registered coverage":einstein_probe_vo,"Swift XRT position photometry":swift_both}
    for t,rad in (("chanmaster",22*60),("swiftmastr",17*60),("xmmmaster",18*60),("xmmssc",12),("rosmaster",60*60)):
        tasks[t+" native"]=lambda t=t,rad=rad:bare_heasarc(t,rad)
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(lambda item:xr.guard(*item),tasks.items()))
    with (OUT/"xray-extra-fluxes.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=["archive","identifier","mjd","end_mjd","flux","error","kind","notes"])
        w.writeheader()
        w.writerows(sorted(results,key=lambda r:(xr.mjd(r["mjd"]) or 1e7,r["archive"],str(r["identifier"]))))
    if xr.hashes:
        (OUT/"xray-extra-service-inputs.tsv").write_text("input\tsha256\n"+"".join(p+"\t"+digest+"\n" for p,digest in sorted(xr.hashes)))


def swift_both():
    # Upper-limit submission installs a request recorder shared by both SDK calls.
    xr.guard("UKSSDC position limits",uk_upper_limits)
    xr.guard("Swift per-visit product",corrected_curves)


def csc_vizier():
    for host,cat in (("vizier.cfa.harvard.edu","IX/57A"),("vizier.cds.unistra.fr","IX/57"),
                     ("vizier.cfa.harvard.edu","IX/57")):
        try:
            raw=fetch("xr-extra-CDS-CSC-"+host+"-"+cat.split("/")[-1],
                      "https://"+host+"/viz-bin/votable",params={"-source":cat,"-c":f"{HOST_RA} {HOST_DEC}",
                      "-c.rs":4,"-out.all":"","-out.max":100},timeout=22)
            count=0
            for tab in parse(io.BytesIO(raw)).iter_tables():
                t=tab.to_table(use_names_over_ids=True)
                if len(t)==0:
                    continue
                rows=[{k:xr.clean(entry[k]) for k in t.colnames} for entry in t]
                ident=re.sub(r"[^a-zA-Z0-9_-]","-",str(tab.name))
                dump("csc-cds-"+ident,rows)
                fields={f.name: str(f.description or "") for f in tab.fields}
                dump("csc-cds-fields-"+ident,fields)
                for r in rows:
                    if not any("125851.4+235526" in str(v) for v in r.values()):
                        continue
                    tagged={}
                    for k,v in r.items():
                        f=xr.number(v)
                        descr=fields.get(k,"").lower()
                        if f is None or f<=0 or "flux" not in descr or any(w in descr for w in (
                                "error","bound","limit","photon","absorbed","power-law","blackbody","model","minimum","maximum","("+"1")):
                            continue
                        if "aperture" not in descr and "aperture" not in k.lower():
                            continue
                        unit=t[k].unit
                        if unit is None or not unit.is_equivalent(u.erg/(u.cm**2*u.s)):
                            continue
                        for label,bandname in (("b","broad"),("s","soft"),("m","medium")):
                            if re.search(r"\b"+bandname+r"\b",descr):
                                tagged["fluxaper"+label]=(f*unit).to_value(u.erg/(u.cm**2*u.s))
                    if add_csc_master(tagged or r,str(tab.name)):
                        count+=1
            if count:
                status("CSC CDS mirror","measured",f"{count} exact source aperture flux rows in {cat}.")
                return
        except Exception as exc:
            status("CSC CDS "+host+" "+cat,"unavailable",str(exc)[:180])


def einstein_probe_vo():
    rows=xr.tap("extra-ep-vo","https://dc.g-vo.org/tap/sync",
        "SELECT DISTINCT res_title,standard_id,access_url FROM rr.resource NATURAL JOIN rr.capability "
        "NATURAL JOIN rr.interface WHERE 1=ivo_nocasematch(res_title,'%Einstein Probe%')")
    dump("ep-registered-services",rows)
    for i,r in enumerate(rows):
        ident="EP VO "+str(r["res_title"])
        url=r["access_url"]
        std=str(r["standard_id"]).lower()
        try:
            if "std/tap" in std:
                root=url.rstrip("/")
                endpoint=root if root.endswith("/sync") else root+"/sync"
                rows2=xr.tap("extra-ep-registry-tap-"+str(i),endpoint,
                    "SELECT TOP 100 * FROM ivoa.obscore WHERE 1=CONTAINS("
                    f"POINT('ICRS',{HOST_RA},{HOST_DEC}),s_region)")
            elif "std/sia" in std or "std/ssa" in std:
                if "#query-2" in std or "#query-aux-2" in std:
                    params={"POS":f"CIRCLE {HOST_RA} {HOST_DEC} 0.0001"}
                else:
                    params={"POS":f"{HOST_RA},{HOST_DEC}","SIZE":0.0002}
                    if "std/ssa" in std:
                        params["REQUEST"]="queryData"
                raw=fetch("extra-ep-vo-coverage-"+str(i),url,params=params,timeout=30)
                tables=list(parse(io.BytesIO(raw)).iter_tables())
                if not tables:
                    raise ValueError("No recognized VOTable")
                tab=tables[0].to_table(use_names_over_ids=True)
                rows2=[{k.lower():xr.clean(v[k]) for k in tab.colnames} for v in tab]
            else:
                continue
            dump("ep-vo-nearby-"+str(i),rows2)
            for v in rows2:
                mission=xr.first(v,"obs_collection","facility_name","instrument_name")
                if mission and not re.search(r"einstein|\bep\b|wxt|fxt",str(mission),re.I):
                    continue
                event(ident,"X-ray",xr.mjd(xr.first(v,"t_min","mjd","date_obs","time")),kind="footprint",
                      source=str(xr.first(v,"obs_id","obsid","obs_publisher_did","title") or ""),
                      end_mjd=xr.mjd(xr.first(v,"t_max","date_end")),
                      note="Registered EP service at nuclear coordinates. No count/flux limit from footprint alone.")
            status(ident,"coverage_only" if rows2 else "no_match",f"{len(rows2)} position-specific service rows; "+str(r["standard_id"]))
        except Exception as exc:
            status(ident,"unavailable",str(exc)[:250])
