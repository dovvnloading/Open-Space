#!/usr/bin/env python3
"""Bounded Euclid DR1-F HATS/Parquet fallback after IRSA TAP timeouts.

All catalogue bytes are accessed on GitHub Actions via public NASA IRSA DR1 S3.
The same cones and colour thresholds are read from the prospectively committed
freeze.json. HATS partition metadata, source ETags, projected selected rows,
and processing receipts are deposited; this is not a full 1900 deg² survey.
"""
import csv,hashlib,io,json,math,os,sys,time
from xml.etree import ElementTree as ET
from urllib.parse import quote
from pathlib import Path
from collections import defaultdict
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.fs as pafs
from astropy import units as u
from astropy_healpix import HEALPix
import requests
import run as frozen

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"s3_output"
CONF=frozen.CFG
BASE="https://nasa-irsa-euclid-dr1.s3.us-east-1.amazonaws.com/"
BUCKET="nasa-irsa-euclid-dr1"
PREFIX="contributed/dr1/hats/reformatted_catalogs/"
CACHE={}
USERAGENT={"User-Agent":"Open-Space-Euclid-DR1-F-HATS-frozen-pilot/1.0"}
FIELDS=["ra","dec","object_id","tileid","vis_det","point_like_prob",
"spurious_prob","ellipticity","mumax_minus_mag","det_quality_flag",
"flag_y","flag_j","flag_h","flag_vis",
"flux_vis_3fwhm_aper","fluxerr_vis_3fwhm_aper",
"flux_y_3fwhm_aper","fluxerr_y_3fwhm_aper",
"flux_j_3fwhm_aper","fluxerr_j_3fwhm_aper",
"flux_h_3fwhm_aper","fluxerr_h_3fwhm_aper",
"flux_g_ext_decam_3fwhm_aper","fluxerr_g_ext_decam_3fwhm_aper",
"flux_r_ext_decam_3fwhm_aper","fluxerr_r_ext_decam_3fwhm_aper",
"flux_i_ext_decam_3fwhm_aper","fluxerr_i_ext_decam_3fwhm_aper",
"flux_z_ext_decam_3fwhm_aper","fluxerr_z_ext_decam_3fwhm_aper",
"flux_i_ext_panstarrs_3fwhm_aper","fluxerr_i_ext_panstarrs_3fwhm_aper",
"flux_z_ext_panstarrs_3fwhm_aper","fluxerr_z_ext_panstarrs_3fwhm_aper"]
def sha(raw):return hashlib.sha256(raw).hexdigest()
def path_for(kind):
    collection="mer_catalogue_"+kind
    return PREFIX+collection+"/euclid_dr1_reformatted_catalogs_"+collection+"-hats/"
def get(url,timeout=60):
    for i in range(3):
        try:
            r=requests.get(url,headers=USERAGENT,timeout=timeout)
            r.raise_for_status()
            return r
        except requests.RequestException:
            if i==2:raise
            time.sleep(2**i)
def healpix_pixels(order,field):
    hp=CACHE.get(order)
    if hp is None:
        hp=HEALPix(nside=2**order,order="nested",frame="icrs")
        CACHE[order]=hp
    return set(int(x) for x in hp.cone_search_lonlat(
        field["ra"]*u.deg,field["dec"]*u.deg,field["r"]*u.deg))
def inside(ra,dec,field):
    d1=np.radians(dec);d2=math.radians(field["dec"])
    dx=np.radians(ra-field["ra"])
    sn=np.sin((d1-d2)/2)**2+np.cos(d1)*math.cos(d2)*np.sin(dx/2)**2
    return sn<=math.sin(math.radians(field["r"])/2)**2
def broad(r):
    n=frozen.num
    j=n(r.get("flux_j_3fwhm_aper"));h=n(r.get("flux_h_3fwhm_aper"))
    y=n(r.get("flux_y_3fwhm_aper"))
    je=n(r.get("fluxerr_j_3fwhm_aper"))
    he=n(r.get("fluxerr_h_3fwhm_aper"))
    spur=n(r.get("spurious_prob"))
    return (n(r.get("vis_det"))==0 and j is not None and j>.65
            and h is not None and h>0 and y is not None and y<.85*j
            and j>.65*h and je is not None and je>0 and he is not None
            and he>0 and (spur is None or spur<.7))

class HTTPRangeFile(io.RawIOBase):
    """Parquet fallback: exact S3 HTTP ranges with SHA256 per byte-range."""
    def __init__(self,url,size,etag):
        self.url,self.size,self.etag,self.pos=url,size,etag,0
        self.ranges=[]
    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.pos
    def seek(self,offset,whence=0):
        if whence==0:self.pos=offset
        elif whence==1:self.pos+=offset
        elif whence==2:self.pos=self.size+offset
        else:raise ValueError(whence)
        return self.pos
    def read(self,n=-1):
        if n is None or n<0:n=self.size-self.pos
        n=min(n,self.size-self.pos)
        if n<=0:return b""
        a,b=self.pos,self.pos+n-1
        headers={"Range":"bytes=%d-%d"%(a,b),"Accept-Encoding":"identity"}
        if self.etag:headers["If-Match"]=self.etag
        rr=requests.get(self.url,headers=headers,timeout=65)
        rr.raise_for_status()
        expected="bytes %d-%d/%d"%(a,b,self.size)
        if rr.status_code!=206 or rr.headers.get("Content-Range")!=expected or len(rr.content)!=n:
            raise IOError("Unexpected S3 HTTP range response")
        self.ranges.append([a,b,sha(rr.content)])
        self.pos+=n
        return rr.content

def scalar(x):
    if isinstance(x,bytes):return x.decode("utf-8","replace")
    if isinstance(x,(np.integer,)):return int(x)
    if isinstance(x,(np.floating,)):return float(x)
    return x
def main():
    OUT.mkdir(exist_ok=True)
    raws=[];records=[];found={};controls=defaultdict(list);failures=[];field_completeness={}
    # The DR1 HATS metadata are official publisher objects, not inferred paths.
    fs=None
    try:
        proxy=os.environ.get("HTTPS_PROXY") or os.environ.get("https_proxy")
        fs=pafs.S3FileSystem(anonymous=True,region="us-east-1",
                             proxy_options=proxy or None,
                             request_timeout=60,connect_timeout=12)
    except Exception as e:
        failures.append({"stage":"S3FileSystem","error":str(e)[:400]})
    for kind in ("wide_survey","deep_survey"):
        fields=[f for f in CONF["fields"] if kind in f["table"]]
        parent=path_for(kind);meta_url=BASE+parent+"partition_info.csv"
        try:
            response=get(meta_url)
            raw=response.content
            (OUT/("partition_info_"+kind+".csv")).write_bytes(raw)
            source={"kind":kind,"uri":meta_url,"bytes":len(raw),"sha256":sha(raw),
                    "etag":response.headers.get("ETag"),"retrieved_utc":frozen.utc()}
            raws.append(source)
            items=list(csv.DictReader(io.StringIO(raw.decode("utf-8-sig"))))
            assert items and set(("Norder","Npix","Dir","num_rows")).issubset(items[0]),"HATS partition schema changed"
        except Exception as e:
            failures.append({"stage":"partition_info","kind":kind,"uri":meta_url,"error":str(e)[:400]})
            for f in fields:field_completeness[f["id"]]={"status":"METADATA_UNAVAILABLE"}
            continue
        intersect={}
        for norder in sorted({int(r["Norder"]) for r in items}):
            pixels=set()
            for f in fields:pixels|=healpix_pixels(norder,f)
            intersect[norder]=pixels
        matching=[r for r in items if int(r["Npix"]) in intersect[int(r["Norder"])]]
        matching.sort(key=lambda r:(int(r["Norder"]),int(r["Npix"])))
        original_count=len(matching)
        # Bound payload and resource use prospectively; truncation is explicit.
        maxparts=24
        matching=matching[:maxparts]
        for f in fields:
            eligible=[r for r in items if int(r["Npix"]) in healpix_pixels(int(r["Norder"]),f)]
            selected=set((int(r["Norder"]),int(r["Npix"])) for r in matching)
            fully=all((int(r["Norder"]),int(r["Npix"])) in selected for r in eligible)
            field_completeness[f["id"]]={"status":"PARTITIONS_ATTEMPTED" if fully else "TRUNCATED_PARTITIONS",
                 "n_intersecting_partitions":len(eligible),
                 "n_scheduled_partitions":sum((int(r["Norder"]),int(r["Npix"])) in selected for r in eligible)}
        print("HATS",kind,"source_partitions",len(items),"intersecting",original_count,
              "bounded",len(matching),flush=True)
        for pr in matching:
            order,pixel=int(pr["Norder"]),int(pr["Npix"])
            pixel_dir=parent+("dataset/Norder=%d/Dir=%s/Npix=%d/" % (order,pr["Dir"],pixel))
            listing_url=BASE+"?list-type=2&prefix="+quote(pixel_dir,safe="/")
            entry={"kind":kind,"Norder":order,"Npix":pixel,"Dir":pr["Dir"],
                   "advertised_rows":int(pr["num_rows"]),"prefix_listing_url":listing_url,
                   "read_timestamp_utc":frozen.utc()}
            try:
                listing=get(listing_url)
                entry["listing_sha256"]=sha(listing.content)
                xml=ET.fromstring(listing.content)
                keys=[x.find("{*}Key").text for x in xml.findall("{*}Contents")
                      if x.find("{*}Key") is not None and x.find("{*}Key").text.endswith(".parquet")]
                if len(keys)!=1:
                    raise ValueError("Expected exactly one Parquet object in HATS pixel: "+repr(keys[:8]))
                fullpath=keys[0]
                url=BASE+fullpath
                entry["object_uri"]=url
                head=requests.head(url,timeout=25)
                head.raise_for_status()
                entry["content_length"]=int(head.headers.get("Content-Length","0"))
                entry["s3_etag"]=head.headers.get("ETag")
                try:
                    if fs is None:raise RuntimeError("S3 filesystem unavailable")
                    stream=fs.open_input_file(BUCKET+"/"+fullpath)
                    entry["transport"]="s3-native"
                except Exception as e:
                    entry["native_transport_error"]=str(e)[:240]
                    stream=HTTPRangeFile(url,entry["content_length"],entry["s3_etag"])
                    entry["transport"]="s3-http-ranges"
                pf=pq.ParquetFile(stream)
                available=set(pf.schema_arrow.names)
                columns=[c for c in FIELDS if c in available]
                essential=set(FIELDS[:5]+FIELDS[14:22])
                if not essential.issubset(available):
                    raise ValueError("Missing essential HATS columns: "+str(sorted(essential-available)))
                entry["projected_columns"]=columns
                entry["projected_schema_sha256"]=sha(("\n".join(columns)).encode())
                entry["num_rows"]=pf.metadata.num_rows
                entry["num_row_groups"]=pf.num_row_groups
                if entry["num_rows"]!=entry["advertised_rows"]:
                    raise ValueError("Publisher row count mismatch")
                source_fields=[f for f in fields if pixel in healpix_pixels(order,f)]
                entry["field_ids"]=[f["id"] for f in source_fields]
                count_broad=0;count_probes=0
                # Read row group by row group to keep RAM bounded.
                for i in range(pf.num_row_groups):
                    tab=pf.read_row_group(i,columns=columns,use_threads=False)
                    rr=np.asarray(tab["ra"].to_numpy(),dtype=float)
                    dd=np.asarray(tab["dec"].to_numpy(),dtype=float)
                    finite=np.isfinite(rr)&np.isfinite(dd)
                    cone=np.zeros(len(rr),dtype=bool)
                    for f in source_fields:cone|=finite&inside(rr,dd,f)
                    # Source-level CSV signed fluxes remain distinguishable from NULL.
                    near=np.zeros(len(rr),dtype=bool)
                    for c in CONF["known"]:
                        near |= finite & (inside(rr,dd,{"ra":c["ra"],"dec":c["dec"],"r":5/3600}))
                    idx=np.flatnonzero(cone|near)
                    if not len(idx):continue
                    subset=tab.take(pa.array(idx)).to_pylist()
                    for j,x in zip(idx,subset):
                        x={k:scalar(v) for k,v in x.items()}
                        near_controls=[]
                        for c in CONF["known"]:
                            dist=frozen.sep_arcsec(rr[j],dd[j],c["ra"],c["dec"])
                            if dist<5:near_controls.append((c,dist))
                        if near_controls:
                            for c,dist in near_controls:
                                controls[c["name"]].append({"sep_arcsec":round(dist,4),
                                "object_id":str(x["object_id"]),"J_uJy":frozen.num(x.get("flux_j_3fwhm_aper")),
                                "pass_broad":broad(x),"pass_strict":not frozen.reasons({k:str(v) if v is not None else "" for k,v in x.items()})})
                            count_probes+=1
                        if not cone[j] or not broad(x):continue
                        match=[f["id"] for f in source_fields if inside(np.array([rr[j]]),np.array([dd[j]]),f)[0]]
                        x={k:str(v) if v is not None else "" for k,v in x.items()}
                        x["_table"]="euclid_dr1_mer_catalogue_"+kind
                        x["_field"]=",".join(match)
                        key=kind+"/"+x["object_id"]
                        found.setdefault(key,x)
                        count_broad+=1
                if entry["transport"]=="s3-http-ranges":
                    entry["range_sha256_receipts"]=stream.ranges
                stream.close()
                entry["broad_cone_rows_before_dedup"]=count_broad
                entry["positive_control_objects"]=count_probes
                entry["status"]="READ"
                print("PART",kind,order,pixel,"rows",entry["num_rows"],
                      "broad",count_broad,"controls",count_probes,flush=True)
            except Exception as e:
                entry["status"]="FAILED"
                entry["error"]=str(e)[:550]
                failures.append({"stage":"parquet","kind":kind,"Norder":order,
                                 "Npix":pixel,"error":entry["error"]})
                print("PART FAILED",kind,order,pixel,entry["error"],flush=True)
            records.append(entry)
    track=frozen.toy_tracks()
    strict=[];reject={}
    for r in found.values():
        bad=frozen.reasons(r)
        if bad:
            for a in set(bad):reject[a]=reject.get(a,0)+1
        else:strict.append(frozen.score(r,track))
    strict.sort(key=lambda x:(-x["J_SNR"],x["id"]))
    shortlist=strict[:CONF["shortlist"]]
    # Independence: no TAP-based crosschecks claimed when TAP service is down.
    for row in shortlist:
        row["spectra_status"]={v:"NOT_CHECKED_service_or_cone_coverage_unverified"
                for v in ["DESI_DR1","SDSS","ESO","JWST_MAST"]}
        row["radio_xray_status"]={v:"NOT_CHECKED"
                for v in ["eRASS1","VLASS","LoTSS","RACS"]}
        row["unWISE_W1_W2"]="NOT_CHECKED_TAP_UNAVAILABLE"
        row["Gaia_CatWISE_motion"]="NOT_CHECKED_TAP_UNAVAILABLE"
        row["calibrated_P_quasar_vs_dwarf"]=None
        row["classification"]="UNVETTED_COLOUR_PROXY"
    (OUT/"shortlist.json").write_text(json.dumps(shortlist,indent=2)+"\n")
    (OUT/"partitions.json").write_text(json.dumps(records,indent=2)+"\n")
    (OUT/"source_receipts.json").write_text(json.dumps(raws,indent=2)+"\n")
    # CSV of all broad colour rows (not raw HATS bytes); audit of all selected rows.
    columns=["object_id","ra","dec","_table","_field","vis_det","flux_y_3fwhm_aper",
       "fluxerr_y_3fwhm_aper","flux_j_3fwhm_aper","fluxerr_j_3fwhm_aper",
       "flux_h_3fwhm_aper","fluxerr_h_3fwhm_aper","ellipticity","spurious_prob",
       "flag_j","flag_h"]
    with (OUT/"broad_rows.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,columns);w.writeheader()
        for r in sorted(found.values(),key=lambda x:(x["_table"],x["object_id"])):
            w.writerow({k:r.get(k) for k in columns})
    controls_result=[]
    for c in CONF["known"]:
        obj=sorted(controls[c["name"]],key=lambda x:x["sep_arcsec"])
        controls_result.append({"name":c["name"],"literature_z":c["z"],
            "nearest":obj[:3],"catalog_recovered":bool(obj),
            "dropout_cut_recovered":any(x["pass_broad"] and x["pass_strict"] for x in obj)})
    count_good=sum(r["status"]=="READ" for r in records)
    complete_field_count=sum(v["status"]=="PARTITIONS_ATTEMPTED" for v in field_completeness.values())
    # A "fully scheduled" field still requires all those partitions to have status READ.
    for f in CONF["fields"]:
        details=field_completeness.get(f["id"])
        if not details or details["status"]!="PARTITIONS_ATTEMPTED":continue
        needed=[r for r in records if f["id"] in r.get("field_ids",[])]
        if not needed or any(r["status"]!="READ" for r in needed):
            details["status"]="FAILED_OR_UNREAD_PARTITION"
        else:details["status"]="PARTITIONS_READ_NOT_SKY_COMPLETENESS"
    summary={"time_utc":frozen.utc(),"release":CONF["release"],"release_date":CONF["release_date"],
       "wide_catalog_doi":"10.26131/IRSA673","deep_catalog_doi":"10.26131/IRSA675",
       "original_freeze_sha256":sha((ROOT/"freeze.json").read_bytes()),
       "s3_script_sha256":sha(Path(__file__).read_bytes()),
       "source_metadata":raws,"queried_cones":CONF["fields"],
       "nominal_requested_cones_deg2_NOT_VERIFIED_IMAGE_FOOTPRINT":2.014698,
       "field_partitions":field_completeness,
       "partitions_scheduled":len(records),"partitions_read":count_good,
       "partitions_failed":len(records)-count_good,
       "sources_satisfying_broad_NIR_VIStable_colours":len(found),
       "strict_colour_shape_proxy_count":len(strict),
       "strict_shortlist_count":len(shortlist),
       "strict_unvetted_beyond_shortlist":max(0,len(strict)-len(shortlist)),
       "rejection_flags_not_disjoint":reject,
       "known_quasars":controls_result,
       "discovery_level":False,"confirmed_new_quasars":0,
       "status":"S3_READ" if count_good and not failures else "S3_PARTIAL_OR_BLOCKED",
       "failure_details":failures,
       "limitations":["not full DR1 footprint; requested cones not measured effective area",
           "point-like classification NOT established by ellipticity",
           "toy colour likelihood not a calibrated probability",
           "no confirmed spectra, external WISE/radio/Xray or image/artefact controls",
           "full known z>6.5 quasar population not cross-checked"]}
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    hashes={str(p.relative_to(ROOT)):{"bytes":p.stat().st_size,"sha256":sha(p.read_bytes())}
           for p in sorted(OUT.rglob("*")) if p.is_file() and p.name!="hashes.json"}
    (OUT/"hashes.json").write_text(json.dumps(hashes,indent=2)+"\n")
    print("S3 FINAL",json.dumps({"status":summary["status"],"partitions":count_good,
         "broad":len(found),"strict":len(strict),"shortlist":len(shortlist)}),flush=True)
if __name__=="__main__":
    main()
