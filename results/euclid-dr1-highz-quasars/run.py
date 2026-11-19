#!/usr/bin/env python3
"""Frozen Euclid DR1-Foundation cone pilot. Requests only, plus Python standard library.
Every upstream response is hashed and retained. Absence is never a certified nondetection.
"""
import csv,datetime,hashlib,io,json,math,sys
from pathlib import Path
import requests
BASE=Path(__file__).resolve().parent
OUT=BASE/"output"
RAW=OUT/"raw"
CFG=json.loads((BASE/"freeze.json").read_text())
ENDPOINT="https://irsa.ipac.caltech.edu/TAP/sync"
HEAD={"User-Agent":"Open-Space-DR1-Foundation-Quasar-Pilot/1.0"}
COLS=("object_id ra dec tileid vis_det point_like_prob spurious_prob "
      "det_quality_flag ellipticity mumax_minus_mag deblended_flag blended_prob "
      "gaia_id flag_vis flag_y flag_j flag_h "
      "flux_vis_3fwhm_aper fluxerr_vis_3fwhm_aper "
      "flux_y_3fwhm_aper fluxerr_y_3fwhm_aper "
      "flux_j_3fwhm_aper fluxerr_j_3fwhm_aper "
      "flux_h_3fwhm_aper fluxerr_h_3fwhm_aper "
      "flux_g_ext_decam_3fwhm_aper fluxerr_g_ext_decam_3fwhm_aper "
      "flux_r_ext_decam_3fwhm_aper fluxerr_r_ext_decam_3fwhm_aper "
      "flux_i_ext_decam_3fwhm_aper fluxerr_i_ext_decam_3fwhm_aper "
      "flux_z_ext_decam_3fwhm_aper fluxerr_z_ext_decam_3fwhm_aper "
      "flux_i_ext_panstarrs_3fwhm_aper fluxerr_i_ext_panstarrs_3fwhm_aper "
      "flux_z_ext_panstarrs_3fwhm_aper fluxerr_z_ext_panstarrs_3fwhm_aper").split()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest(b):return hashlib.sha256(b).hexdigest()
def num(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None
def mag(x):return 23.9-2.5*math.log10(x) if x is not None and x>0 else None
def sep_arcsec(a,b,c,d):
    da=math.radians(a-c);db=math.radians(b-d)
    s=math.sin(db/2)**2+math.cos(math.radians(b))*math.cos(math.radians(d))*math.sin(da/2)**2
    return 7200*math.degrees(math.asin(min(1,math.sqrt(max(0,s)))))
def query_sql(field,controls=False):
    radius=field["r"]
    cols="object_id,ra,dec,vis_det,flux_j_3fwhm_aper" if controls else ",".join(COLS)
    clause=("1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',%.7f,%.7f,%.7f))" %
            (field["ra"],field["dec"],radius))
    if not controls:
        clause+=(" AND vis_det=0 AND flux_j_3fwhm_aper>0.65"
                 " AND flux_h_3fwhm_aper>0 AND fluxerr_j_3fwhm_aper>0"
                 " AND fluxerr_h_3fwhm_aper>0"
                 " AND flux_y_3fwhm_aper<0.85*flux_j_3fwhm_aper"
                 " AND flux_j_3fwhm_aper>0.65*flux_h_3fwhm_aper"
                 " AND (spurious_prob<0.7 OR spurious_prob IS NULL)")
    return "SELECT TOP %d %s FROM %s WHERE %s" % (CFG["cap"],cols,field["table"],clause)
def fetch(label,endpoint,sql,receipts,timeout=90):
    t={"label":label,"endpoint":endpoint,"utc":utc(),"adql":sql,"adql_sha256":digest(sql.encode())}
    try:
        rr=requests.post(endpoint,data={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv",
                          "QUERY":sql},headers=HEAD,timeout=timeout)
        t.update({"http_code":rr.status_code,"response_bytes":len(rr.content),
                  "response_sha256":digest(rr.content)})
        path=RAW/(label+".response")
        path.write_bytes(rr.content)
        t["response_file"]=str(path.relative_to(BASE))
        rr.raise_for_status()
        txt=rr.content.decode("utf-8-sig","replace")
        if txt.lstrip().startswith("<"):
            raise ValueError("XML/error rather than CSV: "+txt[:170])
        reader=csv.DictReader(io.StringIO(txt))
        if not reader.fieldnames or "object_id" not in [x.lower().strip() for x in reader.fieldnames]:
            raise ValueError("Missing object_id header: "+txt[:170])
        rows=[{k.lower().strip():v for k,v in r.items()} for r in reader]
        t["row_count"]=len(rows)
        return rows
    except Exception as e:
        t["error"]=str(e)[:250]
        return None
    finally:
        t["done_utc"]=utc();receipts.append(t)
def band(r,b):
    f=num(r.get("flux_"+b+"_3fwhm_aper"))
    e=num(r.get("fluxerr_"+b+"_3fwhm_aper"))
    return f,e,(f/e if f is not None and e and e>0 else None)
def metrics(r):
    f={b:band(r,b) for b in ("vis","y","j","h")}
    fm={b:mag(f[b][0]) for b in f}
    yj=fm["y"]-fm["j"] if fm["y"] is not None and fm["j"] is not None else None
    jh=fm["j"]-fm["h"] if fm["j"] is not None and fm["h"] is not None else None
    y_limit=mag(max(2*f["y"][1],f["y"][0] or -1))-fm["j"] if f["y"][1] and fm["j"] is not None else None
    return f,yj,jh,y_limit
def reasons(r):
    f,yj,jh,yl=metrics(r);bad=[]
    if f["j"][2] is None or f["j"][2]<6:bad.append("J_SNR")
    if f["h"][2] is None or f["h"][2]<4:bad.append("H_SNR")
    if f["vis"][2] is None or abs(f["vis"][2])>=2:bad.append("VIS_MISSING_OR_DETECTION")
    if (yj is None and (yl is None or yl<.40)) or (yj is not None and yj<.40):bad.append("Y_J")
    if jh is None or jh>.45:bad.append("J_H")
    if num(r.get("ellipticity")) is None or num(r["ellipticity"])>=.3:bad.append("NIR_SHAPE_PROXY")
    for b in ("j","h"):
        if r.get("flag_"+b) not in ("0","0.0"):bad.append(b+"_FLAG")
    return bad
def toy_tracks():
    """Approximate QSO IGM-break top-hat tracks; NOT real Euclid throughput."""
    ans=[]
    for i in range(21):
        z=7+i*.1
        flux={}
        for b,(lo,hi) in {"y":(.97,1.19),"j":(1.19,1.54),"h":(1.54,1.99)}.items():
            vals=[lo+(k+.5)*(hi-lo)/160 for k in range(160)]
            flux[b]=sum(v**.5 for v in vals if v>.1216*(1+z))/160
        ans.append({"z":round(z,2),"yj":-2.5*math.log10(flux["y"]/flux["j"]) if flux["y"]>0 else 5,
                    "jh":-2.5*math.log10(flux["j"]/flux["h"])})
    return ans
def toy_p(yj,jh,track):
    if yj is None or jh is None:return None
    def gaussian(x,y,m,n,s,t):
        return math.exp(-.5*((x-m)/s)**2-.5*((y-n)/t)**2)/(2*math.pi*s*t)
    q=sum(gaussian(yj,jh,t["yj"],t["jh"],.4,.35) for t in track)/len(track)
    d=sum(gaussian(yj,jh,m,n,.43,.42) for m,n in ((.75,.55),(.15,-.15),(1.,-.35)))/3
    return round(q/(q+1000*d),8) if q+1000*d else None
def score(r,track):
    f,yj,jh,yl=metrics(r)
    known=[]
    for c in CFG["known"]:
        if sep_arcsec(float(r["ra"]),float(r["dec"]),c["ra"],c["dec"])<3:
            known.append(c["name"])
    return {
        "id":r["object_id"],"ra":num(r["ra"]),"dec":num(r["dec"]),"field":r["_field"],
        "table":r["_table"],"Y_minus_J_AB":yj,"Y_minus_J_2sig_limit":yl,"J_minus_H_AB":jh,
        "VIS_SNR":f["vis"][2],"J_SNR":f["j"][2],"H_SNR":f["h"][2],
        "ellipticity":num(r.get("ellipticity")),"point_like_prob":num(r.get("point_like_prob")),
        "toy_P_Q_given_Q_or_dwarf_UNCALIBRATED":toy_p(yj,jh,track),
        "known_published":known,"status":"KNOWN" if known else "UNCONFIRMED_COLOUR_SELECTION",
        "nir_point_source_validated":False,
        "ground_optical":{b:{"flux_ujy":num(r.get("flux_"+b+"_3fwhm_aper")),
                    "err_ujy":num(r.get("fluxerr_"+b+"_3fwhm_aper"))}
                    for b in ("g_ext_decam","r_ext_decam","i_ext_decam","z_ext_decam",
                              "i_ext_panstarrs","z_ext_panstarrs")}
    }
def external(item,receipts):
    ra,dec=item["ra"],item["dec"]
    resp={}
    for table,cols in (
       ("unwise_2019","ra,dec,flux_1,dflux_1,flux_2,dflux_2,fracflux_1,fracflux_2,flags_unwise_1,flags_unwise_2"),
       ("catwise_2020","ra,dec,w1mpro,w2mpro,pmra,pmdec,sigpmra,sigpmdec")):
        sql=("SELECT TOP 10 "+cols+" FROM "+table+
             " WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',%.7f,%.7f,%.7f))" %
             (ra,dec,3/3600))
        rows=fetch("x_"+item["id"]+"_"+table,ENDPOINT,sql,receipts,timeout=28)
        if rows is None:resp[table]={"status":"QUERY_ERROR"};continue
        nearby=sorted((sep_arcsec(ra,dec,float(x["ra"]),float(x["dec"])),x)
                      for x in rows if num(x.get("ra")) is not None and num(x.get("dec")) is not None)
        resp[table]={"status":"CATALOG_MATCH" if nearby else "NO_CATALOG_ROW_NOT_A_FLUX_LIMIT",
                     "nearest_arcsec":nearby[0][0] if nearby else None,
                     "row":nearby[0][1] if nearby else None,
                     "caution":"unWISE blended beam; nearest catalog detection not forced target photometry"}
    q=("SELECT TOP 10 source_id,ra,dec,pmra,pmdec,parallax,parallax_over_error FROM gaiadr3.gaia_source"
       " WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',%.7f,%.7f,%.7f))" %
       (ra,dec,2/3600))
    g=fetch("x_"+item["id"]+"_gaia","https://gea.esac.esa.int/tap-server/tap/sync",
            q,receipts,timeout=28)
    resp["Gaia_DR3"]={"status":"QUERY_ERROR" if g is None else "CATALOG_RETURNED",
                       "rows":g if g is not None else None}
    resp["spectra"]={k:"PENDING" for k in ("DESI_DR1","SDSS","ESO","JWST_MAST")}
    resp["high_energy_radio"]={k:"PENDING" for k in ("eRASS1","VLASS","LoTSS","RACS")}
    resp["independent_optical"]={k:"PENDING" for k in ("Legacy_Surveys","DES","PanSTARRS")}
    return resp
def run():
    RAW.mkdir(parents=True,exist_ok=True)
    receipts=[];fields=[];allrows={}
    for f in CFG["fields"]:
        rows=fetch(f["id"],ENDPOINT,query_sql(f),receipts,timeout=95)
        fields.append({"id":f["id"],"table":f["table"],"ra":f["ra"],"dec":f["dec"],
                       "radius_deg":f["r"],"requested_area_deg2":round(2*math.pi*(180/math.pi)**2*
                              (1-math.cos(math.radians(f["r"]))),6),
                       "status":"FAILED" if rows is None else "RETURNED",
                       "rows":len(rows) if rows is not None else None,
                       "hit_cap":rows is not None and len(rows)>=CFG["cap"]})
        if rows is not None:
            for x in rows:
                x["_table"]=f["table"];x["_field"]=f["id"]
                allrows.setdefault(f["table"]+"/"+x["object_id"],x)
        print("Field:",fields[-1],flush=True)
    track=toy_tracks();strict=[];reject={}
    for r in allrows.values():
        bad=reasons(r)
        if bad:
            for b in set(bad):reject[b]=reject.get(b,0)+1
        else:strict.append(score(r,track))
    strict.sort(key=lambda x:(-x["J_SNR"],x["id"]))
    shortlist=strict[:CFG["shortlist"]]
    for c in shortlist:c["crosschecks"]=external(c,receipts)
    controls=[]
    for c in CFG["known"]:
        trials=[]
        for table in ("euclid_dr1_mer_catalogue_wide_survey","euclid_dr1_mer_catalogue_deep_survey"):
            f={"ra":c["ra"],"dec":c["dec"],"r":5/3600,"table":table}
            label=("ctrl_"+c["name"].replace("+","p").replace(".","_")+"_"+table.split("_")[-2])
            rows=fetch(label,ENDPOINT,query_sql(f,controls=True),receipts,timeout=30)
            trials.append({"table":table,"status":"ERROR" if rows is None else "FOUND" if rows else "NO_MATCH",
                           "n":len(rows) if rows is not None else None})
        controls.append({"name":c["name"],"z_literature":c["z"],"trials":trials,
                         "note":"Catalog lookup not a positive recovery of frozen dropout cut"})
    rawfields=["id","ra","dec","field","status","Y_minus_J_AB","J_minus_H_AB",
               "VIS_SNR","J_SNR","H_SNR","ellipticity","toy_P_Q_given_Q_or_dwarf_UNCALIBRATED"]
    with (OUT/"shortlist.csv").open("w",newline="") as h:
        w=csv.DictWriter(h,rawfields);w.writeheader()
        for item in shortlist:w.writerow({k:item.get(k) for k in rawfields})
    (OUT/"shortlist.json").write_text(json.dumps(shortlist,indent=2)+"\n")
    (OUT/"tracks_toy.json").write_text(json.dumps({"NOT_CALIBRATED":True,"qso":track,
       "toy_dwarf_LTY_YJ_JH":[[.75,.55],[.15,-.15],[1,-.35]],"prior_Q_to_dwarf":1/1000},indent=2)+"\n")
    summary={
      "time_utc":utc(),"release":CFG["release"],"release_date":CFG["release_date"],
      "release_footprint_deg2_approx":CFG["release_footprint_sqdeg_approx"],
      "freeze_sha256":digest((BASE/"freeze.json").read_bytes()),
      "script_sha256":digest(Path(__file__).read_bytes()),
      "fields":fields,"sum_requested_cone_area_deg2_NOT_VALID_EFFECTIVE_AREA":sum(x["requested_area_deg2"] for x in fields),
      "n_queries_failed":sum(x["status"]=="FAILED" for x in fields),
      "n_queries_capped":sum(x["hit_cap"] for x in fields),
      "unique_broad_rows":len(allrows),"strict_proxy_rows":len(strict),
      "shortlist_count":len(shortlist),"strict_rows_not_independently_vetted":max(0,len(strict)-len(shortlist)),
      "reject_counts_overlapping":reject,"known_controls":controls,
      "status":"PARTIAL" if any(x["status"]=="FAILED" or x["hit_cap"] for x in fields) else "BOUNDED_CONES_COMPLETED",
      "new_quasar_confirmed":False,"new_spectrum_confirmed":False,"discovery_claim":False,
      "probability_status":"Toy track versus toy L/T/Y. Not empirical or astrophysical posterior.",
      "pending":["independent Legacy/DES/PanSTARRS flux limits","NIR PSF shape and image artifacts",
      "dwarf spectral-library bandpass calibration","red galaxy and injection controls",
      "Gaia/CatWISE significant motion tests","DESI DR1 SDSS ESO JWST spectra",
      "eRASS1 VLASS LoTSS RACS footprint-aware checks","measured effective footprint"]
    }
    (OUT/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    (OUT/"requests.json").write_text(json.dumps(receipts,indent=2)+"\n")
    hashes={str(p.relative_to(BASE)):{"sha256":digest(p.read_bytes()),"bytes":p.stat().st_size}
            for p in sorted(OUT.rglob("*")) if p.is_file()}
    (OUT/"hashes.json").write_text(json.dumps(hashes,indent=2)+"\n")
    print(json.dumps({"status":summary["status"],"rows":len(allrows),
                      "strict":len(strict),"shortlist":len(shortlist)}),flush=True)
if __name__=="__main__":run()
