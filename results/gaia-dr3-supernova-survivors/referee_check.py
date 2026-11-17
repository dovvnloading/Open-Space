"""Independent deposited-evidence checks for PR #145; no network is needed.

Run after (or before) the authoring scripts. The separate workflow actually reruns
the source calculations and checks that they regenerate exactly the same files.
"""
import base64
import csv
import hashlib
import json
import math
import struct
from collections import Counter
from pathlib import Path
from statistics import NormalDist
from xml.etree import ElementTree

ROOT = Path(__file__).resolve().parent
K = 4.74047
Z95 = NormalDist().inv_cdf(.95)
FLOOR = .020
checks = []

def table(path):
    with (ROOT/path).open(newline="") as file:
        return list(csv.DictReader(file))

def getjson(path):
    return json.loads((ROOT/path).read_text())

def eq(value, expected, label):
    assert value == expected, (label, value, expected)
    checks.append(f"{label}: {value}")

def near(value, expected, tol, label):
    assert abs(value-expected) <= tol, (label, value, expected)
    checks.append(f"{label}: {value:.6f}")

def clean_hashes():
    manifest=ROOT/"SHA256SUMS"
    rows=[line.split("  ", 1) for line in manifest.read_text().splitlines() if line.strip()]
    for expected, rel in rows:
        path=ROOT/rel
        assert path.is_file(), ("missing", rel)
        actual=hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual==expected, ("SHA256 mismatch",rel,actual,expected)
    checks.append(f"SHA256SUMS: {len(rows)} checked-in files verified")
    n=0
    for section in ("search","control"):
        for receipt in (ROOT/section).glob("*_receipt.json"):
            record=json.loads(receipt.read_text())
            fname=receipt.name.replace("_receipt.json", ".csv")
            if not isinstance(record,dict) or "sha256" not in record:
                continue
            raw=receipt.parent/fname
            if not raw.exists():
                continue
            assert hashlib.sha256(raw.read_bytes()).hexdigest()==record["sha256"], fname
            if "nbytes" in record:
                assert len(raw.read_bytes())==record["nbytes"], fname
            if "bytes" in record:
                assert len(raw.read_bytes())==record["bytes"], fname
            if "rows" in record:
                assert len(table(f"{section}/{fname}"))==record["rows"], fname
            n+=1
    validation=getjson("search/parent_validation.json")
    assert hashlib.sha256((ROOT/"search/parent_validation.xml").read_bytes()).hexdigest()==validation["sha256"]
    eq(validation["tap_statuses"],["OK"],"Gaia result status")
    eq(validation["distinct_ids"],12,"Gaia independent ID count")
    origin=getjson("origin/gaia_receipt.json")
    for record in origin:
        p=ROOT/"origin"/record["file"]
        assert hashlib.sha256(p.read_bytes()).hexdigest()==record["sha256"], p
        assert len(list(csv.DictReader(p.open())))==record["rows"]
    checks.append(f"query receipts: {n} deposited CSV replies + VOTable + {len(origin)} orbit batches verified")

def union_check():
    overlap=table("union/overlap.csv")
    eq(len(overlap),15,"distinct nonempty survey masks")
    eq({int(r["mask"]) for r in overlap},set(range(1,16)),"all 15 survey combinations")
    for field, expected in (("unique",113166),("eligible_snr3",75684),("eligible_snr5",67433)):
        eq(sum(int(r[field]) for r in overlap),expected,"overlap "+field)
    summary=getjson("union/summary.json")
    eq(summary["unique"],113166,"census distinct DR3 IDs")
    for i,(label,bit,expected_ids,expected5) in enumerate([
        ("DESI",1,63152,33322),("SDSS_DR19",2,10364,7843),
        ("SDSS_DR20",4,50862,36400),("LAMOST",8,13850,8956)
    ]):
        eq(sum(int(r["unique"]) for r in overlap if int(r["mask"])&bit),expected_ids,label+" IDs")
        eq(sum(int(r["eligible_snr5"]) for r in overlap if int(r["mask"])&bit),expected5,label+" eligible")
    for cut,n in [(3,75684),(5,67433)]:
        reported=next(r for r in summary["results"] if r["sample"]=="UNION" and r["corrected_snr_min"]==cut)
        eq(reported["eligible"],n,f"summary at S/N={cut}")
        eq(reported["above1000"],0,f"above 1000 at S/N={cut}")
    top=table("union/top10.csv")
    for cut in (3,5):
        rows=[r for r in top if int(r["threshold"])==cut]
        eq(len(rows),10,f"top-ten list at {cut}")
        values=[]
        for r in rows:
            value=K*math.hypot(float(r["pmra"]),float(r["pmdec"]))/(float(r["parallax"])-float(r["zpt_mas"]))
            near(value,float(r["vt_kms"]),.00001,f"speed recompute {r['source_id'][:5]} cut {cut}")
            values.append(value)
        near(max(values),704.5551144887545,.00001,f"max vt at {cut}")
        eq(max(rows,key=lambda r:float(r["vt_kms"]))["source_id"],"3688712561723372672","max Gaia ID")
    eq(round(1e6*(-math.expm1(math.log(.05)/67433)),1),44.4,"conditional binomial upper bound per million")

def score(r,zp,v):
    pma,pmd=[float(r[k]) for k in ("pmra","pmdec")]
    mu=math.hypot(pma,pmd)
    pa=float(r["parallax"])-zp
    sp,sma,smd=[float(r[k+"_error"]) for k in ("parallax","pmra","pmdec")]
    ra,rb,rc=[float(r[k]) for k in ("parallax_pmra_corr","parallax_pmdec_corr","pmra_pmdec_corr")]
    a,b=K*pma/(v*mu),K*pmd/(v*mu)
    sigma2=(sp*sp+a*a*sma*sma+b*b*smd*smd
            -2*a*ra*sp*sma-2*b*rb*sp*smd+2*a*b*rc*sma*smd+FLOOR**2)
    return (K*mu/v-pa)/math.sqrt(sigma2)

def gaia_check():
    eq(int(table("search/eligible_count.csv")[0]["n"]),1023709,"all-sky pre-speed count")
    eq(int(table("search/count.csv")[0]["n"]),12,"all-sky raw threshold count")
    parent=table("search/parent.csv")
    eq(len(parent),12,"Gaia enumerated candidates")
    eq(len({r["source_id"] for r in parent}),12,"Gaia unique source IDs")
    g={r["source_id"]:r for r in table("search/gaia.csv")}
    reported={r["source_id"]:r for r in table("search/kinematics.csv")}
    eq(set(g),set(reported),"candidate Gaia/kinematics ID join")
    controls={r["source_id"]:r["name"] for r in table("control/controls.csv")}
    for r in parent:
        mu=float(r["pm"])
        par=float(r["parallax"])
        eq(bool(mu>35 and float(r["bp_rp"])<1.2 and par/float(r["parallax_error"])>3
            and K*mu/par>1000 and float(r["phot_g_mean_mag"])+5*math.log10(par)-10>3.5),
            True,"all-sky raw cuts "+r["source_id"])
    recovered=Counter()
    lead=None
    for sid,r in g.items():
        k=reported[sid]
        zero=float(k["zp_mas"])
        value=K*math.hypot(float(r["pmra"]),float(r["pmdec"]))/(float(r["parallax"])-zero)
        near(value,float(k["corrected_vt"]),.00001,"corrected vt "+sid)
        low,hi=10.,10000.
        for _ in range(70):
            mid=(low+hi)/2
            if score(r,zero,mid)>Z95: low=mid
            else: hi=mid
        near((low+hi)/2,float(k["vt_lower95"]),.0001,"95% lower vt "+sid)
        if float(k["vt_lower95"])>1000: recovered[controls.get(sid,"OTHER")]+=1
        if sid=="6094201393481158656": lead=k
    eq(dict(recovered),{"D6-1":1,"D6-3":1},"Gaia 95% survivors")
    near(float(lead["corrected_vt"]),1117.447756,.001,"unresolved target point speed")
    near(float(lead["vt_lower95"]),810.750319,.001,"unresolved target 95% lower")
    near(float(lead["sep_arcsec"]),1.349198,.001,"unresolved target neighbour")
    eq(lead["neighbour_3arcsec"],"6094201393488527744","neighbour source ID")
    eq(lead["gspphot_teff"]!="",True,"cool-star Gaia atmosphere fit present")
    others=[float(r["vt_lower95"]) for sid,r in reported.items() if sid not in controls]
    near(max(others),810.750319,.001,"max non-control 95% lower")

def controls_check():
    control=table("control/controls.csv")
    gaia={r["source_id"]:r for r in table("control/gaia.csv")}
    rec={r["source_id"]:r for r in table("control/recovery.csv")}
    tally=Counter()
    failures=Counter()
    for star in control:
        g=gaia[star["source_id"]]
        p,e=float(g["parallax"]),float(g["parallax_error"])
        mu=math.hypot(float(g["pmra"]),float(g["pmdec"]))
        first=mu>35 and float(g["bp_rp"])<1.2
        sn=first and p/e>3
        fast=sn and K*mu>1000*p and (float(g["phot_g_mean_mag"])+5*math.log10(p)-10>3.5)
        eq(str(fast),rec[star["source_id"]]["pass_fast"],"control selection "+star["name"])
        if fast: tally[star["class"]]+=1
        if star["class"]=="D6" and not fast:
            failures["D6_failed_snr" if first and not sn else "D6_other"]+=1
    eq((tally["D6"],tally["LP40"],tally["He"]),(3,0,0),"control recoveries D6 LP40 He")
    eq(failures["D6_failed_snr"],8,"D6 missed at parallax S/N >3 gate")
    eq(len(control),20,"published controls total")
    summary=getjson("control/recovery_summary.json")
    for clas,n in [("D6",3),("LP40",0),("He",0)]:
        eq(summary["control_recovery"][clas]["pass_fast_quality"],n,clas+" quality recovery")
    eq(sum(1 for r in control if r["class"]=="D6" and float(gaia[r["source_id"]]["phot_g_mean_mag"])>18),8,"faint D6 controls (G>18)")

def injection_check():
    rows=table("search/injection.csv")
    eq(len(rows),64,"injection bins")
    count={int(float(r["g_lo"])):int(r["templates"]) for r in rows}
    eq(sum(count.values()),744,"injection Gaia templates")
    for r in rows:
        eq(int(r["trials"]),int(r["templates"])*256,"injection trials "+r["g_lo"]+" "+r["distance_kpc"]+" "+r["true_3d_kms"])
    for d,value in [(1,.562302),(2,.362970),(4,.018880)]:
        subset=[r for r in rows if float(r["g_lo"])==17 and float(r["g_hi"])==18
                and int(r["distance_kpc"])==d and int(r["true_3d_kms"])==1500]
        eq(len(subset),1,"unique injection stratum")
        near(float(subset[0]["strong_fraction"]),value,.000001,f"injection d={d}kpc")

def orbit_check():
    rows=table("origin/remnant_checks.csv")
    eq(len({r["source_id"] for r in rows}),7,"remnant orbit controls")
    assoc=[r for r in rows if r["name"]=="D6-2" and r["remnant"]=="G70.0-21.5"]
    eq(len(assoc),1,"D6-2 remnant record")
    rec=assoc[0]
    eq(rec["hit"],"True","D6-2 shell overlap")
    near(float(rec["time_yr"]),93750,0,"remnant lookback years")
    near(float(rec["separation_deg"]),.90469,.00001,"remnant angular separation")
    assert float(rec["separation_deg"]) < float(rec["max_radius_deg"])
    states=table("origin/states.csv")
    for name in ("D6-1","D6-2","D6-3"):
        g=next(r for r in states if r["name"]==name and r["distance_model"]=="bj_geo")
        eq(int(g["draws"]),4096,name+" posterior draws")
        eq((float(g["p_escape_mc17"]),float(g["p_escape_gala"])),(1.0,1.0),name+" unbound draws")
    for name,expected in [("D6-1",1991.22),("D6-2",1125.553),("D6-3",2303.288)]:
        g=next(r for r in states if r["name"]==name and r["distance_model"]=="bj_geo")
        near(float(g["vgal_kms_50"]),expected,.001,name+" galactocentric speed")

if __name__ == "__main__":
    clean_hashes()
    union_check()
    gaia_check()
    controls_check()
    injection_check()
    orbit_check()
    print(json.dumps({"status":"PASS","tested":len(checks),"headlines":[s for s in checks if any(k in s for k in (
          "overlap ","all-sky pre-speed","Gaia enumerated","control recoveries","D6 missed",
          "injection d=","remnant angular","max vt at"))]},indent=2))
