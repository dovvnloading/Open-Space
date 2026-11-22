"""Runner-only SIMBAD identity and bibliography check for two explicitly named XRBs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import csv
import gzip
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import re
from urllib.parse import quote, urlencode
import xml.etree.ElementTree as ET

from astropy.coordinates import SkyCoord
from astropy import units as u

from literature import Archive, SIMBADS, dump, literal, norm, sky_distance, uniq

TARGETS = {
    "1E1145.1-6141":["1E 1145.1-6141","V* V830 Cen","V830 Cen","1E 1145.1-6141.0","SWIFT J1147.2-6156B"],
    "H1820-303":["4U 1820-30","4U 1820-303","H 1820-303","1A 1820-303","NGC 6624 X-1","X Sgr X-4","Sgr X-4","XB 1820-303"],
}
MIRRORS = ("https://simbad.cds.unistra.fr/simbad/","https://simbad.cfa.harvard.edu/simbad/")
BINARY_TYPES = {"hxb","lxb","hmxb","lmxb","xrb","highmassxraybinary","lowmassxraybinary"}


def decoded_rows(raw):
    reader=csv.DictReader(io.StringIO(raw.decode("utf-8-sig")))
    if not reader.fieldnames or any("<" in v for v in reader.fieldnames):
        raise ValueError("Response was not TAP CSV")
    return [{k.lower():v for k,v in row.items()} for row in reader]


def perform(archive, checks, kind, url, params, parser):
    try:
        raw,meta=archive.fetch(url,params)
        out=parser(raw)
        checks.append({"method":kind,"url":meta["url"],"sha256":meta["sha256"],
                       "returned_rows":len(out)})
        return out
    except (ValueError,KeyError,TypeError,UnicodeError,ET.ParseError,csv.Error) as exc:
        checks.append({"method":kind,"url":url+("?"+urlencode(params) if params else ""),
                       "error":str(exc)[:450]})
        return None


def tap(archive, checks, kind, query, limit=128, nonempty=True):
    params={"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","MAXREC":str(limit),"QUERY":query}
    for url in SIMBADS:
        rows=perform(archive,checks,kind,url,params,decoded_rows)
        if rows is not None and (rows or not nonempty):
            return rows
    return []


def choose(rows, candidate, variants):
    names={norm(x) for x in variants}
    valid={}
    for row in rows:
        try:
            offset=sky_distance(candidate["ra"], candidate["dec"], float(row["ra"]), float(row["dec"]))
        except (ValueError,TypeError,KeyError):
            continue
        typed=re.sub("[^a-z]","",row.get("otype","").lower()) in BINARY_TYPES
        exact=bool(names & {norm(row.get("main_id","")),norm(row.get("query_id",""))})
        cluster=norm(row.get("main_id",""))=="ngc6624" or row.get("otype","").lower() in {"glc","cl*"}
        if offset<=1.0 and not cluster and (exact or typed):
            key=row.get("oid") or norm(row["main_id"])
            valid[key]={**row,"offset_arcmin":float(offset),"exact_target_alias":exact}
    exact=[v for v in valid.values() if v["exact_target_alias"]]
    matches=exact or list(valid.values())
    return matches[0] if len(matches)==1 else None


def sesame_rows(raw):
    root=ET.fromstring(raw)
    records=[]
    for r in root.iter():
        if r.tag.split("}")[-1]!="Resolver" or "SIMBAD" not in r.attrib.get("name","").upper():
            continue
        fields={x.tag.split("}")[-1]:(x.text or "").strip() for x in r}
        if all(fields.get(k) for k in ("oname","jradeg","jdedeg")):
            records.append({"main_id":fields["oname"],"ra":fields["jradeg"],"dec":fields["jdedeg"],
                            "otype":fields.get("otype","")})
    return records


def ascii_rows(raw):
    text=raw.decode(errors="replace")
    name=re.search(r"(?m)^\s*Object\s+(.+?)\s*---\s*(.+?)\s*---",text)
    coord=re.search(r"(?mi)^\s*Coordinates\s*\(ICRS[^:]*:\s*(\d{1,2}\s+\d{1,2}\s+[\d.]+)\s+([+-]\d{1,2}\s+\d{1,2}\s+[\d.]+)",text)
    if not name or not coord:
        return []
    pos=SkyCoord(coord.group(1),coord.group(2),unit=(u.hourangle,u.deg),frame="icrs")
    return [{"main_id":name.group(1).strip(),"otype":name.group(2).strip(),
             "ra":float(pos.ra.deg),"dec":float(pos.dec.deg)}]


def html_rows(raw):
    from bs4 import BeautifulSoup
    text=BeautifulSoup(raw,"html.parser")
    # HTML layouts vary. Accept a record only if it exposes a machine-readable
    # coordinate and a canonical object label; a query echo alone proves nothing.
    nodes=text.select("meta[property='og:title'], title")
    names=[x.get("content") or x.get_text(" ",strip=True) for x in nodes]
    coords=None
    for tr in text.find_all("tr"):
        row=" ".join(tr.stripped_strings)
        if not re.search(r"ICRS",row,re.I):
            continue
        m=re.search(r"(\d{1,2}\s+\d{1,2}\s+[\d.]+)\s+([+-]\d{1,2}\s+\d{1,2}\s+[\d.]+)",row)
        if m:
            coords=SkyCoord(m.group(1),m.group(2),unit=(u.hourangle,u.deg))
            break
    if coords is None or not names:
        return []
    label=re.sub(r"(?i)\s*[-|:]\s*SIMBAD.*$","",names[0]).strip()
    return [{"main_id":label,"ra":float(coords.ra.deg),"dec":float(coords.dec.deg),"otype":""}]


def lookup(target, input_path, directory):
    original=json.loads(input_path.read_text())
    candidate=original["candidate"]
    old=original.get("simbad") or {}
    receipt=old.get("query") or {}
    digest=receipt.get("sha256","")
    archived=[]
    if re.fullmatch("[0-9a-f]{64}",digest):
        blob=input_path.parent/"responses"/(digest+".gz")
        if blob.exists():
            raw=gzip.decompress(blob.read_bytes())
            if hashlib.sha256(raw).hexdigest()!=digest:
                raise ValueError("Mismatch in original SIMBAD source-response checksum")
            archived=decoded_rows(raw)
    archived_context=[{**r,"offset_arcmin":sky_distance(candidate["ra"],candidate["dec"],float(r["ra"]),float(r["dec"]))}
                      for r in archived if r.get("ra") and r.get("dec")]
    aliases=uniq(TARGETS[target]+[x for x in candidate.get("aliases",[]) if norm(x)!="ngc6624"])
    stamp=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H")
    token=f"{literal(stamp)}={literal(stamp)}"
    archive=Archive(directory/target)
    checks=[]
    base="SELECT DISTINCT b.oid,b.main_id,b.ra,b.dec,b.otype"
    exact=(base+",i.id AS query_id FROM basic AS b JOIN ident AS i ON i.oidref=b.oid WHERE i.id IN ("+
           ",".join(map(literal,aliases))+") AND "+token)
    rows=tap(archive,checks,"identifier TAP",exact)
    obj=choose(rows,candidate,aliases)
    method="identifier TAP"
    if obj is None:
        cone=(base+" FROM basic AS b WHERE 1=CONTAINS(POINT('ICRS',b.ra,b.dec),"+
              f"CIRCLE('ICRS',{candidate['ra']:.8f},{candidate['dec']:.8f},0.016666667)) AND "+token)
        rows=tap(archive,checks,"coordinate TAP",cone,1000)
        obj=choose(rows,candidate,aliases)
        method="coordinate TAP"
    if obj is None:
        obj=choose(archived,candidate,aliases)
        if obj:
            method="archived original SIMBAD query, resolved by target identifiers and coordinates"
    if obj is None:
        for alias in TARGETS[target][:3]:
            url="https://cds.unistra.fr/cgi-bin/nph-sesame/-oxp/S?"+quote(alias,safe="")
            rows=perform(archive,checks,"CDS Sesame SIMBAD resolver",url,None,sesame_rows) or []
            for row in rows:
                row["query_id"]=alias
            obj=choose(rows,candidate,aliases)
            if obj:
                method="CDS Sesame SIMBAD resolver"
                break
    if obj is None:
        for alias in TARGETS[target][:2]:
            for host in MIRRORS:
                for fmt,parse in (("ASCII",ascii_rows),("HTML",html_rows)):
                    rows=perform(archive,checks,"SIMBAD sim-id "+fmt,host+"sim-id",
                                 {"Ident":alias,"output.format":fmt},parse) or []
                    for row in rows:
                        row["query_id"]=alias
                    obj=choose(rows,candidate,aliases)
                    if obj:
                        method="SIMBAD sim-id "+fmt
                        break
                if obj: break
            if obj: break
    link=MIRRORS[0]+"sim-id?"+urlencode({"Ident":(obj or {}).get("main_id") or aliases[0]})
    bib=None
    xids=[]
    if obj:
        if obj.get("oid"):
            xids=tap(archive,checks,"canonical object cross-identifiers","SELECT id FROM ident WHERE oidref="+
                     str(int(obj["oid"]))+" AND "+token,3000)

        if obj.get("oid"):
            scope=" FROM has_ref AS h WHERE h.oidref="+str(int(obj["oid"]))+" AND "+token
            join=" FROM has_ref AS h JOIN ref AS r ON r.oidbib=h.oidbibref WHERE h.oidref="+str(int(obj["oid"]))+" AND "+token
        else:
            scope=" FROM has_ref AS h JOIN basic AS b ON h.oidref=b.oid WHERE b.main_id="+literal(obj["main_id"])+" AND "+token
            join=(" FROM has_ref AS h JOIN ref AS r ON r.oidbib=h.oidbibref JOIN basic AS b ON h.oidref=b.oid "+
                  "WHERE b.main_id="+literal(obj["main_id"])+" AND "+token)
        count=tap(archive,checks,"attached references count","SELECT COUNT(*) AS nref"+scope,1)
        count_kind="attached reference rows"
        if not count and obj.get("oid"):
            count=tap(archive,checks,"SIMBAD basic bibliography count",
                      f"SELECT nbref AS nref FROM basic WHERE oid={int(obj['oid'])} AND "+token,1)
            count_kind="SIMBAD basic.nbref"
        if not count:
            all_codes=tap(archive,checks,"direct bibliographic-code inventory","SELECT DISTINCT r.bibcode"+join,30000)
            if all_codes and len(all_codes)<30000:
                count=[{"nref":len(all_codes)}]
                count_kind="retrieved distinct attached bibcodes"
        refs=tap(archive,checks,"recent attached bibliography","SELECT TOP 20 r.bibcode,r.title"+join+" ORDER BY r.bibcode DESC",20)
        ref_url=SIMBADS[0]+"?"+urlencode({"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv",
                                          "QUERY":"SELECT r.bibcode,r.title"+join+" ORDER BY r.bibcode DESC"})
        bib={"reported_count":int(count[0]["nref"]) if count and str(count[0].get("nref","")).isdigit() else None,
             "count_basis":count_kind,"recent_records":refs,"recent_limit":20,"full_bibliography_query":ref_url}
    result={"source_id":target,"lookup_status":"identified" if obj else "unresolved", "method":method if obj else None,
            "record":obj,"original_simbad_response":receipt,"original_ambiguity_context":archived_context,
            "attached_identifiers":[r["id"] for r in xids if r.get("id")],"input_position_deg":[candidate["ra"],candidate["dec"]],"queried_identifiers":aliases,
            "record_url":link,"bibliography":bib,"checks":checks,"evidence_relative_to_survey_evidence":f"literature/{input_path.name}",
            "original_evidence_sha256":hashlib.sha256(input_path.read_bytes()).hexdigest(),
            "retrievals_relative_to_output":target+"/queries.jsonl"}
    dump(directory/target/"result.json",result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("evidence",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    files={}
    for path in (args.evidence/"literature").glob("candidate-*.json"):
        candidate=json.loads(path.read_text()).get("candidate",{})
        sid=candidate.get("source_id")
        if sid in TARGETS:
            files[sid]=path
    if set(files)!=set(TARGETS):
        raise ValueError(f"Missing original target evidence: {set(TARGETS)-set(files)}")
    with ThreadPoolExecutor(max_workers=2) as pool:
        outputs=list(pool.map(lambda sid:lookup(sid,files[sid],args.output),TARGETS))
    dump(args.output/"summary.json",{"scope":"SIMBAD source identity and attached bibliography; a failed service query is retained as unresolved and cannot establish novelty.",
                                    "targets":outputs})
    for row in outputs:
        print(row["source_id"],row["lookup_status"],row["method"],(row.get("bibliography") or {}).get("reported_count"),flush=True)


if __name__=="__main__":
    main()
