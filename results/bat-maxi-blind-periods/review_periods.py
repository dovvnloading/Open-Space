"""Index existing evidence and compare candidate frequencies with published clocks.

Run this postselection audit in Actions. Optional --manual JSON is a list of
{"source_id": "...", "period_days": <candidate>, "disposition": "...",
 "rationale": "...", "references": [{"title":"...", "url":"..."}],
 "published": [{"period_days":..., "uncertainty_days":..., "title":"...", "url":"..."}]}.
Omit the candidate period_days to supply published values to all periods of a source.
Dispositions require a candidate period, rationale and reference.
"""

import argparse
import gzip
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.parse import quote


# P, half-range in P (days), largest positive window harmonic.
WINDOWS = {"solar day":(1.,0.,128), "sidereal day":(.99726957,0.,128),
           "synodic Moon":(29.53059,.02,4), "sidereal Moon":(27.32166,.02,4),
           "anomalistic Moon":(27.55455,.02,4), "draconic Moon":(27.21222,.02,4),
           "Swift precession":(53.,2.,4), "ISS precession":(70.,2.,4),
           "year":(365.256,.25,4)}
HARMONICS = tuple(1/n for n in range(8,1,-1)) + tuple(float(n) for n in range(1,9))
ATOM = "{http://www.w3.org/2005/Atom}"
CLOCK = re.compile(r"orbit|superorb|period|spin|pulsat|recurr|modulat|precess|cycle|ephemer",re.I)
VALUE = re.compile(r"(?<![\w.])(\d+(?:\.\d+)?)(?:\((\d+)\)|\s*(?:±|\+/-)\s*(\d+(?:\.\d+)?))?\s*"
                   r"(?:-\s*)?(days?\b|d\b|hours?\b|hrs?\b|h\b|minutes?\b|mins?\b|"
                   r"seconds?\b|secs?\b|s\b|years?\b|yrs?\b)", re.I)


def norm(x):
    return re.sub(r"[^a-z0-9+-]", "", str(x).lower().replace("−","-").replace("–","-"))


def lines(path):
    fn = gzip.open if path.suffix==".gz" else open
    with fn(path,"rt") as f:
        for line in f:
            if line.strip():
                yield json.loads(line)


def abstracts(item, directory, cache):
    pubs = item.get("bibliography",{})
    for name, query in pubs.get("queries",{}).items():
        ref = query.get("query",{}).get("sha256")
        if not ref:
            continue
        if ref not in cache:
            raw = gzip.decompress((directory/"responses"/(ref+".gz")).read_bytes())
            if hashlib.sha256(raw).hexdigest()!=ref:
                raise ValueError("Literature input SHA-256 mismatch: "+ref)
            if name.startswith("arXiv"):
                entries = [(e.findtext(ATOM+"id",""),e.findtext(ATOM+"summary",""))
                           for e in ET.fromstring(raw).findall(ATOM+"entry")]
            elif name.startswith("Crossref"):
                entries = [(w.get("DOI",""),w.get("abstract",""))
                           for w in json.loads(raw)["message"]["items"]]
            elif name.startswith("ADS"):
                entries = [(w["bibcode"],w.get("abstract","")) for w in json.loads(raw)["response"]["docs"]]
            else:
                entries = []
            cache[ref] = entries
        for key, abstract in cache[ref]:
            paper = pubs.get("papers",{}).get(key)
            if paper:
                yield paper, abstract, ref


def values(item, directory, cache, manual):
    for table in item.get("catalogues",[]):
        for p in table.get("periods",[]):
            if not p.get("period_days"):
                continue
            number = p["period_days"]
            # Published catalogue uncertainties are used when provided as standard
            # e_* fields. In their absence the precision of the printed value is
            # kept as a rounding interval, not treated as a measured error.
            f = p["field"]
            v = float(p["value"])
            if not math.isfinite(v) or v<=0:
                continue
            unit = number/v
            es = table.get("row",{}).get("e_"+f)
            try:
                sigma = abs(float(es))*unit
            except (TypeError,ValueError):
                sigma = 0.
            kind = ("catalogue value with unit under review" if
                    table["vizier"]=="J/A+A/684/A124" and f=="Period" else "catalogue")
            yield dict(period_days=number, uncertainty_days=sigma, rounding_days=half_digit(p["value"])*unit,
                       kind=kind, title=table["family"]+" / "+table["table"]+" / "+f,
                       url="https://vizier.cds.unistra.fr/viz-bin/VizieR?-source="+quote(table["table"],safe=""),
                       context=p["description"], association="identifier" if table["matched_identifiers"] else "coordinate only")
    aliases = [norm(v) for v in item["candidate"].get("aliases",[])+[item["candidate"]["source"]]]
    for paper, abstract, digest in abstracts(item,directory,cache):
        title = paper["title"]
        text = re.sub(r"<[^>]+>", " ", title+" "+(abstract or ""))
        text = " ".join(text.replace("\\pm","±").replace("$","").replace("−","-").split())
        if not paper.get("alias_hits") and not any(a in norm(text) for a in aliases if len(a)>=5):
            continue
        for m in VALUE.finditer(text):
            context = text[max(0,m.start()-100):m.end()+100]
            if not CLOCK.search(context):
                continue
            unit = m[4].lower()[0]
            scale = {"d":1.,"h":1/24,"m":1/1440,"s":1/86400,"y":365.25}[unit]
            q = half_digit(m[1])
            err = float(m[3]) if m[3] else 2*q*int(m[2]) if m[2] else 0.
            v = float(m[1])
            if v > 0:
                yield dict(period_days=v*scale,uncertainty_days=err*scale,rounding_days=q*scale,
                           title=title,url=paper.get("url",""),kind="abstract lead",association="source in title or abstract",
                           context=" ".join(text[:m.start()].split()[-7:])+" "+m[0]+" "+
                                   " ".join(text[m.end():].split()[:7]),query_sha256=digest)
    pubs = item.get("bibliography",{})
    for scan in pubs.get("manuscript_checks",[]):
        ident = scan["manuscript"]
        paper = next((p for k,p in pubs.get("papers",{}).items() if k.endswith("/"+ident)),{})
        title = paper.get("title",ident)
        for page in scan.get("source_pages",[]):
            reference = dict(kind="manuscript lead",title=f"{title} (p. {page['page']})",
                             url="https://arxiv.org/abs/"+ident,
                             pdf_sha256=scan.get("retrieval",{}).get("sha256"))
            for line in page.get("source_context",[]):
                for match in VALUE.finditer(line):
                    scale = {"d":1.,"h":1/24,"m":1/1440,"s":1/86400,"y":365.25}[match[4].lower()[0]]
                    val,rounded = float(match[1])*scale,half_digit(match[1])*scale
                    err = float(match[3])*scale if match[3] else 2*rounded*int(match[2]) if match[2] else 0.
                    if val>0:
                        yield dict(reference,period_days=val,rounding_days=rounded,uncertainty_days=err,
                                   association="same printed line as source",context=line)
            for hit in page.get("unit_bearing_periods_on_same_page",[]):
                if hit.get("period_days",0)>0:
                    yield dict(reference,period_days=hit["period_days"],rounding_days=0.,uncertainty_days=0.,
                               association="same printed page; ownership requires verification",context=hit["text"])
    for row in manual:
        for p in row.get("published",[]):
            if not p.get("title") or not p.get("url") or not p.get("period_days"):
                raise ValueError("Manual published period needs a primary reference and numerical value")
            yield dict(kind="curated publication",association="manual primary-paper identification",
                       title=p["title"],url=p["url"],period_days=float(p["period_days"]),
                       rounding_days=float(p.get("rounding_days",0) or 0),
                       uncertainty_days=float(p.get("uncertainty_days",0) or 0),context=p.get("context",""))



def half_digit(value):
    mantissa,*exponent = str(value).lower().split("e")
    decimals = len(mantissa.split(".")[1]) if "." in mantissa else 0
    return .5*10**((int(exponent[0]) if exponent else 0)-decimals)


def frequency_bounds(period,half_width):
    if not math.isfinite(period+half_width) or period<=0 or half_width<0 or half_width>=period:
        return None
    return 1/(period+half_width),1/(period-half_width)


def overlaps(interval,bands):
    a,b = interval
    return all(a<=hi and b>=lo for lo,hi in bands)


def rectified(a,b):
    if a<=0<=b:
        return (0.,max(abs(a),abs(b)))
    return tuple(sorted((abs(a),abs(b))))


def window_models():
    # Generate every specified daily clock harmonic, even for a high-frequency
    # known clock whose alias folds into the daily product's long-period band.
    yield "intrinsic",0,0.,0.,0.
    for name,(period,width,max_harmonic) in WINDOWS.items():
        lo,hi = frequency_bounds(period,width)
        for n in range(1,max_harmonic+1):
            for sign in (-1,1):
                k = n*sign
                a,b = sorted((k*lo,k*hi))
                yield name,k,k/period,a,b


def compare(p,bands,f):
    period = p["period_days"]
    pub = frequency_bounds(period,max(3*p.get("uncertainty_days",0),p.get("rounding_days",0)))
    if pub is None:
        return []
    matches = []
    for h in HARMONICS:
        for name,k,shift,left,right in window_models():
            predicted = rectified(h*pub[0]+left, h*pub[1]+right)
            if not overlaps(predicted,bands):
                continue
            central = abs(h/period+shift)
            residual = abs(central-f)
            width = predicted[1]-predicted[0]
            matches.append(dict(category=("direct" if h==1 else "harmonic") if k==0 else "window alias",
                                carrier_harmonic=h,window=name,window_harmonic=k,central_alias_per_day=central,
                                alias_range_per_day=predicted,difference_per_day=residual,
                                width_over_candidate_frequency=width/f,
                                alias_period_days=1/central if central>0 else None))
    return matches


def escape(v):
    return str(v).replace("|","\\|").replace("\n"," ")


def write_index(record):
    c = record["candidate"]
    lines = [f"### {escape(c['source'])} · {c['period_days']:.8g} days", ""]
    if record.get("manual_decision"):
        m = record["manual_decision"]
        refs = ", ".join(f"[{escape(r['title'])}]({r['url']})" for r in m["references"])
        lines += [f"**{escape(m['disposition'])}.** {escape(m['rationale'])} {refs}", ""]
    checks = ", ".join(record["incomplete_checks"]) or "no missing catalogue groups or literature services flagged by the retrieval checks"
    lines += [f"Retrieval coverage — {escape(checks)}. [Source evidence]({record['detail_link']})", ""]
    if not record["reported_matches"]:
        return lines+["No source-linked numeric value matched the frequency intervals in this audit.", ""]
    lines += ["| Prior value (days) | Comparison | Source |", "| ---: | --- | --- |"]
    for x in record["reported_matches"]:
        p,m = x["publication"],x["match"]
        if m["category"]=="window alias":
            comp=f"{m['carrier_harmonic']:.3g} × published frequency {m['window_harmonic']:+d} × {m['window']} frequency"
        else:
            comp=f"{m['category']}, {m['carrier_harmonic']:.3g} × published frequency"
        if m["width_over_candidate_frequency"]>.1:
            comp+=" (broad prior interval)"
        lines.append(f"| {p['period_days']:.9g} | {escape(comp)} | {escape(p['kind'])} · [{escape(p['title'])}]({p['url']}) |")
    return lines+[""]


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evidence",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--manual",type=Path)
    args=ap.parse_args()
    root,dest=args.evidence,args.output
    dest.mkdir(parents=True,exist_ok=True)
    lit=root/"literature"
    manual=json.loads(args.manual.read_text()) if args.manual else []
    if not isinstance(manual,list):
        raise ValueError("--manual expects a JSON list")
    for m in manual:
        if m.get("decision"):
            m["disposition"] = m["decision"]
            m["rationale"] = m["explanation"]
        for ref in m.get("references",[]):
            if "label" in ref and "title" not in ref:
                ref["title"] = ref["label"]
        periods = list(m.get("published",[]))
        periods.extend(ref for ref in m.get("references",[]) if "period_days" in ref)
        m["published"] = periods
    seen_manual=set()
    spans={r["roster_id"]:{k:b.get("span_days") for k,b in r.get("instruments",{}).items()}
           for path in root.glob("search-*.jsonl.gz") for r in lines(path)}
    controls_path=root/"validation"/"controls.json"
    controls=json.loads(controls_path.read_text()) if controls_path.is_file() else []
    heading=["# Published-period and alias review", "",
             "Every row is a lead for source-by-source review. The ledger preserves all matching catalogue and abstract values. A source name and a number in the same abstract may refer to different components, objects or timescales.", "",
             "A match must overlap each instrument's refined frequency interval (one Rayleigh width plus its three-standard-error interval). Published frequencies use the larger of three reported standard errors and half the final printed digit. The test includes direct carriers, source harmonics 2–8, source subharmonics 1/2–1/8, solar and sidereal daily beats through order 128, and lunar, precession and annual sidebands through order 4. Each combination uses one window family at a time. Very wide published intervals or high-order sidebands readily generate chance matches. Coherent leakage strength needs the actual exposures and source spectrum.", ""]
    index=list(heading)
    cache={}
    with gzip.open(dest/"matches.jsonl.gz","wt") as records:
        for item in lines(lit/"triage.jsonl"):
            detail=json.loads((lit/item["detail"]).read_text())
            cand=detail["candidate"]
            f=1/cand["period_days"]
            sm=[]
            for n,m in enumerate(manual):
                if m["source_id"]!=cand["source_id"]:
                    continue
                if "period_days" in m and not math.isclose(float(m["period_days"]),cand["period_days"],rel_tol=1e-8):
                    continue
                sm.append(m)
                seen_manual.add(n)
            decisions=[m for m in sm if m.get("disposition")]
            if len(decisions)>1 or any("period_days" not in m or not m.get("references") or not m.get("rationale") for m in decisions):
                raise ValueError("Each period needs at most one reference-backed, frequency-specific manual disposition")
            bands=[]
            for band in ("bat","maxi"):
                b=cand[band]
                span=spans[cand["source_id"]][band]
                interval=frequency_bounds(b["period_days"],3*b["uncertainty_days"])
                if not interval or not span or span<=0:
                    raise ValueError("Missing finite frequency interval for "+cand["source"])
                bands.append([max(0,interval[0]-1/span),interval[1]+1/span])
            refs={}
            for p in values(detail,lit,cache,sm):
                refs[(p["kind"],p["title"],p["period_days"],p.get("uncertainty_days",0),p.get("rounding_days",0))]=p
            for c in controls:
                if c.get("roster_id")!=cand["source_id"]:
                    continue
                val=c["reference_period_days"]
                prior=dict(kind="published control",association="named literature control",title=c["reference"],
                           url=c["url"],period_days=val,uncertainty_days=0.,rounding_days=half_digit(val),
                           context=c["known_cycle"]+" cycle of "+c["source"])
                refs[(prior["kind"],prior["title"],val,0.,prior["rounding_days"])]=prior
            matches=[]
            for p in refs.values():
                for m in compare(p,bands,f):
                    matches.append(dict(publication=p,match=m))
            cat_order={"direct":0,"harmonic":1,"window alias":2}
            matches.sort(key=lambda v:(cat_order[v["match"]["category"]],v["publication"]["kind"]=="catalogue value with unit under review",v["publication"]["kind"] in ("abstract lead","manuscript lead"),
                         v["match"]["width_over_candidate_frequency"]>.1,
                         v["match"]["carrier_harmonic"]!=1,abs(v["match"]["window_harmonic"]),
                         v["match"]["difference_per_day"]))
            short=[]
            used=set()
            for row in matches:
                p,m=row["publication"],row["match"]
                # Show a concise view of the distinct evidence, retaining the full
                # source, interval and order for every record in the compressed ledger.
                key=(m["category"],round(p["period_days"],8),p["kind"])
                if key not in used and sum(v["match"]["category"]==m["category"] for v in short)<3:
                    short.append(row)
                    used.add(key)
            rec=dict(candidate={k:cand[k] for k in ("source","source_id","period_days","uncertainty_days")},
                     detail_link=str(Path("..")/"literature"/item["detail"]),manual_decision=decisions[0] if decisions else None,
                     incomplete_checks=detail.get("incomplete",[]),band_frequency_intervals_per_day=dict(zip(("bat","maxi"),bands)),
                     numeric_values_considered=len(refs),publications=list(refs.values()),reported_matches=short,matches=matches)
            records.write(json.dumps(rec,allow_nan=False)+"\n")
            index.extend(write_index(rec))
    if set(range(len(manual)))!=seen_manual:
        raise ValueError("Unmatched manual records: "+str(set(range(len(manual)))-seen_manual))
    (dest/"index.md").write_text("\n".join(index)+"\n")


if __name__=="__main__":
    main()
