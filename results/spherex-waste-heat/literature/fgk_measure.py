"""Measure the CatWISE offsets and final-release NEOWISE visits at FGK targets."""
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from statistics import median
from fgk_baseline import ROOT, table
from fgk_fetch import position, read_result

def rows(name):
    with (ROOT / name).open() as stream:
        return list(csv.DictReader(stream))

def load_verified(kind):
    archive = ROOT / f"fgk_{kind}_raw.jsonl.gz"
    receipts = {x["source_id"]: x for x in json.loads((ROOT / f"fgk_{kind}_receipts.json").read_text())}
    results = {}
    for line in gzip.decompress(archive.read_bytes()).decode().splitlines():
        packet = json.loads(line)
        sid = packet["source_id"]
        data = packet.get("response", packet.get("csv")).encode()
        if hashlib.sha256(data).hexdigest() != receipts[sid]["sha256"]:
            raise ValueError(f"Input checksum mismatch {kind} {sid}")
        results[sid] = read_result(kind, data)
        if len(results[sid]) != receipts[sid]["rows"]:
            raise ValueError(f"Input row count mismatch {kind} {sid}")
    return results

def sep(ra, dec, center):
    ra0, dec0 = center
    return 3600 * math.hypot(((float(ra) - ra0 + 180) % 360 - 180) * math.cos(math.radians(dec0)), float(dec) - dec0)

def number(raw):
    try:
        return float(raw)
    except (ValueError, TypeError):
        return float("nan")

def rounded(value):
    return round(value, 5) if math.isfinite(value) else ""

def catalog(target, baseline, neighbors):
    center = position(target, 2015.405)
    pairs = sorted((sep(r["ra_pm"], r["dec_pm"], center), i, r) for i,r in enumerate(neighbors)
                   if math.isfinite(number(r["ra_pm"])) and math.isfinite(number(r["dec_pm"])))
    near = [(dist, r) for dist, _, r in pairs if dist <= 2.]
    summary = dict(source_id=target["source_id"], n_catwise_12arcsec=len(neighbors),
                   n_catwise_2arcsec=len(near), n_catwise_6arcsec=sum(p[0] <= 6 for p in pairs),
                   catwise_source="", offset_arcsec="", ab_flags="", cc_flags="")
    obs = near[0][1] if len(near)==1 else {}
    if obs:
        summary.update(catwise_source=obs["source_name"], offset_arcsec=rounded(near[0][0]), ab_flags=obs["ab_flags"], cc_flags=obs["cc_flags"])
    for b,zp in [(1,.004),(2,.032)]:
        pref = f"w{b}"
        mag = number(obs.get(pref+"mpro_pm"))
        err = number(obs.get(pref+"sigmpro_pm"))
        good = obs.get("ab_flags","")[b-1:b]=="0" and obs.get("cc_flags","")[b-1:b]=="0" and number(obs.get(pref+"snr_pm")) >=5 and err>0
        pred = float(baseline[f"pred_w{b}_native"]) - zp
        unw = float(baseline[f"unwise_w{b}_native"]) - zp
        summary.update({f"catwise_{pref}": rounded(mag), f"catwise_{pref}_error": rounded(err),
                        f"{pref}_good_match": bool(good), f"pred_minus_catwise_{pref}": rounded(pred - mag),
                        f"unwise_minus_catwise_{pref}": rounded(unw - mag)})
    return summary

def read_visits(target, detections):
    frames = Counter()
    points = defaultdict(list)
    accepted = []
    for r in detections:
        mjd = number(r["mjd"])
        center = position(target, 2000 + (mjd - 51544.5) / 365.25)
        if sep(r["ra"], r["dec"], center) <=2.:
            frames[r["scan_id"], r["frame_num"]]+=1
            accepted.append(r)
    for r in accepted:
        if frames[r["scan_id"], r["frame_num"]] != 1 or number(r["qual_frame"]) <= 5 or number(r["qi_fact"]) <=.5 or number(r["saa_sep"]) <=5 or number(r["nb"]) != 1:
            continue
        for b in (1, 2):
            pos = b-1
            pref=f"w{b}"
            mag=number(r[pref+"mpro"])
            err=number(r[pref+"sigmpro"])
            if r["moon_masked"][pos:pos+1]!="0" or r["cc_flags"][pos:pos+1]!="0" or r["ph_qual"][pos:pos+1] not in ("A","B"):
                continue
            if all([number(r[pref+"sat"])==0, number(r[pref+"snr"])>=5, number(r[pref+"rchi2"])<=3, err>0, math.isfinite(mag)]):
                points[b].append((number(r["mjd"]), mag, err))
    epochs = []
    for b, observations in points.items():
        visits = []
        for point in sorted(observations):
            if not visits or point[0]-visits[-1][-1][0] > 30:
                visits.append([])
            visits[-1].append(point)
        for visit in visits:
            if len(visit)<4:
                continue
            vals = [r[1] for r in visit]
            m = median(vals)
            sigma = max(1.2533*median([r[2] for r in visit])/math.sqrt(len(vals)),
                        1.2533*1.4826*median([abs(v-m) for v in vals])/math.sqrt(len(vals)))
            epochs.append(dict(source_id=target["source_id"], band=b, n=len(visit), mjd_first=round(visit[0][0],5),
                               mjd_last=round(visit[-1][0],5), median_magnitude=rounded(m), sigma_mag=rounded(math.hypot(sigma,.03))))
    return epochs

def main():
    astro = rows("fgk_astrometry.csv")
    baseline = {r["source_id"]: r for r in rows("fgk_baseline.csv")}
    cat, neo = load_verified("catwise"), load_verified("neo")
    if set(cat) != set(neo) or set(neo) != set(baseline) or len(astro) != 53:
        raise ValueError("Incomplete target set")
    catrows, all_visits, flags = [], [], []
    for t in astro:
        sid=t["source_id"]
        catrows.append(catalog(t, baseline[sid], cat[sid]))
        visits=read_visits(t,neo[sid])
        all_visits += visits
        record = dict(source_id=sid,n_raw_neo=len(neo[sid]))
        for b in (1,2):
            v = [x for x in visits if x["band"]==b]
            record[f"w{b}_visits"] = len(v)
            record[f"w{b}_good_exposures_in_visits"] = sum(x["n"] for x in v)
            delta = v[-1]["median_magnitude"]-v[0]["median_magnitude"] if len(v)>1 else float("nan")
            sig = math.hypot(v[-1]["sigma_mag"], v[0]["sigma_mag"]) if len(v)>1 else float("nan")
            record[f"w{b}_last_minus_first_mag"] = rounded(delta)
            record[f"w{b}_delta_significance"] = rounded(delta/sig)
        flags.append(record)
    table(ROOT/"fgk_catwise_comparison.csv",catrows)
    table(ROOT/"fgk_neo_visits.csv",all_visits)
    table(ROOT/"fgk_neo_coverage.csv",flags)
    clean=[r for r in catrows if r["w1_good_match"] and r["w2_good_match"]]
    isolated=[r for r in clean if r["n_catwise_6arcsec"]==1]
    def summaries(data):
        return dict(n=len(data), n_both_brighter_0p10=sum(all(r[f"pred_minus_catwise_w{b}"]>.10 for b in (1,2)) for r in data),
                    median_pred_minus_catwise=[round(median(r[f"pred_minus_catwise_w{b}"] for r in data),3) for b in (1,2)],
                    median_unwise_minus_catwise=[round(median(r[f"unwise_minus_catwise_w{b}"] for r in data),3) for b in (1,2)])
    def count(b,n):
        return sum(r[f"w{b}_visits"]>=n for r in flags)
    stats = dict(clean_catwise=summaries(clean), isolated_catwise=summaries(isolated),neo_raw_extractions=sum(map(len,neo.values())),
                 neo_stars_with_visit=[count(b,1) for b in (1,2)],neo_stars_with_two_visits=[count(b,2) for b in (1,2)],
                 neo_first_last_gt5sigma=[sum(abs(number(r[f"w{b}_delta_significance"]))>=5 for r in flags) for b in (1,2)])
    (ROOT/"fgk_statistics.json").write_text(json.dumps(stats,indent=2)+"\n")

if __name__ == "__main__":
    main()
