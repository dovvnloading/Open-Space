"""Independent test of printed NICER MRT calendar times on ordinary bursts."""
import argparse
import csv
from datetime import datetime
import hashlib
from pathlib import Path
import urllib.request

import numpy as np
from astropy.coordinates import get_body_barycentric
from astropy.io import fits
from astropy.time import Time, TimeDelta

HERE = Path(__file__).resolve().parent
OFFSET = 69.184
EXCLUDE = {"1E 1048.1-5937", "CXOU J010043.1-721134"}


def selected(mrt):
    rows = []
    for line in mrt.read_text().splitlines():
        if line[26:36].isdigit():
            rows.append((line[:21].strip(), line[26:36], line[37:60],
                         float(line[61:69]), float(line[120:125]), int(line[22:25])))
    chosen = []
    for source in sorted({row[0] for row in rows} - EXCLUDE):
        target = [r for r in rows if r[0] == source]
        eligible = [r for r in target if 20 <= r[3] <= 300 and sum(
            r[1] == q[1] and abs((datetime.fromisoformat(r[2])-datetime.fromisoformat(q[2])).total_seconds()) <= 150
            for q in target) == 1]
        if eligible:
            chosen.append(max(eligible,key=lambda r:r[4]))
    return chosen


def raw_prediction(label, ra, dec, header, orbit):
    # Earth and spacecraft Rømer delay, projected onto the source direction.
    # The built-in ephemeris is sufficient for a seconds-scale comparison.
    a, d = np.radians([ra, dec])
    n = np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
    epoch = Time(header["MJDREFI"],header["MJDREFF"],format="mjd",scale="tt")
    target = Time(label,scale="tdb")
    guess = (target.tt-epoch).sec
    ot = orbit.data["TIME"] + orbit.header.get("TIMEZERO",0)
    for _ in range(3):
        at = epoch + TimeDelta(guess,format="sec")
        xyz = [np.interp(guess,ot,orbit.data[k]) for k in "XYZ"]
        bary = get_body_barycentric("earth",at).xyz.to_value("m") + xyz
        predicted = at.tdb + TimeDelta(np.dot(bary,n)/299792458.,format="sec")
        guess += (target-predicted).sec
    # The uncorrected first guess can sit outside an orbit file near a
    # visit boundary. Only the converged spacecraft time needs bracketing.
    if not ot[0] <= guess <= ot[-1]:
        raise ValueError("Orbit does not bracket converged prediction")
    return float(guess-header.get("TIMEZERO",0))


def overlap(start, stop, gti):
    return float(np.maximum(0.,np.minimum(gti[:,1],stop)-np.maximum(gti[:,0],start)).sum())


def observations(cache, selected_row):
    src, obs, label, t90, log_flux, number = selected_row
    with fits.open(cache/f"ni{obs}_0mpu7_cl.evt.gz") as h, fits.open(cache/f"ni{obs}.orb.gz") as o:
        tab, hdr = h["EVENTS"].data, h["EVENTS"].header
        pred = raw_prediction(label,hdr["RA_OBJ"],hdr["DEC_OBJ"],hdr,o[1])
        gt = h["GTI"]
        gti = np.stack([gt.data["START"],gt.data["STOP"]],axis=1) + gt.header.get("TIMEZERO",0)-hdr.get("TIMEZERO",0)
        photons = tab["TIME"][(tab["PI"]>=50)&(tab["PI"]<=800)]
        sub = np.sort(photons[(photons>=pred-90)&(photons<=pred+90)])-pred
    # Completely blind maximum-count 250 ms window in ±90 s, with a
    # deterministic earliest-first tie break. No alignment to either prediction.
    width = 0.25
    starts = sub[sub<=90-width]
    counts = np.searchsorted(sub,starts+width,"left")-np.searchsorted(sub,starts,"left")
    peak = int(np.argmax(counts)) if len(counts) else None
    peak_start = round(float(starts[peak]),6) if peak is not None else ""
    peak_counts = int(counts[peak]) if peak is not None else 0
    def nt(t):
        return int(np.count_nonzero((sub>=t-.4)&(sub<t+.4)))
    return dict(source=src,obsid=obs,burst=number,start_label=label,t90_ms=t90,logFBB=log_flux,
        target_ra=hdr["RA_OBJ"],target_dec=hdr["DEC_OBJ"],raw_literal_s=round(pred,6),
        offset_peak_s=peak_start,peak_counts=peak_counts,counts_literal=nt(0),
        counts_shifted=nt(OFFSET),gti_literal_s=round(overlap(pred-.4,pred+.4,gti),6),
        gti_shifted_s=round(overlap(pred+OFFSET-.4,pred+OFFSET+.4,gti),6),
        gti_search_s=round(overlap(pred-90,pred+90,gti),6))


def replay(cache):
    hashes = list(csv.DictReader((HERE/"inputs.csv").open()))
    cache.mkdir(exist_ok=True,parents=True)
    for row in hashes:
        dest = cache/row["file"]
        if not dest.exists():
            with urllib.request.urlopen(row["url"],timeout=120) as r:
                data = r.read()
                if hashlib.sha256(data).hexdigest()!=row["sha256"]:
                    raise ValueError(f"Archive hash changed: {dest.name}")
                dest.write_bytes(data)
        if hashlib.sha256(dest.read_bytes()).hexdigest()!=row["sha256"]:
            raise ValueError(f"Archive hash changed: {dest.name}")
    return [observations(cache,r) for r in selected(cache/"pubtable.txt")]


def gbm_matches(cache, rows):
    from lxml import html
    events = []
    for tr in html.parse(str(cache/"gbm-catalog.html")).xpath("//tr[td]"):
        c = [" ".join(td.text_content().split()) for td in tr.xpath("td")]
        if len(c)>9:
            try:
                utc = datetime.strptime(c[1], "%y%m%d %H:%M:%S.%f")
            except ValueError:
                continue
            events.append((utc,c[2],c[8],float(c[0])))
    assert len(events)==1254
    matches = []
    for r in rows:
        obs = r["obsid"]
        hdr = fits.getheader(cache/f"ni{obs}_0mpu7_cl.evt.gz",1)
        epoch = Time(hdr["MJDREFI"],hdr["MJDREFF"],format="mjd",scale="tt")
        raw = float(r["raw_literal_s"])+hdr["TIMEZERO"]
        literal = (epoch+TimeDelta(raw,format="sec")).utc.datetime
        shifted = (epoch+TimeDelta(raw+OFFSET,format="sec")).utc.datetime
        for event,trigger,source,met in events:
            direct = (event-literal).total_seconds()
            diff = (event-shifted).total_seconds()
            if min(abs(direct),abs(diff))>1:
                continue
            # Inter-spacecraft propagation: nearer the source arrives earlier.
            with fits.open(cache/f"ni{obs}.orb.gz") as o:
                od = o[1].data
                ots = od["TIME"]+o[1].header.get("TIMEZERO",0)
                nicer_xyz = np.array([np.interp(raw+OFFSET,ots,od[k]) for k in "XYZ"])
            with fits.open(cache/f"glg_trigdat_all_bn{trigger}_v01.fit") as f:
                data = f["EVNTRATE"].data
                k = int(np.argmin(np.abs(data["TIME"]-met)))
                fermi_xyz = np.asarray(data["EIC"][k],float)*1000
            a,d = np.radians([r["target_ra"],r["target_dec"]])
            n = np.array([np.cos(a)*np.cos(d),np.sin(a)*np.cos(d),np.sin(d)])
            geometric = float(np.dot(n,nicer_xyz-fermi_xyz)/299792458)
            matches.append(dict(obsid=obs,source=source,trigger=trigger,gbm_event_utc=event.isoformat(timespec="milliseconds"),
                nicer_shifted_utc=shifted.isoformat(timespec="milliseconds"),
                fermi_minus_literal_s=round(direct,6),fermi_minus_shifted_s=round(diff,6),
                expected_fermi_minus_nicer_s=round(geometric,6),
                residual_after_geometry_s=round(diff-geometric,6)))
    return matches


def external_table(mrt):
    entries = [(line[26:36],int(line[22:25]),datetime.fromisoformat(line[37:60]))
               for line in mrt.read_text().splitlines() if line.startswith("SGR 1935+2154")]
    with (HERE/"published-times.csv").open() as f:
        for event in csv.DictReader(f):
            epoch = datetime.fromisoformat(event["epoch_tdb"])
            candidates = [(abs((epoch-t).total_seconds()),(epoch-t).total_seconds(),num)
                          for obs,num,t in entries if obs==event["obsid"]]
            best = min(candidates, default=None)
            if best and best[0]<=90:
                print(f"Borghese {event['epoch_tdb']} minus Chu #{best[2]}: {best[1]:.3f} s")
            else:
                print(f"Borghese {event['epoch_tdb']}: no same-observation match within 90 s")


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cache",type=Path,help="Directory for downloaded orbit, event and MRT files")
    args=parser.parse_args()
    out=replay(args.cache)
    with (HERE/"results.csv").open() as f:
        reference=list(csv.DictReader(f))
    external_table(args.cache/"pubtable.txt")
    external = gbm_matches(args.cache,out)
    with (HERE/"gbm-matches.csv").open() as f:
        saved = list(csv.DictReader(f))
    for match, frozen in zip(external,saved,strict=True):
        if not all(str(v)==frozen[k] for k,v in match.items()):
            raise ValueError(f"GBM timing comparison changed for {match['obsid']}")
        print(match["source"],match["trigger"],match["residual_after_geometry_s"])
    for row, frozen in zip(out,reference,strict=True):
        for k,v in row.items():
            if isinstance(v,(int,float)):
                if abs(v-float(frozen[k]))>0.003:
                    raise ValueError(f"{k} differs for {row['obsid']}")
            elif v!=frozen[k]:
                raise ValueError(f"{k} differs for {row['obsid']}")
        print(row['source'],row['obsid'],row['offset_peak_s'],row['counts_literal'],row['counts_shifted'])
