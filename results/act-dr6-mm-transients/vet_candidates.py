#!/usr/bin/env python3
"""Timestamped public counterpart-service *attempts* for every provisional ACT peak.

Output is deliberately unvetted: HTTP success is NOT an association or a proof
of absence. JPL SB_IDENT covers small bodies, limited Horizons query a bright
planet, NVSS radio alone does not cover VLASS/RACS, and TNS/GCN HTML is not an
authenticated all-alerts crossmatch. Actual FITS time pixels are required.
"""
from __future__ import annotations
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

from astropy.coordinates import SkyCoord
import astropy.units as u

from run_preregistered_field import BASE, selection, time_at_position, KNOWN_RA, KNOWN_DEC

MAX_BYTES = 160_000
TIMEOUT = 10
USER_AGENT = "Open-Space-ACT-vetting-provisional-field/1.0"


def probe(url):
    try:
        with urlopen(Request(url, headers={"User-Agent": USER_AGENT}), timeout=TIMEOUT) as response:
            raw = response.read(MAX_BYTES+1)
            if len(raw) > MAX_BYTES:
                return {"status": "too_large_not_interpreted", "http_status": response.status, "url": url}
            ans = {"status": "http_response_not_independently_interpreted",
                   "http_status": int(response.status), "url": url, "bytes": len(raw),
                   "response_sha256": sha256(raw).hexdigest()}
            # Never infer that an empty or successful HTTP response is a catalog veto.
            if "ssd-api.jpl.nasa.gov" in url:
                try:
                    j = json.loads(raw)
                    ans["jpl_signature"] = j.get("signature")
                    ans["small_body_second_pass_count"] = j.get("n_second_pass")
                    ans["jpl_api_error"] = j.get("message") or j.get("error")
                    ans["status"] = ("jpl_small_body_api_response" if "signature" in j and
                                    not j.get("message") and not j.get("error") else "jpl_response_unvalidated")
                    if "data_second_pass" in j:
                        ans["first_three_small_body_rows"] = j["data_second_pass"][:3]
                except (ValueError, TypeError):
                    ans["status"] = "jpl_non_json_response"
            if "horizons.api" in url:
                try:
                    j=json.loads(raw)
                    body=j.get("result","")
                    ans["horizons_error"] = j.get("error")
                    ans["horizons_has_ephemeris"] = "$$SOE" in body and "$$EOE" in body
                    if ans["horizons_has_ephemeris"]:
                        ans["horizons_ephemeris_lines"] = body.split("$$SOE")[1].split("$$EOE")[0].strip().splitlines()[:2]
                except (ValueError,TypeError):
                    ans["status"]="horizons_non_json_response"
            return ans
    except Exception as e:
        return {"status":"request_failed", "url":url, "error":f"{type(e).__name__}: {str(e)[:180]}"}


def urls(ra,dec,iso_time):
    pos = SkyCoord(ra=ra*u.deg, dec=dec*u.deg)
    ra_sexa = pos.ra.to_string(unit=u.hourangle, sep="-", precision=2, pad=True)
    dec_sexa = pos.dec.to_string(unit=u.deg, sep="-", precision=2, pad=True, alwayssign=True)
    coord = f"{ra:.8f} {dec:+.8f}"
    radial = f"{ra:.8f},{dec:.8f}"
    cone = f"CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',{ra:.8f},{dec:.8f},0.001388889))=1"
    gaia_q = ("SELECT TOP 8 source_id,ra,dec,phot_g_mean_mag,parallax FROM gaiadr3.gaia_source WHERE "
              +cone)
    lookups = {
        "gaia_dr3_tap_5arcsec": "https://gea.esac.esa.int/tap-server/tap/sync?"+
               urlencode({"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"json","QUERY":gaia_q}),
        "nvss_vizier_90arcsec": "https://vizier.cds.unistra.fr/viz-bin/asu-tsv?"+
               urlencode({"-source":"VIII/65/nvss","-c":coord,"-c.rs":"90","-out.max":"30"}),
        "simbad_2arcmin": "https://simbad.cds.unistra.fr/simbad/sim-coo?"+
               urlencode({"Coord":coord,"Radius":"2","Radius.unit":"arcmin","output.format":"ASCII"}),
        "gcn_public_html_not_conclusive": "https://gcn.nasa.gov/search?"+urlencode({"query":radial}),
        "tns_public_html_not_conclusive": "https://www.wis-tns.org/search?"+
               urlencode({"ra":str(ra),"decl":str(dec),"radius":"2","coords_unit":"arcmin","format":"html"})
    }
    if iso_time:
        dt = datetime.fromisoformat(iso_time).astimezone(timezone.utc)
        stamp = dt.strftime("%Y-%m-%d_%H:%M:%S")
        lookups["jpl_sb_ident_topocentric_2pass"] = "https://ssd-api.jpl.nasa.gov/sb_ident.api?"+urlencode({
            "lon":"-67.7875","lat":"-22.9586","alt":"5.19",
            "obs-time":stamp,"fov-ra-center":ra_sexa,"fov-dec-center":dec_sexa,
            "fov-ra-hwidth":"0.20","fov-dec-hwidth":"3.0",
            "two-pass":"true","suppress-first-pass":"true","req-elem":"false"})
        minute = dt.strftime("%Y-%m-%d %H:%M")
        end = dt.replace(second=0, microsecond=0).timestamp()+120
        stop = datetime.fromtimestamp(end,timezone.utc).strftime("%Y-%m-%d %H:%M")
        lookups["jpl_horizons_jupiter_only_not_all_objects"] = (
            "https://ssd.jpl.nasa.gov/api/horizons.api?"+urlencode({
                "format":"json","COMMAND":"'599'","MAKE_EPHEM":"'YES'",
                "EPHEM_TYPE":"'OBSERVER'","CENTER":"'coord@399'",
                "COORD_TYPE":"'GEODETIC'","SITE_COORD":"'-67.7875,-22.9586,5.19'",
                "START_TIME":f"'{minute}'","STOP_TIME":f"'{stop}'",
                "STEP_SIZE":"'1 m'","QUANTITIES":"'1'"}))
    return lookups


def one(item,rows):
    s,p=item
    row=rows.get(str(s["map_id"]))
    out={"map_id":s["map_id"],"ra_deg":p["ra_deg"],"dec_deg":p["dec_deg"],
         "nominal_white_snr":p["snr_white"],
         "input_is_a_provisional_peak_not_a_vetted_candidate":True,
         "association_verdict":"unvetted"}
    dist_arcmin = ((p["ra_deg"]-KNOWN_RA)**2*0.959 +
                   (p["dec_deg"]-KNOWN_DEC)**2)**0.5 *60
    out["angular_distance_to_published_ACT_2018_event_arcmin_approx"]=dist_arcmin
    if row is None:
        out["error"]="Frozen row missing"
        return out
    timing=time_at_position(row,p["ra_deg"],p["dec_deg"])
    out["time_from_actual_mean_time_fits"]=timing
    when=timing.get("utc") if timing.get("status")=="mean_time_fits_pixel_read" else None
    lookups=urls(p["ra_deg"],p["dec_deg"],when)
    out["services"]={k:probe(url) for k,url in lookups.items()}
    out["mandatory_gates_not_resolved"]=[
        "measured_pa4_depth1_beam_and_transfer",
        "scan_correlated_noise_and_trials",
        "multi_epoch_multi_array_independence",
        "Gaia_radio_matches_not_interpreted_as_associations",
        "GCN_TNS_search_pages_are_not_authoritative_broker_crossmatches",
        "JPL_SB_IDENT_and_Jupiter_only_do_not_exclude_all_moving_objects",
        "native_pixel_time_and_coverage_for_event_model",
    ]
    return out


def main():
    receipt=json.loads((BASE/"science-pilot-receipt.json").read_text())
    scans=(receipt.get("matched_filter") or {}).get("scans",[])
    if not scans or receipt.get("blocker"):
        raise ValueError("Cannot vet without completed science-screen receipt")
    tests,_=selection()
    rows={str(r["map_id"]):r for r in tests}
    all_items=[(s,p) for s in scans for p in s.get("provisional_peaks",[])]
    # Exact 9 local peaks in this fixed field; never silently discard peaks.
    if any(s.get("provisional_peaks_omitted_beyond_first_30",0) for s in scans):
        raise ValueError("Upstream truncated peaks; incomplete vetting census")
    with ThreadPoolExecutor(max_workers=5) as pool:
        results=list(pool.map(lambda it:one(it,rows),all_items))
    record={"classification":"provisional_candidate_service_queries_not_associations",
            "utc":datetime.now(timezone.utc).isoformat(),
            "source":"actual depth1 pixel screen: science-pilot-receipt.json",
            "count_of_peaks_examined":len(all_items),"no_peak_omitted":True,
            "actual_timing_pixel_attempted_for_each_peak":True,
            "published_control_is_known":"ACT-T J200758+160954; Naess et al. 2021",
            "horizons_only_jupiter_and_small_bodies_not_an_exhaustive_solar_system_veto":True,
            "vlass_racs_not_machine_crossmatched":True,
            "external_vetting_complete":False,"discovery":False,"rate_limit":None,
            "objects":results}
    (BASE/"candidate-vetting-receipt.json").write_text(json.dumps(record,sort_keys=True,indent=2,allow_nan=False)+"\n")
    print(json.dumps({k:v for k,v in record.items() if k!="objects"},indent=2))
    for item in results:
        print("candidate",item["map_id"],item["ra_deg"],item["dec_deg"],
              item["time_from_actual_mean_time_fits"].get("status"),
              {k:v["status"] for k,v in item.get("services",{}).items()},flush=True)


if __name__=="__main__":
    main()
