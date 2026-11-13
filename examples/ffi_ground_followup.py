"""Recheck the unmatched LMC-4 source against resolved ground-based photometry."""

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import tarfile
from urllib.parse import urlencode

import numpy as np
from astropy.timeseries import LombScargle

from open_space.tess_archive import atomic_bytes, download
from open_space.tess_eb_catalog import separation

SOURCE = "ICRS-77.05502-69.23710"
OGLE = "https://ogle.astrouw.edu.pl/ogle"
MACHO_QUERY = ("SELECT TOP 10000 dateobs,sideofpier,rmag,rerr,rcut,bmag,berr,bcut "
               "FROM public.photometry_view WHERE fieldid=1 AND tileid=4532 AND seqn=22 ORDER BY dateobs")
GAIA_ID = "4661229089008619008"


def tap(service, query):
    return service + "?" + urlencode({"REQUEST":"doQuery", "LANG":"ADQL", "FORMAT":"csv", "QUERY":query})


def read_source(cache, url):
    path = cache / hashlib.sha256(url.encode()).hexdigest()
    if not path.exists():
        atomic_bytes(path, download(url, 20_000_000))
    raw = path.read_bytes()
    return raw, {"url": url, "sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}


def search(time, mag, err, periods, near_range, grid=(.15,8), detrend=False):
    valid = np.isfinite(time) & np.isfinite(mag) & np.isfinite(err) & (err > 0) & (err < .2)
    t, y, e = time[valid], mag[valid], err[valid]
    if detrend:
        y = y - np.array([np.median(y[abs(t-at) < 12]) for at in t])
    ls = LombScargle(t, y, e)
    freq, power = ls.autopower(minimum_frequency=1/grid[1], maximum_frequency=1/grid[0],
                               samples_per_peak=12)
    near = (freq > 1/near_range[1]) & (freq < 1/near_range[0])
    return {"epochs": len(t), "baseline_days": float(np.ptp(t)), "detrended": detrend,
            "search_days": list(grid), "local_search_days": list(near_range),
            "global_peak_days": float(1/freq[np.argmax(power)]), "global_peak_power": float(power.max()),
            "local_peak_days": float(1/freq[near][np.argmax(power[near])]),
            "local_peak_power": float(power[near].max()),
            "power_at_sector_periods": [float(ls.power(1/p)) for p in periods]}


def build(field_result, cache):
    fields = json.loads(field_result)["fields"]
    src = next(s for field in fields for s in field["sources"] if s["source_id"] == SOURCE)
    periods = [src[k]["period_days"] for k in ("discovery", "validation")]
    repeated = lambda values: [search(*values, periods, (1.07,1.16), detrend=d) for d in (False,True)]
    results = {}
    url = OGLE + "/ogle2/dia/ac/lmc_sc11.tar.gz"
    data, receipt = read_source(cache, url)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        for name in ("OGLE050812.41-691410.1.mag", "OGLE050813.83-691422.9.mag"):
            payload = archive.extractfile(name).read()
            values = np.loadtxt(io.BytesIO(payload)).T
            result = repeated(values)
            halves = [search(*values[:,i], periods, (1.07,1.16), detrend=True)
                      for i in np.array_split(np.arange(values.shape[1]), 2)]
            results[name] = {"receipt": receipt, "member":name,
                             "member_sha256": hashlib.sha256(payload).hexdigest(),
                             "analysis": result, "detrended_time_halves": halves}
    name = "OGLE-LMC-LPV-20399"
    url = OGLE + "/ogle3/OIII-CVS/lmc/lpv/phot/I/" + name + ".dat"
    data, receipt = read_source(cache, url)
    results[name] = {"receipt":receipt,
                     "analysis":search(*np.loadtxt(io.BytesIO(data)).T,periods,(1.07,1.16),grid=(.15,500))}
    url = tap("https://machotap.asvo.nci.org.au/ncitap/tap/sync", MACHO_QUERY)
    data, receipt = read_source(cache, url)
    rows = list(csv.DictReader(io.StringIO(data.decode())))
    if not 1 <= len(rows) < 10000:
        raise ValueError("Empty or truncated MACHO response")
    item = {"receipt": receipt, "input_rows": len(rows), "bands": {}}
    for k in ("r","b"):
        good = [r for r in rows if r[k+"cut"]=="1" and r[k+"mag"] and r[k+"err"] and 0<float(r[k+"err"])<.08]
        t, y, e = np.array([[float(r[c]) for r in good] for c in ("dateobs", k+"mag", k+"err")])
        ok = abs(y-np.median(y)) < 6*1.4826*np.median(abs(y-np.median(y)))
        item["bands"][k] = repeated((t[ok],y[ok],e[ok]))
    results["MACHO 1.4532.22"] = item
    url = tap("https://gea.esac.esa.int/tap-server/tap/sync",
              "SELECT * FROM gaiadr3.vari_classifier_result WHERE source_id=" + GAIA_ID)
    data, receipt = read_source(cache, url)
    gaia = {"receipt": receipt, "history":list(csv.DictReader(io.StringIO(data.decode())))}
    url = "https://gea.esac.esa.int/data-server/data?" + urlencode({
        "RETRIEVAL_TYPE":"EPOCH_PHOTOMETRY","ID":GAIA_ID,"FORMAT":"CSV","BAND":"G","RELEASE":"Gaia DR3"})
    data, receipt = read_source(cache, url)
    rows = [r for r in csv.DictReader(io.StringIO(data.decode())) if r["g_transit_flux"]
            and r["g_transit_flux_error"] and r["variability_flag_g_reject"]=="false"
            and r["rejected_by_photometry"]=="false"]
    t,y,e=np.array([[float(r[c]) for r in rows] for c in ("g_transit_time","g_transit_flux","g_transit_flux_error")])
    ls=LombScargle(t,y,e)
    gaia["epoch_receipt"]=receipt
    gaia["weighted_power_at_exact_tess_periods"]=[float(ls.power(1/p)) for p in periods]
    # Same 30-arcsec cone as the independent catalog search.
    conepos = f"{src['ra_deg']:.9f},{src['dec_deg']:.9f},0.0083333333"
    query = ("SELECT main_id,otype,ra,dec FROM basic WHERE 1=CONTAINS("
             f"POINT('ICRS',ra,dec),CIRCLE('ICRS',{conepos}))")
    data, receipt = read_source(cache, tap("https://simbad.cfa.harvard.edu/simbad/sim-tap/sync", query))
    records = list(csv.DictReader(io.StringIO(data.decode())))
    for r in records:
        r["separations_arcsec"] = [separation(src[k]["ra_deg"],src[k]["dec_deg"],float(r["ra"]),float(r["dec"]))
                                    for k in ("discovery","validation")]
    return {"source_id": SOURCE, "field_input_sha256": hashlib.sha256(field_result).hexdigest(),
            "tess_periods_days": periods, "simbad": {"receipt": receipt,"objects":records},
            "gaia_classifier_and_targeted_search": gaia,"ground_based_sources":results,
            "macho_quality": "Per-band survey cut=1, 0<error<0.08 mag, six-MAD clipping",
            "conditioning": "Both raw magnitudes and values minus a centered 24-day running median",
            "local_search": "Post hoc range covering the two TESS periods and possible nearby ground-based maxima"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("field_result",type=Path)
    parser.add_argument("--cache",required=True,type=Path)
    args=parser.parse_args()
    print(json.dumps(build(args.field_result.read_bytes(),args.cache),indent=2,allow_nan=False))


if __name__ == "__main__":
    main()
