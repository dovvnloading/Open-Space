"""Replay independent public catalogue requests, retaining source responses."""
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen, Request

OUT = Path(__file__).resolve().parent / "catalogue_receipts"
IRSA = """SELECT designation, ra, dec, w1mpro, w1sigmpro, w1sat, w2mpro, w2sigmpro, w2sat, w3mpro, w3sigmpro, w4mpro, w4sigmpro, ph_qual, cc_flags, ext_flg
FROM allwise_p3as_psd
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',153.57985,21.07495,0.0030))"""
GAIA = """SELECT source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,phot_g_mean_mag,astrometric_params_solved,ruwe
FROM gaiadr3.gaia_source
WHERE 1=CONTAINS(POINT('ICRS',ra,dec), CIRCLE('ICRS',153.57985,21.07495,0.00834))"""
EXO = """SELECT DISTINCT hostname,pl_name,gaia_dr3_id,hip_name,tic_id FROM ps WHERE hip_name='HIP 50156' OR gaia_dr3_id='Gaia DR3 628827128175694976' OR tic_id='TIC 95339414' OR hostname IN ('MCC 124','HIP 50156','DK Leo','GJ 2079')"""
SERVICES = [
    ("allwise", "https://irsa.ipac.caltech.edu/TAP/sync", IRSA),
    ("twomass", "https://irsa.ipac.caltech.edu/TAP/sync",
     "SELECT designation,ra,dec,j_m,j_msigcom,h_m,h_msigcom,k_m,k_msigcom,ph_qual,cc_flg FROM fp_psc WHERE designation='10141918+2104297'"),
    ("gaia_dr2", "https://gea.esac.esa.int/tap-server/tap/sync",
     "SELECT source_id,ra,dec,parallax,parallax_error,pmra,pmra_error,pmdec,pmdec_error,astrometric_excess_noise FROM gaiadr2.gaia_source WHERE source_id=628827128175694976"),
    ("gaia_dr3", "https://gea.esac.esa.int/tap-server/tap/sync", GAIA),
    ("confirmed_planets", "https://exoplanetarchive.ipac.caltech.edu/TAP/sync", EXO),
]


def get_one(entry):
    name, base, query = entry
    params = dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", QUERY=query)
    url = base + "?" + urlencode(params)
    receipt = dict(query=query, url=url, started_utc=datetime.now(timezone.utc).isoformat())
    try:
        with urlopen(Request(url, headers={"User-Agent": "Open-Space-JWST-research/1.0"}), timeout=45) as r:
            body = r.read()
            receipt.update(http_status=r.status, content_type=r.headers.get("Content-Type"))
            if b"QUERY_STATUS" in body or b"<html" in body.lower():
                raise ValueError("Service returned non-CSV response: " + body[:300].decode(errors="replace"))
            receipt.update(bytes=len(body), sha256=hashlib.sha256(body).hexdigest())
            (OUT / (name + ".csv")).write_bytes(body)
    except Exception as err:
        receipt["error"] = str(err)
    print(name,receipt.get("bytes"),receipt.get("error"))
    return name, receipt


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    with ThreadPoolExecutor(max_workers=3) as pool:
        receipts = dict(pool.map(get_one, SERVICES))
    (OUT / "query.json").write_text(json.dumps(receipts, indent=2) + "\n")
