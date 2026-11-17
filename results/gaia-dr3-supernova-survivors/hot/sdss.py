"""Use the source-level BOSS RV fields, independently of SnowWhite classification."""
import csv
import hashlib
import json
import pathlib
import sys

import numpy as np
from astropy.io import fits
from screen import HERE, catalogue

SIGMA95 = 1.64485362695
INPUTS = {"dr19":("0.6.0","a7489a2f7a0a155fe127badc7c21a5b21eea5791d097300e2f683e0ad6a923a9"),
          "dr20":("0.8.1","a19c7aa645459d83aabe58915693fa95312bd59b9c9559bd4f0459a57b819099")}


def run(cache):
    cat = {r["GaiaEDR3"]:r for r in catalogue() if r["GaiaEDR3"]}
    ids = np.array(sorted(map(int,cat)), dtype=np.int64)
    matches,meta = [],{}
    for rel,(ver,expected) in INPUTS.items():
        name=f"astraAllStarSnowWhite-{ver}.fits.gz"
        with (cache/name).open("rb") as f:
            assert hashlib.file_digest(f,"sha256").hexdigest()==expected
        # Uncompressed, memory-mapped FITS created from these bytes by ../sdss/screen.py.
        found = []
        with fits.open(cache/f"{rel}_snow.fits",memmap=True) as hdus:
            for h in hdus[1:]:
                if not h.header["NAXIS2"]: continue
                for r in h.data[np.isin(h.data["gaia_dr3_source_id"],ids)]:
                    ident=str(r["gaia_dr3_source_id"])
                    def clean(name):
                        v=r[name]
                        if isinstance(v,np.floating) and not np.isfinite(v):return None
                        if isinstance(v,np.generic):return v.item()
                        return v
                    row=dict(release=rel,source_id=ident,name=cat[ident]["Name"],
                             **{k:clean(k) for k in ("sdss_id","telescope","v_rad","e_v_rad","n_good_rvs",
                                               "std_v_rad","snr","classification","min_mjd","max_mjd")})
                    found.append(row)
        good = [r for r in found if r["v_rad"] is not None and r["e_v_rad"] is not None
                and r["e_v_rad"]>0 and r["n_good_rvs"]>0]
        tails=[r for r in good if abs(r["v_rad"])>=700]
        meta[rel]=dict(url=f"https://{rel}.sdss.org/sas/{rel}/spectro/astra/{ver}/summary/{name}",
                       compressed_sha256=expected,matched_rows=len(found),matched_unique_ids=len({r['source_id'] for r in found}),
                       usable_rv_unique_ids=len({r['source_id'] for r in good}),rows_above_700=tails)
        matches.extend(found)
    with (HERE/"sdss_match.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(matches[0]),lineterminator="\n");w.writeheader();w.writerows(matches)
    (HERE/"sdss_summary.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(json.dumps(meta,indent=2))

if __name__=="__main__": run(pathlib.Path(sys.argv[1]))
