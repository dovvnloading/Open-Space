"""Deduplicate white dwarf spectral products before applying one Gaia selection."""
import argparse
import csv
import gzip
import hashlib
import json
import sys
import struct
import zlib
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from zero_point import zpt

HERE = Path(__file__).resolve().parent
RESEARCH = HERE.parent
sys.path.insert(0, str(RESEARCH/"desi"))
from fetch import BATCH, COLUMNS, TAP, digest, fetch_block

SURVEYS = ["DESI", "SDSS_DR19", "SDSS_DR20", "LAMOST"]


def verified_chunks(path, receipts):
    out = []
    for b in receipts:
        p = path/f"gaia_{b['batch']:03}.csv"
        if digest(p.read_bytes()) != b["response_sha256"]:
            raise ValueError(f"Gaia batch checksum mismatch for {p}")
        d = pd.read_csv(p,dtype={"source_id":str})
        if (len(d)!=b["rows"] or list(d.columns)!=COLUMNS or
                digest(",".join(d.source_id).encode())!=b["id_sha256"]):
            raise ValueError(f"Gaia batch identity mismatch for {p}")
        out.append(d)
    return pd.concat(out,ignore_index=True)


def sdss_ids(cache, rel, meta):
    # Need membership, including overlaps between DESI and each SDSS release.
    # Reuse memory-mapped, previously verified uncompressed FITS; access ONLY IDs
    # and labels. The gzipped release is the portable provenance reference.
    path = cache/f"{rel.lower()}_snow.fits"
    version = "0.6.0" if rel=="DR19" else "0.8.1"
    archive = cache/f"astraAllStarSnowWhite-{version}.fits.gz"
    with archive.open("rb") as f:
        if hashlib.file_digest(f,"sha256").hexdigest()!=meta["sha256"]:
            raise ValueError(f"{rel} public FITS archive hash mismatch")
    # The source gzip is hash-verified above. Check that its cached decompressed
    # representation, used by upstream sdss/screen.py, matches its CRC and length.
    with archive.open("rb") as f:
        f.seek(-8,2)
        expected_crc, expected_length = struct.unpack("<II",f.read(8))
    crc=0
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(8*1024**2),b""):
            crc=zlib.crc32(chunk,crc)
    if crc!=expected_crc or path.stat().st_size % 2**32 != expected_length:
        raise ValueError(f"{rel} decompressed FITS differs from its gzip")
    ids=set()
    with fits.open(path,memmap=True) as hdus:
        for h in hdus[1:]:
            r=h.data
            if r is None or not len(r):
                continue
            keep=(np.char.strip(r["classification"])!="") & (r["gaia_dr3_source_id"]>0)
            ids.update(str(i) for i in r["gaia_dr3_source_id"][keep])
    if len(ids)!=meta["unique_ids"]:
        raise ValueError(f"{rel} membership differs from upstream census")
    return ids


def main(cache, lamost, membership_out=None):
    desi_meta=json.loads((RESEARCH/"desi"/"receipt.json").read_text())
    sdss_meta=json.loads((RESEARCH/"sdss"/"summary.json").read_text())
    d = verified_chunks(cache,desi_meta["batches"])
    extra = verified_chunks(cache/"sdss_gaia",json.loads((RESEARCH/"sdss"/"gaia_receipt.json").read_text()))
    prior = pd.concat([d,extra],ignore_index=True)
    if prior.source_id.duplicated().any():
        raise ValueError("Expected the parent caches to be disjoint")
    try:
        with gzip.open(lamost,"rt") as f:
            lmo = {r["gaia_source_id"] for r in csv.DictReader(f) if int(r["gaia_source_id"])>0}
    except FileNotFoundError:
        raise SystemExit("Supply a copy of the public LAMOST archive verified by spectra/lamost_parent.py.")
    lmeta=json.loads((RESEARCH/"spectra"/"lamost_receipt.json").read_text())
    if digest(lamost.read_bytes())!=lmeta["archive_sha256"]:
        raise ValueError("LAMOST catalogue hash mismatch")
    missing=sorted(lmo-set(prior.source_id),key=int)
    dest=cache/"union_lamost_gaia"
    dest.mkdir(exist_ok=True)
    jobs=[(dest,i//BATCH, missing[i:i+BATCH]) for i in range(0,len(missing),BATCH)]
    with ThreadPoolExecutor(6) as pool:
        receipts=list(pool.map(fetch_block,jobs))
    newer=verified_chunks(dest,receipts)
    assert set(newer.source_id)==set(missing)
    all_data=pd.concat([prior,newer],ignore_index=True)
    memberships={"DESI":set(d.source_id), "LAMOST":lmo}
    for rel in ["DR19","DR20"]:
        memberships["SDSS_"+rel]=sdss_ids(cache,rel,sdss_meta["sources"][rel])
    all_ids=set(all_data.source_id)
    assert len(all_ids)==len(all_data) and all_ids==set.union(*memberships.values())
    assert len(memberships["SDSS_DR19"] & memberships["SDSS_DR20"])==sdss_meta["dr19_dr20_shared_ids"]
    assert len((memberships["SDSS_DR19"]|memberships["SDSS_DR20"]) & memberships["DESI"])==sdss_meta["matched_from_desi"]

    receipts_data=dict(tap=TAP,method="GET",utc=datetime.now(timezone.utc).isoformat(),
        query_builder="../desi/fetch.py::fetch_block",new_inputs="LAMOST IDs not present in DESI or SDSS",
        sorted_ids_sha256=digest(",".join(missing).encode()),additional=len(missing),batches=receipts)
    (HERE/"lamost_gaia_receipt.json").write_text(json.dumps(receipts_data,indent=2)+"\n")
    masks=np.zeros(len(all_data),dtype=int)
    for j, s in enumerate(SURVEYS):
        masks += all_data.source_id.isin(memberships[s]).to_numpy().astype(int)*(1<<j)
    all_data["survey_mask"]=masks
    all_data["surveys"]=["|".join(s for j,s in enumerate(SURVEYS) if m & (1<<j)) for m in masks]
    d=all_data
    valid=d.astrometric_params_solved.isin([31,95]) & d.parallax.notna() & (d.parallax_error>0) & d.pmra.notna() & d.pmdec.notna()
    quality=(d.ruwe.gt(0) & d.ruwe.lt(1.4) & d.visibility_periods_used.ge(9) & d.duplicated_source.eq(False)
             & d.ipd_frac_multi_peak.le(2) & d.ipd_gof_harmonic_amplitude.lt(0.1))
    zpt.load_tables()
    cols=["phot_g_mean_mag","nu_eff_used_in_astrometry","pseudocolour","ecl_lat","astrometric_params_solved"]
    d["zpt_mas"]=np.nan
    d.loc[valid,"zpt_mas"]=zpt.get_zpt(*(d.loc[valid,c].to_numpy() for c in cols),_warnings=False)
    corrected=d.parallax-d.zpt_mas
    d["corrected_snr"]=corrected/d.parallax_error
    d["vt_kms"]=4.74047*np.hypot(d.pmra,d.pmdec)/corrected.where(corrected.gt(0))
    d["quality"]=quality
    d["eligible_snr3"]=valid & quality & d.corrected_snr.ge(3)
    d["eligible_snr5"]=valid & quality & d.corrected_snr.ge(5)
    ctrl=pd.read_csv(RESEARCH/"control"/"controls.csv",dtype={"source_id":str})
    detail=["source_id","survey_mask","surveys","parallax","parallax_error","zpt_mas","corrected_snr",
            "pmra","pmdec","quality","vt_kms","eligible_snr3","eligible_snr5"]
    ctrl.merge(d[detail],on="source_id",how="left",validate="one_to_one").to_csv(HERE/"controls.csv",index=False)
    # Source membership and gate flags are reconstructible from the upstream receipts.
    # Hash membership separately for exact regeneration without duplicating 100k rows.
    membership_bytes = d.sort_values("source_id")[["source_id","survey_mask","eligible_snr3","eligible_snr5"]].to_csv(
        index=False,lineterminator="\n").encode()
    membership_sha = digest(membership_bytes)
    if membership_out is not None:
        membership_out.write_bytes(membership_bytes)
    overlap=[]
    for m, grp in d.groupby("survey_mask"):
        overlap.append(dict(mask=int(m),surveys=grp.surveys.iloc[0],unique=len(grp),
            eligible_snr3=int(grp.eligible_snr3.sum()),eligible_snr5=int(grp.eligible_snr5.sum())))
    pd.DataFrame(overlap).to_csv(HERE/"overlap.csv",index=False)
    summaries=[]
    tail=[]
    for cut in (3,5):
        for s in ["UNION"]+SURVEYS:
            sub=d[d[f"eligible_snr{cut}"] & (True if s=="UNION" else d.source_id.isin(memberships[s]))].sort_values("vt_kms",ascending=False)
            top=sub.iloc[0]
            found=sub[sub.vt_kms>1000]
            summaries.append(dict(sample=s,corrected_snr_min=cut,eligible=len(sub),above1000=len(found),above600=int((sub.vt_kms>600).sum()),
                               max_vt_kms=round(top.vt_kms,3),max_source_id=top.source_id,
                               null_binomial95_fraction=(float(-np.expm1(np.log(.05)/len(sub))) if not len(found) else None)))
            if s=="UNION":
                tail.append(sub.head(10)[detail].assign(threshold=cut))
    pd.concat(tail).to_csv(HERE/"top10.csv",index=False)
    summary=dict(unique=len(d),with_solution=int(valid.sum()),with_individual_zpt=int(d.zpt_mas.notna().sum()),
                 quality_and_calibrated=int((quality & d.zpt_mas.notna()).sum()),survey_bit_assignments=dict(zip(SURVEYS,[1,2,4,8])),
                 memberships={s:len(m) for s,m in memberships.items()},membership_and_eligibility_sha256=membership_sha,
                 results=summaries)
    (HERE/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    args=argparse.ArgumentParser()
    args.add_argument("cache",type=Path)
    args.add_argument("lamost",type=Path)
    args.add_argument("--membership-out",type=Path,help="Write exact per-source four-survey masks and both eligibility flags as CSV.")
    opts=args.parse_args()
    main(opts.cache,opts.lamost,opts.membership_out)
