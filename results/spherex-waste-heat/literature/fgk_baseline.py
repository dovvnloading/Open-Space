"""Extract published Gaia/WISE photometry from the authors' hash-pinned table."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import urllib.request

ROOT = Path(__file__).resolve().parent
REV = "fc257c9c3fd54f4edb7e94d61cc5b9a20b27bc2b"
URL = f"https://raw.githubusercontent.com/contardog/NotATechnosignatureSearch/{REV}/53candidates_allcol.csv"
SHA = "63efdd5df2f1b1a4a9d62ded586fa866a73f18b830f1f2cc6c437ca74bcf3c34"
PROPS = {
    "source_id": "Gaia DR3 source ID", "designation_2mass": "2MASS ID", "designation_allwise": "DESIGNATION",
    "unwise_objid": "unwise_objid", "allwise_offset_arcsec": "angular_distance",
    "gaia_g_mag": "phot_g_mean_mag", "gaia_bp_rp": "bp_rp", "gaia_ruwe": "ruwe",
    "allwise_ph_qual": "ph_qual", "allwise_cc_flags": "cc_flags", "allwise_ext_flg": "ext_flag",
    "allwise_w1": "w1mpro", "allwise_w1_error": "w1mpro_error",
    "allwise_w2": "w2mpro", "allwise_w2_error": "w2mpro_error",
    "allwise_w3_or_limit": "w3mpro", "allwise_w3_error": "w3mpro_error",
    "allwise_w4_or_limit": "w4mpro", "allwise_w4_error": "w4mpro_error",
    "unwise_w1_native": "mag_w1_vg", "unwise_w2_native": "mag_w2_vg",
    "unwise_w1_fracflux": "fracflux_w1", "unwise_w2_fracflux": "fracflux_w2",
    "pred_w1_native": "hatW1_pred", "pred_w2_native": "hatW2_pred"
}

def table(path, rows):
    with path.open("w") as stream:
        wr = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        wr.writeheader()
        wr.writerows(rows)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input", type=Path)
    a = p.parse_args()
    data = a.input.read_bytes() if a.input else urllib.request.urlopen(URL, timeout=120).read()
    digest = hashlib.sha256(data).hexdigest()
    if digest != SHA:
        raise ValueError(f"Unexpected bytes SHA256 {digest}")
    inp = list(csv.DictReader(io.StringIO(data.decode())))
    by_id = {r["source_id"]: r for r in csv.DictReader((ROOT / "contardo_hogg.csv").open())}
    if len(inp) != 53 or len(by_id) != 53 or {r["Gaia DR3 source ID"] for r in inp} != set(by_id):
        raise ValueError("Published target sets disagree")
    astrometry = []
    photometry = []
    for r in inp:
        sid = r["Gaia DR3 source ID"]
        astrometry.append(dict(source_id=sid, ra=r["ra_gaia"], dec=r["dec_gaia"], ref_epoch=2016.0,
                               pmra=r["pmra"], pmdec=r["pmdec"], parallax=r["parallax"], ruwe=r["ruwe"]))
        diag = {k: r[v] for k, v in PROPS.items()}
        # The author CSV serialized 0000 as the number 0.
        diag["allwise_cc_flags"] = diag["allwise_cc_flags"].zfill(4)
        for b in (3, 4):
            q = r["ph_qual"][b-1]
            diag[f"w{b}_allwise_kind"] = {"A":"SNR>=10", "B":"3<=SNR<10", "C":"2<=SNR<3", "U":"upper_limit", "X":"no_measurement", "Z":"no_quality_estimate"}.get(q, "other")
        photometry.append(diag)
    table(ROOT / "fgk_astrometry.csv", astrometry)
    table(ROOT / "fgk_baseline.csv", photometry)
    (ROOT / "fgk_source_receipt.json").write_text(json.dumps(dict(
        url=URL, commit=REV, git_blob="e5fd4b6c538c2ab82705c0179e3bb4f6277e21bd",
        bytes=len(data), sha256=digest, rows=len(inp),
        source="Author machine-readable supplement to Contardo & Hogg 2024; includes original crossmatches and trained-model output."), indent=2)+"\n")

if __name__ == "__main__":
    main()
