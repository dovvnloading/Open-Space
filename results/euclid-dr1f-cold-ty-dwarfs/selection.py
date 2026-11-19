#!/usr/bin/env python3
"""Independent T/Y search on every Euclid DR1-F wide HATS Parquet row group.

Reuses ONLY the publisher-index, S3 reader, row-group and checkpoint machinery
from euclid-dr1-wide-highz-vetting. The science mask and controls are independent.
"""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/"output"
PARENT=ROOT.parent/"euclid-dr1-wide-highz-vetting"
import sys
sys.path.insert(0,str(PARENT))
import scan as wide

BATCH_SIZE=24
EPOCH=2024.25
SELECTION_VERSION="cold-TY-NISP-2026-11-19-v1"
# Absolute photometric colours are AB flux ratios; all fluxes from consistent
# 3-FWHM aperture columns, never mixed with Vega CatWISE magnitudes.
def asfloat(tab,name):
    return wide.as_float(tab,name)

def cold_mask(tab):
    j,je=asfloat(tab,"flux_j_3fwhm_aper"),asfloat(tab,"fluxerr_j_3fwhm_aper")
    h,he=asfloat(tab,"flux_h_3fwhm_aper"),asfloat(tab,"fluxerr_h_3fwhm_aper")
    y,ye=asfloat(tab,"flux_y_3fwhm_aper"),asfloat(tab,"fluxerr_y_3fwhm_aper")
    v,ve=asfloat(tab,"flux_vis_3fwhm_aper"),asfloat(tab,"fluxerr_vis_3fwhm_aper")
    sp=asfloat(tab,"spurious_prob")
    ell=asfloat(tab,"ellipticity")
    det=asfloat(tab,"vis_det")
    fj=asfloat(tab,"flag_j")
    fh=asfloat(tab,"flag_h")
    with np.errstate(divide="ignore",invalid="ignore",over="ignore"):
        js=j/je
        vs=np.abs(v/ve)
        ys=y/ye
        hs=h/he
        yj=-2.5*np.log10(y/j)
        jh=-2.5*np.log10(j/h)
        # A low-SNR Y source gives a *lower* limit on Y-J.
        y_lower=-2.5*np.log10((2*ye)/j)
        yj_test=np.where(ys>=2,yj,y_lower)
        # A low-SNR H source gives an *upper* limit on J-H.
        jh_upper=-2.5*np.log10(j/(2*he))
        jh_test=np.where(hs>=2,jh,jh_upper)
        base=((det==0)&(j>0)&(je>0)&(he>0)&(ye>0)&
              (ve>0)&(js>=5)&(vs<3)&
              np.isfinite(ell)&(ell<.35)&
              (fj==0)&(fh==0)&
              ((~np.isfinite(sp))|(sp<.7)))
        strict=(base&(js>=6)&(vs<2)&(ell<.25)&
                ((~np.isfinite(sp))|(sp<.5))&
                np.isfinite(yj_test)&(yj_test>=.35)&
                np.isfinite(jh_test)&(jh_test<=.15))
    return base,strict,yj,jh,y_lower

def cold_row(tab,i,yj,jh,yl):
    r=wide.row_dict_original(tab,i,yj,jh,yl)
    def num(name):
        try:return float(r[name])
        except (KeyError,TypeError,ValueError):return math.nan
    j=num("flux_j_3fwhm_aper")
    h=num("flux_h_3fwhm_aper")
    he=num("fluxerr_h_3fwhm_aper")
    y=num("flux_y_3fwhm_aper")
    ye=num("fluxerr_y_3fwhm_aper")
    if j>0 and he>0:
        low_h=not(h/he>=2)
        r["J_minus_H_for_selection"]=(-2.5*math.log10(j/(2*he)) if low_h
                                         else r["J_minus_H_AB"])
        r["J_minus_H_is_upper_limit"]=low_h
    else:
        r["J_minus_H_for_selection"]=None
        r["J_minus_H_is_upper_limit"]=True
    if j>0 and ye>0:
        low_y=not(y/ye>=2)
        r["Y_minus_J_for_selection"]=(-2.5*math.log10((2*ye)/j) if low_y
                                         else r["Y_minus_J_AB"])
        r["Y_minus_J_is_lower_limit"]=low_y
    else:
        r["Y_minus_J_for_selection"]=None
        r["Y_minus_J_is_lower_limit"]=True
    r["Euclid_epoch_assumed_year"]=EPOCH
    r["point_source_evidence"]="ELLIPTICITY_PROXY_ONLY_NOT_A_NISP_PSF_CLASSIFICATION"
    return r

wide.row_dict_original=wide.row_dict
wide.row_dict=cold_row
wide.mask_rows=cold_mask
wide.ROOT=ROOT
wide.OUT=OUT
wide.BATCHES=OUT/"batches"
wide.BATCH_SIZE=BATCH_SIZE
wide.load_controls=lambda: json.loads((OUT/"control_catalog.json").read_text())

def main():
    p=Path(__file__)
    digest=hashlib.sha256(p.read_bytes()).hexdigest()
    out=OUT/"batches"
    if out.exists():
        for mpath in sorted(out.glob("batch_*_manifest.json")):
            man=json.loads(mpath.read_text())
            if man.get("selection_sha256")!=digest or man.get("selection_version")!=SELECTION_VERSION:
                raise RuntimeError("Checkpoint selection code changed: requires NEW result root, never silent replay")
    wide.main()
    for mpath in sorted(out.glob("batch_*_manifest.json")):
        man=json.loads(mpath.read_text())
        if "selection_sha256" not in man:
            man["selection_sha256"]=digest
            man["selection_version"]=SELECTION_VERSION
            mpath.write_text(json.dumps(man,indent=2,sort_keys=True)+"\n")

if __name__=="__main__":
    main()
