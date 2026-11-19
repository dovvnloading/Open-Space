#!/usr/bin/env python3
"""Frozen, control-tuned v1 eligibility selection. NOT a calibrated quasar classifier."""
import numpy as np
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"euclid-dr1-wide-highz-vetting"))
import scan

TOP_PER_PARTITION=200

def arrays(tab):
    g=lambda k:scan.as_float(tab,k)
    v,ve=g("flux_vis_3fwhm_aper"),g("fluxerr_vis_3fwhm_aper")
    y,ye=g("flux_y_3fwhm_aper"),g("fluxerr_y_3fwhm_aper")
    j,je=g("flux_j_3fwhm_aper"),g("fluxerr_j_3fwhm_aper")
    h,he=g("flux_h_3fwhm_aper"),g("fluxerr_h_3fwhm_aper")
    ell,mu,spur=g("ellipticity"),g("mumax_minus_mag"),g("spurious_prob")
    with np.errstate(divide="ignore",invalid="ignore",over="ignore"):
        vsn=v/ve;jsn=j/je;hsn=h/he
        yj=-2.5*np.log10(y/j)
        yl=-2.5*np.log10(2*ye/j)
        col=np.where(y>0,yj,yl)
        jh=-2.5*np.log10(j/h)
        # This is an inclusion/triage screen, not a strict VIS "dropout".
        eligible=((j>.65)&(h>0)&(je>0)&(he>0)&(ve>0)&
                  (jsn>=6)&(hsn>=4)&(np.abs(vsn)<7)&
                  (col>=0)&(jh>=-.2)&(jh<=.48)&
                  (ell<.4)&(mu>=-1.45)&(mu<=-.75)&
                  (~np.isfinite(spur)|(spur<.7)))
        optical_max=np.full(tab.num_rows,-np.inf)
        for b in scan.EXT_BANDS:
            f,e=g("flux_"+b+"_3fwhm_aper"),g("fluxerr_"+b+"_3fwhm_aper")
            osn=np.where(e>0,f/e,-np.inf)
            optical_max=np.maximum(optical_max,osn)
        fj,fh=g("flag_j"),g("flag_h")
        flags_bad=(~np.isfinite(fj))|(~np.isfinite(fh))|(fj!=0)|(fh!=0)
        opt_gt3=optical_max>3
        # Fixed, auditable priority heuristic. NOT a posterior.
        score=(3*np.clip(col,0,2)+np.clip(jsn,0,25)/7
               -2*np.abs(jh-.12)-np.maximum(np.abs(vsn)-1,0)/4
               -1.5*flags_bad-1.5*opt_gt3)
        score=np.where(np.isfinite(score),score,-1e9)
    return dict(eligible=eligible, score=score, VSNR=vsn, JSNR=jsn,
                HSNR=hsn, YJ=yj, YJ_LIMIT=yl, YJ_USED=col,
                JH=jh, bad_flags=flags_bad,optical_gt3=opt_gt3,
                optical_max=optical_max)
