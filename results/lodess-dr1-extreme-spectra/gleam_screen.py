"""Screen covered LoDeSS positions with simultaneous-in-survey GLEAM subband photometry."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def run(fields_file, matches, output):
    src = pd.read_csv(fields_file, index_col="row")
    for survey in ("gleamegcat","lotss","vlssr","gmrtas150m","nvss"):
        src = src.join(pd.read_csv(matches / f"match_{survey}.csv.gz", index_col="row"))
    # Limit the largest GLEAM beams with a 3-arcmin LoTSS neighbour-flux screen.
    selected = ((src.quality!="poor") & (src.separation_deg<4) & (src.gleamegcat_n_match==1)
                & (src.gleamegcat_sep_as<30) & (src.lotss_sep_as<15) &
                (src.lotss_flux_jy/src.lotss_sum180_jy>.9) & (src.lotss_dc_maj<15) &
                (src.gmrtas150m_sep_as<20) & (src.nvss_sep_as<20))
    bands = pd.read_csv(matches/"gleam_subbands.csv.gz")
    rows = []
    for idx, g in bands[bands.row.isin(src.index[selected])].groupby("row"):
        groups = []
        # Average in four broad groups. Differential calibration can be correlated.
        # Take the largest uncertainty in each group instead of shrinking by sqrt(N).
        for lo, hi in [(75,100),(100,140),(140,178),(178,228)]:
            b = g[(g.freq_mhz>=lo)&(g.freq_mhz<hi)].dropna(subset=["flux_jy","fit_err_jy","rms_jy_beam"])
            if len(b)<2:
                break
            sig = np.maximum(b.fit_err_jy,b.rms_jy_beam)
            flux = b.flux_jy.mean()
            groups.append((np.exp(np.log(b.freq_mhz).mean()),flux,np.hypot(sig.max(),.05*abs(flux))))
        if len(groups)!=4:
            continue
        nu, flux, err = np.array(groups).T
        if np.any(flux<=2*err):
            continue
        # External scale is a common normalization in this slope-only test.
        design = np.column_stack((np.ones(4),np.log(nu/150)))
        w = (flux/err)**2
        cov = np.linalg.inv((design.T*w)@design)
        par = cov@(design.T@(w*np.log(flux)))
        chi = sum(w*(np.log(flux)-design@par)**2)
        se = np.sqrt(cov[1,1]*max(1,chi/2))
        s = src.loc[idx]
        rows.append([idx,s.field,s.quality,s.gleamegcat_name,s.lotss_name,*flux,par[1],se,
                     s.vlssr_flux_jy,s.lotss_flux_jy,s.gmrtas150m_flux_jy,s.nvss_flux_jy])
    out = pd.DataFrame(rows,columns=["row","field","quality","gleam","lotss","group1_jy","group2_jy",
                                    "group3_jy","group4_jy","slope","slope_err","vlssr_jy","lotss_jy","tgss_jy","nvss_jy"])
    out.to_csv(output/"gleam_screen.csv",index=False,float_format="%.8g")
    print(f"{int(selected.sum())} eligible counterparts; {len(out)} fit; {(out.slope>3*out.slope_err).sum()} rising above 3 formal errors")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    for arg in ["fields_file","matches","output"]:
        parser.add_argument(arg,type=Path)
    a=parser.parse_args()
    run(a.fields_file,a.matches,a.output)
