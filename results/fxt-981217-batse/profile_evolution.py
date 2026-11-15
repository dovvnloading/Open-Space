"""Compare the two count profiles without adopting a common spectral response."""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
SOD = 86400
T0 = 51164.079960


def fluence_fraction(counts, covariance, intervals):
    """Propagate shared baseline-fit errors through fractional, summed counts."""
    f = counts[intervals].sum()/counts.sum()
    derivative = (intervals.astype(float)-f)/counts.sum()
    return f, np.sqrt(derivative @ covariance @ derivative)


def compare():
    w = np.genfromtxt(ROOT/"inputs/wfc_lc.csv", delimiter=",", names=True)
    a = np.genfromtxt(ROOT/"inputs/batse_lc.csv", delimiter=",", names=True)
    t = (a["mjd"]-T0)*SOD
    begin, end = (w["mjd_begin"]-T0)*SOD, (w["mjd_end"]-T0)*SOD
    mid = (begin+end)/2
    src = (mid >= 0) & (mid < 111)
    offw = ((mid>=-210)&(mid<-30))|((mid>=141)&(mid<321))
    weights = 1/w["flux_error"][offw]**2
    wbase = np.average(w["flux"][offw], weights=weights)
    duration = (end-begin)[src]
    xcounts = (w["flux"][src]-wbase)*duration
    # Corrected WFC flux and per-bin archive sigma; the same baseline is subtracted from all bins.
    Vw = np.outer(duration,duration)*(np.diag(w["flux_error"][src]**2) + 1/weights.sum())
    duration_mid = (end[src][-1]+begin[src][0])/2
    early = mid[src] <= duration_mid
    wf, we = fluence_fraction(xcounts,Vw,early)
    out = {}
    for order, edge in [(1,210), (2,210), (1,150), (2,150)]:
        off = ((t>=-edge)&(t<-30))|((t>=141)&(t<111+edge))
        X = np.polynomial.polynomial.polyvander((t-55.5)/100,order)
        allcounts = sum(a[f"d{idx}_ch{ch}"] for idx in range(8) for ch in [2,3])
        residual = allcounts-X@np.linalg.lstsq(X[off],allcounts[off],rcond=None)[0]
        gamma = []
        design = []
        var = []
        for lo,hi,dt in zip(begin[src],end[src],duration):
            use = (t>=lo)&(t<hi)
            assert use.sum()>=8
            # Average over covered 1.024-second BATSE bins, then integrate on the WFC grid.
            gamma.append(residual[use].mean()*dt/1.024)
            design.append(X[use].mean(0)*dt/1.024)
            var.append((dt/1.024)**2/use.sum())
        design = np.array(design)
        s2 = np.sum(residual[off]**2)/(off.sum()-order-1)
        Vb = s2*(np.diag(var) + design@np.linalg.solve(X[off].T@X[off],design.T))
        bf, be = fluence_fraction(np.array(gamma), Vb, early)
        out[f"degree{order}_outer{edge}"] = {"hard_fraction":round(float(bf),5),"hard_se":round(float(be),5),
                                            "hard_minus_soft":round(float(bf-wf),5),
                                            "difference_se":round(float(np.hypot(be,we)),5)}
    return {"wfc_bins":int(src.sum()), "early_bins":int(early.sum()),
            "last_early_edge_seconds":round(float(end[src][early][-1]),3),
            "last_source_edge_seconds":round(float(end[src][-1]),3),
            "wfc_fraction":round(float(wf),5),"wfc_se":round(float(we),5),"batse":out}


if __name__ == "__main__":
    print(json.dumps(compare(),indent=2))
