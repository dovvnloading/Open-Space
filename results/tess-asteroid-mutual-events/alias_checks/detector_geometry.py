"""Test phase-locked dimmings against detector sampling and local sky background.

Place measure.py on the Python module path. Run with the source release and the result directory
that holds measurements.csv. Positive depth is fainter in the predicted window.
"""
import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.stats import rankdata, ks_2samp
from measure import residuals, event_rows


def pixel_terms(x, y):
    # All integer lattice directions through second order, plus the third
    # harmonics of the individual pixel axes. Uniformly spaced cells are 1 px.
    terms = []
    for i, j in [(0,1),(0,2),(0,3),(1,-2),(1,-1),(1,0),(1,1),
                 (1,2),(2,-2),(2,-1),(2,0),(2,1),(2,2),(3,0)]:
        p = 2*np.pi*(i*x + j*y)
        terms.extend((np.cos(p), np.sin(p)))
    return np.column_stack(terms)


def background_terms(a):
    # Robust rank normalization keeps occasional bright field stars bounded.
    bg = rankdata(a[:,5])/(len(a)+1) - .5
    rms = rankdata(a[:,6])/(len(a)+1) - .5
    return np.column_stack((bg,rms,bg*bg,rms*rms,bg*rms))


def crossfit(features, values, cycle, away):
    # Each prediction is learned on six other groups of whole cycles, and
    # avoids the entire central half-cycle, including event ingress/egress.
    design = np.column_stack((np.ones(len(values)), features))
    prediction = np.empty_like(values)
    for fold in range(7):
        test = cycle % 7 == fold
        train = (~test) & away
        coefs = np.linalg.lstsq(design[train],values[train],rcond=None)[0]
        prediction[test] = design[test] @ coefs
    return prediction


def bearing(a, period, center, phase_width, residual):
    jd = a[:,0]
    ph = (jd-center)*24/period
    cycles = np.rint(ph).astype(int)
    distance = np.abs(ph-cycles)
    inside, away = distance < phase_width/2, distance > .27
    def stats(v):
        ev = event_rows(jd-center, v, period, 0., phase_width)
        return np.median(ev[:,3]), int((ev[:,3]>0).sum()), len(ev)
    out = dict()
    obs, npos, nevents = stats(residual)
    out.update(original_depth_mag=obs,n_positive=npos,n_events=nevents)
    pixel = pixel_terms(a[:,1],a[:,2])
    bg = background_terms(a)
    for label,X in [("pixel",pixel),("background",bg),("pixel_background",np.column_stack((pixel,bg)))]:
        prediction = crossfit(X, residual, cycles, away)
        depth,pos,total = stats(residual - prediction)
        out.update({label+"_corrected_depth_mag": depth,label+"_corrected_n_positive":pos,
                    label+"_predicted_depth_mag":stats(prediction)[0]})
    # Investigate two independently moving coordinates, and ranks of sky
    # background/RMS. The quoted tests are descriptive because visits cluster.
    for i, col in [(1,'x'),(2,'y'),(4,'formal_error'),(5,'diff_background'),(6,'bg_rms')]:
        raw = a[:,i]
        values = raw % 1 if i in (1,2) else raw
        ranks = rankdata(values)/(len(values)+1)
        statistic = ks_2samp(values[inside],values[away]).statistic
        out[col+"_in_median"] = np.median(values[inside])
        out[col+"_out_median"] = np.median(values[away])
        out[col+"_ks_distance"] = statistic
        out[col+"_in_minus_out_median_rank"] = np.median(ranks[inside])-np.median(ranks[away])
    for keep,label in [(a[:,6]<np.quantile(a[:,6],.75),"exclude_high_rms"),
                       (np.abs(a[:,5]-np.median(a[:,5])) <
                        np.quantile(abs(a[:,5]-np.median(a[:,5])),.75),"exclude_extreme_background")]:
        ev = event_rows(jd[keep]-center,residual[keep],period,0.,phase_width)
        out[label+"_depth_mag"]=np.median(ev[:,3])
        out[label+"_positive"]=int(sum(ev[:,3]>0))
        out[label+"_events"]=len(ev)
    # A near-integer count of native exposures per cycle could repeatedly
    # sample the same time stamp offset in the event profile. Count occupancy
    # of eighths of a native cadence among the predicted event centers.
    cadence = np.median(np.diff(jd)[np.diff(jd)<.023])  # days
    center_ages = (center+np.unique(cycles)*period/24 - jd[0])/cadence
    occ = np.histogram(center_ages % 1,bins=np.linspace(0,1,9))[0]
    out["native_cadence_s"]=cadence*86400
    out["cadences_per_event_spacing"]=period/(24*cadence)
    out["cadence_phase_bins_occupied_of_8"]=int((occ>0).sum())
    return out


def main(release, output, destination):
    meta={r[0]:float(r[3]) for l in (release/"release.merge").read_text().splitlines() if (r:=l.split())}
    with (output/'measurements.csv').open() as f:
        measurements={int(r['mpc']):r for r in csv.DictReader(f)}
    rows=[]
    for num in [29489,3177]:
        r=measurements[num]
        a,y,_=residuals(release,num,meta[str(num)])
        per=float(r['recurrence_h'])
        diag=bearing(a,per,float(r['first_center_jd']),float(r['box_width_h'])/per,y)
        rows.append(dict(mpc=num,n_samples=len(a),**diag))
    destination.mkdir(parents=True,exist_ok=True)
    with (destination/'detector_geometry.csv').open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader()
        w.writerows([{k:round(v,6) if isinstance(v,(float,np.floating)) else v for k,v in r.items()} for r in rows])


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("release",type=Path)
    p.add_argument("results",type=Path)
    p.add_argument("--out",type=Path,default=Path(__file__).resolve().parent)
    a=p.parse_args()
    main(a.release,a.results,a.out)
