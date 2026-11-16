"""Collapse the two nearly simultaneous bands; test sky- and CCD-fixed offsets."""
import csv
from itertools import product
import numpy as np
from cutouts import ROOT, read_cutout, rows

def constant_chi(values, cov):
    weights=np.linalg.inv(cov)
    variance=np.linalg.inv(weights.sum(axis=0))
    answers=[]
    for signs in product([1,-1],repeat=len(values)-1):
        xy=values*np.array([1,*signs])[:,None]
        center=variance@np.einsum("ijk,ik->j",weights,xy)
        residual=xy-center
        answers.append(np.einsum("ij,ijk,ik->",residual,weights,residual))
    return min(answers)

def main():
    with (ROOT/"output"/"fits.csv").open() as stream:
        fits={r["id"]:r for r in csv.DictReader(stream) if r["kind"]=="Science" and r["variant"]=="nominal"}
    records=[]
    for frame in ["sky","detector"]:
        groups={}
        for t in rows("target.csv"):
            r=fits[t["measurement_id"]]
            xy=np.array([float(r["east_arcsec"]),float(r["north_arcsec"])])
            cov=np.array([[float(r["sigma_east"])**2,float(r["cov_east_north"])],
                          [float(r["cov_east_north"]),float(r["sigma_north"])**2]])
            if frame=="detector":
                h=read_cutout(t["measurement_id"],"Science")[0][0]
                mat=np.array([[float(h[f"PC{i}_{j}"])*3600*float(h[f"CDELT{i}"])
                               for j in (1,2)] for i in (1,2)])
                # Retain approximate arcsec units for the same floor in either frame.
                rotation=np.sqrt(abs(np.linalg.det(mat)))*np.linalg.inv(mat)
                xy=rotation@xy
                cov=rotation@cov@rotation.T
            groups.setdefault(t["visit"][:8],[]).append((xy,cov))
        xy=np.array([np.mean([v for v,c in rs],axis=0) for rs in groups.values()])
        # Only the formal random component is averaged down within night; the
        # additional floor below is applied once per night with full correlation.
        cov=np.array([sum(c for v,c in rs)/len(rs)**2 for rs in groups.values()])
        for floor in [0.,.03,.06,.10]:
            records.append(dict(frame=frame,night_floor_arcsec=floor,chi2=constant_chi(xy,cov+floor**2*np.eye(2)),dof=4))
    with (ROOT/"output"/"epoch_checks.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(records[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(records)
if __name__=="__main__":
    main()
