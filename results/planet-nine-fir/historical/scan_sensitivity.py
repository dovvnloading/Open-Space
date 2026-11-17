"""Re-run the all-track temporal test across nine in-scan template windows."""
import argparse
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orbit-dir",type=Path,required=True)
    args = parser.parse_args()
    rows=[]
    with tempfile.TemporaryDirectory() as tmp:
        tmpfile=Path(tmp)/"result.json"
        for width in (1.2,1.5,1.8):
            for width_multiplier in (3.,6.,9.):
                subprocess.run([sys.executable,str(ROOT/"all_scan_tracks.py"),
                    "--orbit-dir",str(args.orbit_dir),"--width",str(width),
                    "--window",str(width_multiplier),"--output",str(tmpfile)],check=True)
                data=json.loads(tmpfile.read_text())
                rows.append(dict(fwhm_arcmin=width,window_halfwidth_arcmin=width*width_multiplier,
                    max_of_minimum_July_and_late_Jy=data["largest_minimum_flux_Jy"],
                    max_late_at_July_ge_0p2_Jy=data["max_late_Jy_at_July_at_least_0p20"],
                    fixed_July_Jy=data["fixed_July_late_Jy"][0],fixed_late_Jy=data["fixed_July_late_Jy"][1]))
    with (ROOT/"scanpi"/"filter_sensitivity.csv").open("w",newline="") as f:
        w=csv.DictWriter(f,rows[0]);w.writeheader();w.writerows(rows)


if __name__=="__main__":
    main()
