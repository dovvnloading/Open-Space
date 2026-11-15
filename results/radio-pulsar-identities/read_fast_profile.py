"""Read the vector trace in the published Liu et al. (2024) Fig. 1, left panel."""
from pathlib import Path
import csv
import hashlib

import fitz

root=Path(__file__).parent
pdf=root/'fast2024.pdf'
assert hashlib.sha256(pdf.read_bytes()).hexdigest()=='fe36ae9d2bda24f18a04e86bbc6c40f05cc0da0f3f358eefbd180a9559fe913e'
with fitz.open(pdf) as doc:
    draws=doc[2].get_drawings()
# Rendered x ticks at 220 and 230 deg, rendered y ticks at 0 and 1 normalized intensity.
x220=draws[39]['rect'].x0
x230=draws[43]['rect'].x0
y0=draws[2]['rect'].y0
y1=draws[10]['rect'].y0
frame=draws[1]['rect']
points=[seg[2] for seg in draws[33]['items'] if seg[0]=='l' and frame.x0<=seg[2].x<=frame.x1]
with (root/'liu2024_fig1_j1611_profile.csv').open('w') as out:
    w=csv.writer(out)
    w.writerow(['pulse_longitude_deg','normalized_intensity'])
    for p in points:
        w.writerow([round(220+(p.x-x220)*10/(x230-x220),3),round((p.y-y0)/(y1-y0),5)])
