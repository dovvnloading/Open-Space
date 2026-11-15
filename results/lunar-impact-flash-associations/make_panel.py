"""Rebuild the panel and measurements from aligned, 1.31-m/pixel 8-bit display crops."""
import csv
from pathlib import Path
import numpy as np
from PIL import Image, ImageDraw, ImageFont
ROOT=Path(__file__).resolve().parent
def main():
    archive=np.load(ROOT/"registered_crops.npz")
    with (ROOT/"frame_sources.csv").open() as f:
        rows=list(csv.DictReader(f))
    font=ImageFont.load_default(size=17)
    title=ImageFont.load_default(size=24)
    panel=Image.new("RGB",(840,350),"white")
    draw=ImageDraw.Draw(panel)
    draw.text((15,12),"One new crater near 29.2322° N, 90.5225° E",font=title,fill="black")
    for i,row in enumerate(rows):
        arr=archive[row["product"]]
        a,b=np.percentile(arr,[2,98])
        im=Image.fromarray(np.clip((arr-a)*255/max(1,b-a),0,255).astype("uint8")).resize((190,190),Image.Resampling.NEAREST)
        x=15+i*207
        panel.paste(im,(x,62))
        draw.text((x,265),row["date"],font=font,fill="black")
        draw.text((x,290),row["product"],font=font,fill="black")
        draw.text((x,318),row["classification"],font=font,fill="black")
    panel.save(ROOT/"crater_time_series.png")
if __name__=="__main__":
    main()
