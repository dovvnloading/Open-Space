"""Reproduce the AstroVink 33-object A/B inspection sheet from the authors' published JPEGs."""
import argparse
import csv
import hashlib
import io
import zipfile
from pathlib import Path
from urllib.request import urlopen

from PIL import Image, ImageDraw
from astrovink import ARCHIVE, fetch

ZIP_URL = "https://zenodo.org/records/17425610/files/q1_missed_targets_AstroVink.zip?download=1"


def gallery(directory, output):
    rows = fetch(directory, "astrovink.csv", ARCHIVE)
    archive = fetch(directory, "astrovink.zip", ZIP_URL)
    if hashlib.md5(archive.read_bytes()).hexdigest() != "d606ea0e49f269d2e87530ae5cb2ba83":
        raise ValueError("Unexpected AstroVink image archive")
    with rows.open() as stream:
        candidates=[r for r in csv.DictReader(stream) if r["grade"] in {"A","B"}]
    w,h=166,176
    gallery=Image.new("RGB",(6*w,6*h),"#15191d")
    labels=ImageDraw.Draw(gallery)
    with zipfile.ZipFile(archive) as z:
        for i,r in enumerate(candidates):
            key=r["id_str"].replace("_-","_NEG")
            prefix="q1_missed_targets_AstroVink/mtf_vis_y_j/"+key+"_"
            members=[p for p in z.namelist() if p.startswith(prefix) and p.endswith(".jpg")]
            if len(members)!=1:
                raise ValueError("Expected one stretch for "+key)
            im=Image.open(io.BytesIO(z.read(members[0]))).convert("RGB")
            im.thumbnail((150,150))
            x,y=i%6*w,i//6*h
            gallery.paste(im,(x+5,y))
            labels.text((x+3,y+153),r["grade"]+"  "+r["name"][4:],fill="white")
    gallery.save(output/"astrovink_gallery.jpg",quality=88)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_dir",type=Path)
    p.add_argument("--output",type=Path,default=Path(__file__).resolve().parent)
    args=p.parse_args()
    gallery(args.input_dir,args.output)
