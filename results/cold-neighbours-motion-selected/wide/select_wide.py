"""Apply the pilot cuts in brighter, independent fields."""
import csv
import sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/"pilot"))
from select_catwise import select,load_corrections
def main():
    offsets=load_corrections();rows=[]
    for file in sorted((ROOT/"inputs").glob("cap*.csv")):
        for r in csv.DictReader(file.open()):
            selected=select(r,offsets)
            if selected and selected["bright_tier"]:
                rows.append(dict(field=file.stem,**selected))
    if rows:
        with (ROOT/"selected.csv").open("w") as f:
            w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
    print(len(rows))
if __name__=="__main__":main()
