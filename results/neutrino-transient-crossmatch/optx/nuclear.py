"""Project Table 1 positions and Table A.1 optical peak dates from arXiv:2511.19016.
Usage: python nuclear.py /path/to/arxiv-html
"""
import re
import sys
from pathlib import Path
from lxml import html
from optical import sex,write


def cells(row):
    return [" ".join(c.text_content().split()) for c in row.xpath("td|th")]


def main():
    x=html.fromstring(Path(sys.argv[1]).read_bytes())
    names={};group=""
    for row in x.xpath('//*[@id="S3.T1.2.1"]//tr'):
        a=cells(row)
        if len(a)<3:continue
        if not a[1]:group=a[0];continue
        if re.match(r"\d\d:\d\d:",a[1]):
            names[a[0].split()[0]]=(sex(a[1],True),sex(a[2]),group)
    rows=[]
    for row in x.xpath('//*[@id="A1.T3.2"]//tr'):
        a=cells(row)
        if len(a)<2 or not a[1].startswith(("2024","2025","2026")):continue
        ra,dec,cls=names[a[0].split()[0]]
        if cls not in {"Candidates TDE","Candidates ANT"}:continue
        rows.append(dict(source_id=a[0].replace(" #","-peak"),ra_deg=round(ra,8),dec_deg=round(dec,8),
                         utc=a[1][:10]+"T12:00:00Z",family="optical",anchor="fitted_peak_date",
                         catalog="Fink/ZTF Quintin",classification=cls,peak_error_days=re.search(r"pm([\d.]+)d",a[1])[1]))
    write("nuclear_events.csv",rows)

if __name__=="__main__":main()
