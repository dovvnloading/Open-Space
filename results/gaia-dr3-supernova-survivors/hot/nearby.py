"""Screen Dawson et al. (2026) 500-pc sample in both RV and Galactic rest-frame speed."""
import csv
import hashlib
import io
import json
import math
import urllib.parse
import urllib.request
from screen import HERE, VIZ


def get(n, cols):
    p={"-source":f"J/A+A/707/A6/table{n}", "-out":"recno,Name,"+cols, "-out.max":"unlimited"}
    url=VIZ+"?"+urllib.parse.urlencode(p)
    with urllib.request.urlopen(url,timeout=90) as response:
        b=response.read()
        receipt=dict(url=url,server_date=response.getheader("Date"),status=response.status,nbytes=len(b),sha256=hashlib.sha256(b).hexdigest())
    words=b.decode().splitlines()
    i=next(i for i,s in enumerate(words) if s.startswith("recno\t"))
    data=[{k:v.strip() for k,v in r.items()} for r in csv.DictReader([words[i],*[s for s in words[i+3:] if s.strip() and not s.startswith("#")]],delimiter="\t")]
    assert {int(a["recno"]) for a in data} == set(range(1,len(data)+1))
    with (HERE/f"nearby_table{n}.csv").open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]),lineterminator="\n");w.writeheader();w.writerows(data)
    return data, receipt


def main():
    spectra,r2=get(2,"Class,RV,Flag")
    kin,r4=get(4,"Uvel,E_Uvel,e_Uvel,Vvel,E_Vvel,e_Vvel,Wvel,E_Wvel,e_Wvel,Flag")
    assert len({r["Name"] for r in spectra}) == len(spectra)
    assert len({r["Name"] for r in kin}) == len(kin)
    rv=sorted(((abs(float(r["RV"])),r["Name"]) for r in spectra if r["RV"]),reverse=True)
    v=sorted(((math.sqrt(sum(float(r[k])**2 for k in ("Uvel","Vvel","Wvel"))),r["Name"])
              for r in kin if all(r[k] for k in ("Uvel","Vvel","Wvel"))),reverse=True)
    assert not any(speed>=700 for speed,name in rv) and not any(speed>=1000 for speed,name in v)
    meta=dict(paper="https://doi.org/10.1051/0004-6361/202558123", table2_receipt=r2,table4_receipt=r4,
              with_rv=len(rv),with_uvw=len(v),max_abs_rv=rv[0],max_galactic_rest_speed=v[0],
              count_abs_rv_ge_700=0,count_galactic_speed_ge_1000=0,
              conditional_binomial_95_fraction=1-0.05**(1/len(v)))
    (HERE/"nearby_summary.json").write_text(json.dumps(meta,indent=2)+"\n")
    print(meta)

if __name__=="__main__": main()
