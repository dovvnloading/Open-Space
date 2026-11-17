"""Test the automatic fitter on the three known visitors with their names blinded."""
import base64
import gzip
import hashlib
import json
import os
from pathlib import Path

from refit import cat, fit

ROOT=Path(__file__).resolve().parent
OBJECTS={"1I":"50788063","2I":"1003639","3I":"1004083"}


def main():
    catalog,_,_=cat()
    by_id={str(row["spkid"]):row for row in catalog}
    output=[]
    for iso,key in OBJECTS.items():
        source=ROOT.parent/"orbits/data"/(iso+".obs.gz")
        raw=gzip.decompress(source.read_bytes())
        # Select ground optical positions and sample uniformly in sorted record rank.
        # Each variant below gets exactly the same time, station and sky measurement.
        optical=sorted((l for l in raw.splitlines() if len(l)>=80 and l[14:15] in (b"C",b"B")
                        and l[15:19].isdigit()),key=lambda l:(l[15:32],l[77:80]))
        n=min(180,len(optical))
        inds=list(dict.fromkeys(round(j*(len(optical)-1)/(n-1)) for j in range(n)))
        chosen=[optical[j][:80]+b"\n" for j in inds]
        assert len(chosen)>=5
        rows={"label":iso,"source_sha256":hashlib.sha256(raw).hexdigest(),"eligible":len(optical),
              "selected":len(chosen),"variants":[]}
        for name,prefix,unpenalized in [("known permanent ID",("000"+iso+" "*7).encode(),False),
                                        ("blinded with default score",b"     Blind01",False),
                                        ("blinded without hyperbola penalty",b"     Blind01",True)]:
            assert len(prefix)==12
            lines=[prefix+l[12:] for l in chosen]
            os.environ["FO_UNPENALIZED"]="1" if unpenalized else "0"
            res=fit((by_id[key],lines))
            rows["variants"].append({"mode":name,"result":res})
            print("BLIND",iso,name,res.get("status"),res.get("e"),res.get("used"),flush=True)
        output.append(rows)
    data=json.dumps(output,indent=2).encode()
    out=ROOT/"evidence/blinded_controls.json"
    out.write_bytes(data)
    print("ARCHIVE",out.name,base64.b64encode(data).decode())


if __name__=="__main__":
    main()
