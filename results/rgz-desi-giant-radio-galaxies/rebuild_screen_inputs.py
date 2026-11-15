"""Replay primary catalog queries and rebuild compact search inputs.

Run in an empty temporary directory with a copy of provenance.csv.
Catalog updates may change download digests, particularly timestamped query headers.
"""
import csv
import io
from pathlib import Path
import re
import tarfile
from urllib.request import Request, urlopen

from lxml import html
import numpy as np
import pandas as pd
from reproduce import crossmatch, separation

output = Path("inputs")
output.mkdir(exist_ok=True)
manifest = pd.read_csv("provenance.csv").set_index("input")
cache = Path("archive_cache")



def retrieve(key):
    cache.mkdir(exist_ok=True)
    path = cache/key.replace("/", "_")
    if not path.exists():
        with urlopen(Request(manifest.loc[key,"url"]),timeout=120) as response:
            path.write_bytes(response.read())
    return path


def catalog_positions(key):
    data = retrieve(key).read_text()
    if key == "Lao2026":
        return np.array([(float(s[47:56]),float(s[57:66])) for s in data.splitlines() if s.startswith("J")]).T
    if key == "Manik2026":
        rows = html.fromstring(data).xpath('//figure[@id="Sx7.T2"]//tr')[2:]
        positions = []
        for row in rows:
            cells = ["".join(c.itertext()).strip().replace("−","-") for c in row.xpath("./td")]
            positions.append((float(cells[2]),float(cells[3])))
        return np.array(positions).T
    if key == "Gordon2023":
        lines = [x for x in data.splitlines() if x.strip() and not x.startswith("#")]
        table = pd.read_csv(io.StringIO("\n".join([lines[0], *lines[3:]])),sep="\t",skipinitialspace=True)
        return table[["_RA","_DE"]].to_numpy().T
    coords = re.findall(r"^\s*([\d.]+)\s*\t\s*([+-]?[\d.]+)\s*$",data,re.M)
    return np.array(coords,dtype=float).T


def query_batches(prefix):
    frames=[]
    for key in manifest.index:
        if not key.startswith(prefix+"/"):
            continue
        frame = pd.read_csv(retrieve(key))
        if len(frame):
            frame["batch"]=Path(key).stem
            frames.append(frame)
    return pd.concat(frames)


def main():
    with tarfile.open(retrieve("rgz.tgz")) as tar:
        radio = pd.read_csv(tar.extractfile("DR1_FIRST_radio_classifications.csv"))
        optical = pd.read_csv(tar.extractfile("DR1_FIRST_host_properties.csv"))
    joined = radio.merge(optical,left_on="CatID",right_on="#CatID")
    select = joined.query("LAE >= 100 and N_comp >= 2 and CL >= 0.65 and Host_RA > 0")
    hosts = select.sort_values(["LAE","CatID"],ascending=[False,True]).drop_duplicates(["Host_RA","Host_Dec"])
    hosts = hosts.rename(columns={"RGZID_x":"rgz","Host_RA":"ra","Host_Dec":"dec"})
    hosts[["CatID","rgz","ra","dec","CL","LAE","N_comp","z_phot0"]].to_csv(output/"selection.csv",index=False)
    spectra = query_batches("rgz")
    spectra.to_csv(output/"desi.csv",index=False)
    query_batches("controls").to_csv(output/"shifted_desi.csv",index=False)
    matches = crossmatch(hosts,spectra)
    matched_hosts = hosts.set_index("CatID").loc[matches.CatID.unique()]
    catalog_keys = [key for key in manifest.index if key.startswith("k_") or key in ("cdsgrq","Lao2026","Manik2026","Gordon2023")]
    neighbors=[]
    for key in catalog_keys:
        ra,dec = catalog_positions(key)
        for host in matched_hosts.itertuples():
            sep = separation(host.ra,host.dec,ra,dec)
            i = np.argmin(sep)
            neighbors.append((host.Index,key,len(ra),ra[i],dec[i],round(sep[i],3)))
    pd.DataFrame(neighbors,columns=["CatID","catalog","catalog_entries","nearest_ra","nearest_dec","separation_arcsec"]).to_csv(output/"catalog_neighbors.csv",index=False)


if __name__ == "__main__":
    main()
