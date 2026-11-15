"""Reproduce the homogeneous period screen used to select the Kepler-359 test."""
import argparse
import csv
import hashlib
import itertools
import math
from collections import defaultdict
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parent
SOURCE = "https://content.cld.iop.org/journals/2632-3338/5/6/152/revision1/psjad0e6et1_mrt.txt"
SHA = "1fe1360313268ee305dead10cffbf79e1ebf18a76b0032e64e11870fe3af0521"


def main(path=None):
    body = Path(path).read_bytes() if path else urlopen(SOURCE, timeout=90).read()
    if hashlib.sha256(body).hexdigest() != SHA:
        raise ValueError("Table 1 checksum differs from the source snapshot")
    hosts = defaultdict(list)
    for line in body.decode().splitlines():
        if len(line) < 570 or not line[:8].strip().isdigit() or line[547] != "P":
            continue
        period = float(line[32:46])
        if period > 0:
            hosts[line[:8].strip()].append((period,line[17:31].strip() or line[9:16].strip()))
    indices = [(a,b) for a,b in itertools.product(range(1,7),repeat=2)
               if math.gcd(a,b)==1 and a+b<=10]
    hits=[]
    triple_count=0
    searched_hosts=set()
    for kic, planets in hosts.items():
        for inner,middle,outer in itertools.combinations(sorted(planets),3):
            p,q,s=inner[0],middle[0],outer[0]
            if max(q/p,s/q)>4:
                continue
            searched_hosts.add(kic)
            triple_count+=1
            for a,b in indices:
                speed = a/p-(a+b)/q+b/s
                detuning=abs(speed)*q/(a+b)
                hits.append((detuning,kic,inner[1],middle[1],outer[1],a,-a-b,b,speed))
    hits.sort()
    upto=next(i for i,v in enumerate(hits,1) if v[2:5]==('Kepler-359 b','Kepler-359 c','Kepler-359 d') and v[5:8]==(1,-6,5))
    print(f'{len(searched_hosts)} hosts; {triple_count} eligible triplets; {len(indices)} vectors; {len(hits)} tests; Kepler-359 rank {upto}')
    with (ROOT/'period_screen.csv').open('w') as handle:
        writer=csv.writer(handle,lineterminator='\n')
        writer.writerow(['rank','scaled_frequency_offset','kic','inner','middle','outer','k1','k2','k3','cycles_per_day'])
        for rank,v in enumerate(hits[:upto],1):
            writer.writerow((rank,*v))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table',help='Local copy of the original machine-readable table')
    main(parser.parse_args().table)
