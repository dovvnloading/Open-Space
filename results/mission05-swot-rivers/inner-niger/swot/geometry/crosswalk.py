"""Retrieve official reach ID crosswalk for Africa by HTTP Range, verify ZIP CRC, and extract study rows."""
import sys, binascii, csv, hashlib, io, json, pathlib, struct, urllib.request, zipfile, zlib
from fetch_geom import reaches, out

nodes=len(sys.argv)>1
url='https://zenodo.org/records/22259077/files/SWORD_v17b_to_v16_Translation.zip'
def byte_range(lo,hi):
    req=urllib.request.Request(url,headers={'Range':f'bytes={lo}-{hi}'})
    with urllib.request.urlopen(req,timeout=90) as resp:
        assert resp.status==206 and resp.headers['Content-Range'].startswith(f'bytes {lo}-{hi}/')
        return resp.read()
size=265809613
end=byte_range(size-65536,size-1)
with zipfile.ZipFile(io.BytesIO(end)) as zipf:
    info=zipf.getinfo('v17b/AF_NodeIDs_v17b_vs_v16.csv' if nodes else 'v17b/AF_ReachIDs_v17b_vs_v16.csv')
start_dir=struct.unpack('<4s4H2LH',end[-22:])[6]
delta=start_dir-(len(end)-22-struct.unpack('<4s4H2LH',end[-22:])[5])
start=info.header_offset+delta
v=byte_range(start,start+29)
fields=struct.unpack('<4s5H3L2H',v)
assert fields[0]==b'PK\x03\x04'
first=start+30+fields[-2]+fields[-1]
member_cache=pathlib.Path(__file__).absolute().parent/(info.filename.split('/')[-1]+'.deflate')
if not member_cache.exists():
    member_cache.write_bytes(byte_range(first,first+info.compress_size-1))
data=zlib.decompress(member_cache.read_bytes(),-15)
assert binascii.crc32(data)==info.CRC
rows=csv.DictReader(io.TextIOWrapper(io.BytesIO(data),encoding='utf-8-sig'))
selected=[]
count=0
for row in rows:
    count+=1
    if count==1: print(row,flush=True)
    if nodes:
        target_ids=set(sys.argv[1:])
        if row['v16_node_id'] in target_ids: selected.append(row)
    elif row['v16_reach_id'] in reaches:
        selected.append(row)

meta={'source_url':url,'member':info.filename,'member_crc32':hex(info.CRC),'member_sha256':hashlib.sha256(data).hexdigest(),'member_rows':count,'record_doi':'10.5281/zenodo.22259077','selected':selected}
(out/('sword_node_crosswalk.json' if nodes else 'sword_crosswalk.json')).write_text(json.dumps(meta,indent=2)+'\n')
print(meta)
