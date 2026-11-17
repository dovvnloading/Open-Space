"""Fail-closed primitives for immutable microlensing replay campaigns."""
import hashlib,json,os,platform,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parent
REPO=ROOT.parents[1]
FROZEN={'kmt/pySIS-2024.tsv':'863f215cf1a908c19148c8a4de15cd94d13474fe8b55fa76712ee23517e7391f','kmt/pySIS-2025.tsv':'d6e3edfdf28b29c43b62df29a4e4b240db9c2e1cf4046249ee1fddde1b2f23ed','ogle/catalog.csv':'7e9be993664be633f5225d72726ffe2bcafc557a4b7920c81e88f09b9d2edd1c','kmt/catalog-2024.tsv':'77187fe914f5da33a0ae8169ebab7dcc3fb3f3977f2231e302d9b1fedc6a9bb5','kmt/catalog-2025.tsv':'e2947ff7578590cd48192fdbb5fdc236bf17d9207835f50b7af1b5ba45bbb791','tables/cadence.csv.gz':'9d2d56ceebb1512da151edb1b23acbae02d7f3a02decafb68e24cd0ec31b2bdb','tables/peak-comparisons.json.gz':'3e87ce4e395e73a70b39b0f2d2b5e8247e5fa3c0d046b7fc1a6c04c426376567','other/gaia.csv':'52a394f0df36e6baf6b3c31909dd43d9c5abf0a2a994fa5c28b420580fdde80a','other/gaia_index_receipt.json':'3b86fabcf9470a26ea89275641b5dce8b910730c538d15b58861171e8330a8f9','ogle/tables.json':'ee71ddb912e1d60b8f2eee228f58926e238ee93c78f83904fb726b4a3da0811a','ogle/photometry.csv':'4cdf77ca3ca3cc98bf180084d8097b5838f95a6e10dc85ea3ce2f8df49c6dd6a','kmt/index-receipts.json':'6c60b82229d10a228bc0a2e737312f00a1f1dea1ff0cef1cb13eff5c57afd37b','other/index_receipts.csv':'2a632070d11636b52a06717a8ca47b40119b122ccd40e4e7bcd692395b8ae075','other/photometry.csv':'ad458410df492307f0b7beb73645142118ab07f5403c5e69d05dff749dfeb6d0','tables/grid-fits.csv.gz':'dafc904395f8a0771c2beda9c820b1f3cef216a95b734035a0b841b6ece234af'}
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def external(path):
 p=Path(path)
 if not p.is_absolute():raise ValueError('absolute D: path required')
 p=p.resolve()
 allowed=(p.drive.upper()=='D:') if os.name=='nt' else (p==Path('/mnt/d') or Path('/mnt/d') in p.parents)
 if not allowed or p==REPO or REPO in p.parents:raise ValueError('external D: path required')
 return p
def fresh(path):
 p=external(path)
 if p.exists():raise FileExistsError('new exclusive output required: '+str(p))
 p.mkdir(parents=True,exist_ok=False);return p
def verified(path,expected):
 p=Path(path)
 if not p.is_file() or sha(p)!=expected:raise ValueError('frozen input mismatch: '+str(p))
 return p
def identity(path,event,cap,code_hashes,runtime,parameters=None):
 return {'event':event,'input_sha256':sha(path),'cap':bool(cap),'code_sha256':dict(code_hashes),'runtime':runtime,'parameters':parameters}
def check_checkpoint(record,expected):
 if record.get('_identity')!=expected:raise ValueError('stale or unbound checkpoint')
 return record

def verify_receipt(name):
 return verified(ROOT/name,FROZEN[name])

def fresh_for_inputs(out,*inputs):
 dest=external(out)
 for item in inputs:
  src=external(item)
  if dest==src or dest in src.parents or src in dest.parents:raise ValueError('input/output overlap')
 return fresh(dest)
