import io, json, sys, zipfile, urllib.request, urllib.parse, hashlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ra,dec=15.7324167,-72.8732778
base='https://mast.stsci.edu/tesscut/api/v0.1/'
def fetch(s):
    url=base+'astrocut?'+urllib.parse.urlencode(dict(ra=ra,dec=dec,sector=s,x=11,y=11,units='px'))
    fn=Path(f's{s:03}.fits')
    if not fn.exists():
        resp=urllib.request.urlopen(url,timeout=130).read()
        with zipfile.ZipFile(io.BytesIO(resp)) as z:
            [name]=[n for n in z.namelist() if n.endswith('.fits')]
            fn.write_bytes(z.read(name))
    print(s,fn.stat().st_size, hashlib.sha256(fn.read_bytes()).hexdigest(), flush=True)
if __name__=='__main__':
    with ThreadPoolExecutor(2) as ex: list(ex.map(fetch, [int(i) for i in sys.argv[1:]]))
