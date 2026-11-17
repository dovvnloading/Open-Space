"""Re-run frozen ADQL inputs via NOIRLab TAP. Positional output defaults to cwd."""
import argparse, gzip
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import urlopen
HERE=Path(__file__).resolve().parent
def main():
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);p.add_argument('adql',nargs='*',type=Path);args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    for qfile in args.adql or sorted(HERE.glob('*.adql')):
        q=qfile.read_text().strip()
        url='https://datalab.noirlab.edu/tap/sync?'+urlencode({'REQUEST':'doQuery','LANG':'ADQL','FORMAT':'csv','MAXREC':200000,'QUERY':q})
        with urlopen(url,timeout=600) as f: result=f.read()
        assert not result.lstrip().startswith(b'<'),result[:200]
        if qfile.stem=='control':
            (args.output/'control.csv.gz').write_bytes(gzip.compress(result,mtime=0))
        else:(args.output/(qfile.stem+'.csv')).write_bytes(result)
if __name__=='__main__':main()
