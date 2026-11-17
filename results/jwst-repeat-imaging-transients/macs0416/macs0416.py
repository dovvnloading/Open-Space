"""Search a matched F150W/F444W pair around the Warhol and Spock arcs."""
import argparse
import importlib.util
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
LEGACY = HERE.parent / "el-gordo"
CENTER = (64.0375, -24.0686)
EPOCHS = ("jw01176-o213_t009", "jw06882-o054_t054")
BANDS = ("f150w", "f444w")

def module(name, base):
    spec = importlib.util.spec_from_file_location(name, base / (name + ".py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("action",choices=("fetch","measure"));p.add_argument("cache",type=Path)
    p.add_argument("--base",type=Path,default=LEGACY)
    a=p.parse_args();a.cache.mkdir(exist_ok=True, parents=True)
    if a.action=="fetch":
        from concurrent.futures import ThreadPoolExecutor
        sys.argv=[sys.argv[0],str(a.cache)]
        m=module("fetch_elgordo",a.base);m.CENTER=CENTER;m.SIZE_ARCSEC=39
        with ThreadPoolExecutor(max_workers=2) as pool:
            rows=list(pool.map(m.fetch,[(e,b) for e in EPOCHS for b in BANDS]))
        (HERE/"input_ranges.json").write_text(json.dumps(rows,indent=1)+"\n")
    else:
        m=module("analyze_elgordo",a.base);m.HERE=HERE;m.CENTER=CENTER;m.SIDE=601
        # The original screen expects exactly two files per filter. Isolate the
        # discovery pair when auxiliary context epochs share its input cache.
        from tempfile import TemporaryDirectory
        with TemporaryDirectory() as tmp:
            for epoch in EPOCHS:
                for band in BANDS:
                    name=f"{epoch}_nircam_clear-{band}.fits"
                    (Path(tmp)/name).symlink_to((a.cache/name).resolve(strict=True))
            sys.argv=[sys.argv[0],tmp]
            m.main()

if __name__=="__main__":main()
