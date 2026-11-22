"""Keep the BAT timing and count-rate definitions alongside the frozen inputs."""
import io
import json
from pathlib import Path
import sys

from astropy.io import fits
from ingest import BAT, fetch


def main(folder):
    url = BAT + "Crab.lc.fits"
    try:
        data, receipt = fetch(url)
        with fits.open(io.BytesIO(data), memmap=False) as hdus:
            h = next(h for h in hdus if isinstance(h, fits.BinTableHDU))
            out = dict(receipt=receipt, columns=[{"name":c.name,"unit":c.unit,"format":c.format} for c in h.columns],
                       timing={k:h.header[k] for k in ("MJDREF","MJDREFI","MJDREFF","TIMESYS","TIMEUNIT","TIMEZERO",
                               "TSTART","TSTOP","TELESCOP","INSTRUME") if k in h.header})
    except Exception as exc:
        out = dict(url=url, error=str(exc)[:450])
    (folder / "bat-product-metadata.json").write_text(json.dumps(out,indent=2)+"\n")


if __name__=="__main__":
    main(Path(sys.argv[1]))
