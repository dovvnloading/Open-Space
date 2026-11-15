"""Independently verify released measurements against public archive products."""
import csv
import hashlib
import io
from pathlib import Path
from urllib.request import Request, urlopen
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
HARPS_URL = ("https://raw.githubusercontent.com/3fon3fonov/HARPS_RVBank/"
             "09626a652ddc93df8771ce3d818ee04a292f7585/HARPS_RVBank_ver02.csv.zip")
ARCHIVE_SHA = "9218ebd833f8971dcf304c7a6d9835de1c988dc5faae131f3eb939c7b9682586"
KEYWORDS = {"bjd":"QC BJD", "rv":"QC CCF RV", "rv_err":"QC CCF RV ERROR",
            "bis":"QC CCF BIS SPAN", "fwhm":"QC CCF FWHM",
            "contrast":"QC CCF CONTRAST", "berv":"QC BERV"}


def read_csv(path):
    with path.open(newline="") as source:
        return list(csv.DictReader(source))


def verify_harps():
    with urlopen(HARPS_URL, timeout=90) as response:
        data = response.read()
    assert hashlib.sha256(data).hexdigest() == ARCHIVE_SHA
    with ZipFile(io.BytesIO(data)) as archive:
        with archive.open("HARPS_RVBank_ver02.csv") as source:
            rows = [row for row in csv.DictReader(io.TextIOWrapper(source))
                    if row["target"] == "HD32564"]
    local = read_csv(ROOT/"harps.csv")
    assert len(rows) == len(local)
    for archived, released in zip(rows, local):
        for name, value in released.items():
            if name in {"target", "TIMEID", "PROGID"}:
                assert archived[name] == value
            else:
                assert abs(float(archived[name]) - float(value)) < 1e-9


def verify_espresso():
    for row in read_csv(ROOT/"espresso.csv"):
        url = "https://dataportal.eso.org/dataPortal/file/" + row["product"]
        with urlopen(Request(url, headers={"Range":"bytes=0-172799"}), timeout=60) as source:
            header = source.read(172800)
        cards = [header[i:i+80] for i in range(0, len(header), 80)]
        end = next(i for i, card in enumerate(cards) if card.startswith(b"END     "))
        cards = cards[:end+1]
        assert hashlib.sha256(b"".join(cards)).hexdigest() == row["fits_header_sha256"]
        fields = {}
        for card in cards:
            name, sep, value = card.decode("ascii").partition("=")
            if sep:
                fields[name.strip()] = value.split("/",1)[0].strip().strip("'").strip()
        for name, keyword in KEYWORDS.items():
            assert abs(float(row[name]) - float(fields["HIERARCH ESO " + keyword])) < 1e-9
        assert row["programme"] == fields["HIERARCH ESO OBS PROG ID"]


if __name__ == "__main__":
    verify_harps()
    verify_espresso()
    print("Verified 205 HARPS rows and 25 ESPRESSO FITS headers.")
