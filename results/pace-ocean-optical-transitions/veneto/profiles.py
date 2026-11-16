"""Select the shallowest valid reading from each station's November ARPAV cast."""

from hashlib import sha256
from io import BytesIO
from pathlib import Path
from urllib.request import urlopen
import sys
import pandas as pd

URL = "https://www.arpa.veneto.it/dati-ambientali/open-data/file-e-allegati/soaml/profili-sonda/sonda-mare/sonda_mare_2025.csv/@@download/file"
HASH = "7a16ca21bae9d9ba34b1c09467ff2ad2187dbbba4787fe707df99064035f4aa0"
COLUMNS = {
    "Cod_Naz": "station",
    "Data": "date",
    "Orario": "local_time_as_reported",
    "Quota misura (m)": "sample_depth_m",
    "Latitudine_GBO": "GBO_northing",
    "Longitudine_GBO": "GBO_easting",
    "Temperatura acqua (gradi C)": "water_temp_C",
    "Salinita' (PSU)": "salinity_PSU",
    'Clorofilla "a" (ug/l)': "chl_a_ug_L",
    "Torbidita' (FTU)": "turbidity_FTU",
}


def main():
    if len(sys.argv) > 1:
        raw = Path(sys.argv[1]).read_bytes()
    else:
        with urlopen(URL, timeout=60) as response:
            raw = response.read()
    if sha256(raw).hexdigest() != HASH:
        raise ValueError("The ARPAV source differs from the analyzed 2025 release.")
    data = pd.read_csv(BytesIO(raw)).rename(columns=COLUMNS)[list(COLUMNS.values())]
    data["date"] = pd.to_datetime(data.date, format="%d/%m/%Y")
    in_month = data[data.date.dt.month == 11].dropna(
        subset=["sample_depth_m", "chl_a_ug_L", "turbidity_FTU"]
    )
    surface = in_month.sort_values(
        ["date", "station", "sample_depth_m"]
    ).drop_duplicates(["station", "date"])
    surface.to_csv(Path(__file__).with_name("veneto_casts.csv"), index=False)


if __name__ == "__main__":
    main()
