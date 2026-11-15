"""Rebuild the Wind/Elenin data packet with python + numpy, pandas and scipy.

Run: python replay_wind.py --output OUTPUT --cache CACHE
Raw public API responses are cached. The VDF is aggregated only after explicit
shape and fill checks. No missing species moment is converted to a density.
"""
import argparse
import hashlib
import io
import json
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import minimize_scalar


AU_KM = 149597870.7
HAPI = "https://cdaweb.gsfc.nasa.gov/hapi/data?"
HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api?"
MOMENTS = "WI_L2_3MIN_SMS-STICS-NVT-SOLARWIND"
VDF = "WI_L2-3MIN_SMS-STICS-VDF-SOLARWIND"
MOMENT_VARS = "n_tc_hplus,n_tc_oplus,n_tc_o6plus,n_dc_heplus,n_dc_oplus"
SWE_VARS = "Proton_VX_nonlin,Proton_VY_nonlin,Proton_VZ_nonlin"
CHANNELS = ["h_tc", "o_tc", "o6_tc", "he_dc", "o_dc"]


class Inputs:
    def __init__(self, cache):
        self.cache = cache
        self.manifest = []
        cache.mkdir(parents=True, exist_ok=True)

    def fetch(self, name, url):
        path = self.cache / name
        if not path.exists():
            with urllib.request.urlopen(url, timeout=180) as response:
                data = response.read()
            path.write_bytes(data)
        data = path.read_bytes()
        self.manifest.append(dict(file=name, url=url, bytes=len(data),
                                  sha256=hashlib.sha256(data).hexdigest()))
        return data

    def hapi(self, name, dataset, variables, start, stop):
        query = dict(id=dataset, parameters=variables, format="csv")
        query.update({"time.min": start.isoformat().replace("+00:00", "Z"),
                      "time.max": stop.isoformat().replace("+00:00", "Z")})
        data = self.fetch(name, HAPI + urllib.parse.urlencode(query))
        if not data.strip():
            return pd.DataFrame()
        frame = pd.read_csv(io.BytesIO(data), header=None)
        frame.index = pd.to_datetime(frame.pop(0), utc=True)
        return frame

    def ephemeris(self, name, command):
        params = dict(COMMAND=command, EPHEM_TYPE="VECTORS", CENTER="500@10",
                      START_TIME="2011-09-15", STOP_TIME="2011-10-10",
                      STEP_SIZE="6 h", OUT_UNITS="AU-D", VEC_TABLE="2",
                      VEC_CORR="NONE", REF_PLANE="ECLIPTIC", CSV_FORMAT="YES",
                      OBJ_DATA="NO")
        query = {"format": "text", **{k: "'" + v + "'" for k,v in params.items()}}
        raw = self.fetch(name + ".txt", HORIZONS + urllib.parse.urlencode(query)).decode()
        block = raw.split("$$SOE")[1].split("$$EOE")[0].strip()
        states = np.array([[float(line.split(",")[i]) for i in (0,2,3,4,5,6,7)]
                           for line in block.splitlines()])
        return CubicHermiteSpline(states[:, 0], states[:, 1:4], states[:, 4:7],
                                  extrapolate=False)


def utc(date):
    return pd.Timestamp(date, tz="UTC")


def moments(source, start, stop):
    days = pd.date_range(utc(start), utc(stop), inclusive="left")
    pieces = []
    for day in days:
        frame = source.hapi(str(day.date()) + ".csv", MOMENTS, MOMENT_VARS,
                            day, day + pd.Timedelta(days=1))
        if len(frame):
            frame.columns = CHANNELS
            pieces.append(frame)
    result = pd.concat(pieces).sort_index()
    if result.index.has_duplicates:
        raise ValueError("overlapping STICS moments")
    return result


def window_counts(mom, hours):
    # A positive moment is a detection in a native ~3-minute interval. Report
    # those separately from the integer instrument counts in the full VDF.
    binned = mom.gt(0).resample(hours).sum()
    binned.insert(0, "records", mom.resample(hours).size())
    return binned


def miss(time, gse_velocity, wind, earth, comet):
    # UTC to TDB via TT, appropriate to <2 ms in this 2011 interval.
    jd = time.to_julian_date() + 66.184 / 86400
    x = -earth(jd)
    x /= np.linalg.norm(x)
    y = np.cross([0., 0., 1.], x)
    y /= np.linalg.norm(y)
    z = np.cross(x, y)
    v = gse_velocity @ np.stack([x, y, z]) * 86400 / AU_KM
    # Wind/SWE already subtracts orbital aberration. V is rotated only.
    position = wind(jd)
    def distance(flight_days):
        delta = position - flight_days * v - comet(jd - flight_days)
        return np.linalg.norm(delta)
    closest = minimize_scalar(distance, bounds=(0, 3), method="bounded")
    return closest.fun * AU_KM, closest.x


def vector_records(source, start, stop, filename):
    vdf = source.hapi(filename, VDF, "counts_tc_oplus,counts_dc_oplus", start, stop)
    energy = source.hapi(start.strftime("%Y_%m_%d_%H")+"-energy.csv", VDF, "eoq",
                         start, start+pd.Timedelta(minutes=30))
    if not np.allclose(energy.iloc[0],energy.iloc[-1]):
        raise ValueError("inconsistent energy axes")
    values = vdf.to_numpy().reshape(len(vdf), 2, 16, 3, 32)
    positive = np.argwhere(values > 0)
    events = []
    for row, mode, sector, telescope, channel in positive:
        n = values[row, mode, sector, telescope, channel]
        if not np.isclose(n,round(n)):
            raise ValueError("non-integral VDF counts")
        events.append(dict(time_utc=vdf.index[row].isoformat(), mode=("TC","DC")[mode],
                           sector_0based=sector, telescope_1based=telescope+1,
                           energy_step_1based=channel+1, energy_kev_per_e=energy.iloc[0,channel],
                           counts=int(round(n))))
    return pd.DataFrame(events), vdf.index


def run(cache, output):
    source = Inputs(cache)
    output.mkdir(parents=True,exist_ok=True)
    encounter = moments(source, "2011-09-19", "2011-10-09")
    hours = window_counts(encounter,"1h")
    swe = source.hapi("all_swe.csv","WI_H1_SWE",SWE_VARS,utc("2011-09-19"),utc("2011-10-09"))
    swe.columns = ["vx_gse", "vy_gse", "vz_gse"]
    swe = swe.where(swe.abs() <= 1800).dropna()
    hours = hours.join(swe.resample("1h").median())
    hours["swe_records"] = swe.resample("1h").size()
    wind, earth, comet = [source.ephemeris(n,c) for n,c in [
        ("wind","-8"),("earth","399"),("elenin","DES=2010 X1;CAP;")]]
    for t,row in hours.iterrows():
        if row["swe_records"] >= 10:
            v = row[["vx_gse","vy_gse","vz_gse"]].to_numpy(dtype=float)
            distance, flight = miss(t + pd.Timedelta(minutes=30),v,wind,earth,comet)
            hours.loc[t, "miss_km"] = distance
            hours.loc[t, "flight_days"] = flight
    hours.to_csv(output/"wind_hourly.csv",index_label="start_utc",float_format="%.5f")
    bins = window_counts(encounter,"3h")
    bins["rate"] = bins.o_dc/bins.records
    bins["median_miss_km"] = hours.miss_km.resample("3h").median()
    bins.to_csv(output/"wind_three_hour.csv",index_label="start_utc",float_format="%.6f")
    event = utc("2011-09-29 06:00")
    valid = bins[bins.records>=45]
    controls = valid.drop(event)
    selected = valid.loc[event]
    edges = np.arange(0,len(bins))
    rates = bins.rate.where(bins.records>=45).to_numpy()
    close = (bins.median_miss_km <= 2e6).to_numpy()
    # Circular phase shifts on the fixed calendar, preserving all gaps,
    # close-geometry bin intervals and relative timing within the series.
    scores=[]
    for shift in edges:
        region = np.roll(rates,shift)[close]
        scores.append(float(np.nanmax(region)) if np.isfinite(region).any() else np.nan)
    report=dict(window_start_utc=event.isoformat(), minimum_STICS_records_per_3h=45,
                selected_window=selected.to_dict(),
                covered_comparison_bins=int(len(controls)),
                controls_at_least_selected_rate=int((controls.rate>=selected.rate).sum()),
                geometry_bins=int(close.sum()), shifts=len(scores),
                finite_shifts=int(np.isfinite(scores).sum()),
                shifts_at_least_observed_max=int((np.array(scores)>=scores[0]).sum()),
                observed_max_rate_inside_geometry=scores[0])
    # All 2011 native O+ moments, retaining calendar gaps and three bin origins.
    annual = []
    for month in pd.date_range(utc("2011-01-01"), utc("2011-12-01"), freq="MS"):
        frame = source.hapi(month.strftime("%Y-%m")+".csv", MOMENTS,
                            "n_tc_oplus,n_dc_oplus", month, month+pd.offsets.MonthBegin(1))
        frame.columns = ["o_tc", "o_dc"]
        annual.append(frame)
    annual = pd.concat(annual).sort_index()
    if annual.index.has_duplicates:
        raise ValueError("overlap between monthly rows")
    ranked = []
    for offset in (0,1,2):
        native = annual.gt(0).resample("3h",offset=pd.Timedelta(hours=offset))
        stats = native.sum()
        stats.insert(0,"records",native.size())
        stats["rate"] = stats.o_dc / stats.records
        if offset == 0:
            stats.to_csv(output/"wind_2011_three_hour.csv",index_label="start_utc",float_format="%.6f")
        for minimum in (1,30,45,50):
            table = stats[stats.records>=minimum].sort_values("rate",ascending=False)
            top = table.iloc[0]
            ranked.append(dict(offset_hours=offset,minimum_records=minimum,
                               covered_bins=len(table),top_bin=table.index[0].isoformat(),
                               occupied=int(top.o_dc), records=int(top.records), fraction=float(top.rate),
                               tied_maxima=int((table.rate==top.rate).sum())))
    report["annual_ranks"] = ranked

    # Nearby spacecraft's contemporaneous (unshifted) plasma at the Wind
    # release point is a sensitivity check, not a second detection of oxygen.
    ace = source.hapi("ACE_GSE.csv","AC_H0_SWE","V_GSE",utc("2011-09-29"),utc("2011-09-30"))
    ace.columns = ["vx_gse","vy_gse","vz_gse"]
    ace = ace.where(ace.abs()<=1800).dropna()
    med = ace.resample("1h").median()
    med["records"] = ace.resample("1h").size()
    for t,row in med.iterrows():
        if row.records >= 10:
            dist,lag = miss(t+pd.Timedelta(minutes=30),row.iloc[:3].values,wind,earth,comet)
            med.loc[t,"miss_at_wind_km"] = dist
            med.loc[t,"flight_days"] = lag
    med.to_csv(output/"ace_flow_sensitivity.csv",index_label="start_utc",float_format="%.5f")
    pesa = source.hapi("3DP.csv","WI_PM_3DP","P_VELS",utc("2011-09-29 05:00"),
                       utc("2011-09-29 11:00"))
    pesa.columns = ["vx_gse", "vy_gse", "vz_gse"]
    pesa = pesa.where(pesa.abs()<=1800).dropna()
    med = pesa.resample("1h").median()
    med["records"] = pesa.resample("1h").size()
    for t,row in med.iterrows():
        v = row.iloc[:3].to_numpy(dtype=float)
        epoch = t+pd.Timedelta(minutes=30)
        med.loc[t,"miss_as_delivered_km"] = miss(epoch,v,wind,earth,comet)[0]
        med.loc[t,"miss_minus_29p78_y_km"] = miss(epoch,v-[0,29.78,0],wind,earth,comet)[0]
    med.to_csv(output/"pesa_flow_sensitivity.csv",index_label="start_utc",float_format="%.5f")
    ev,times = vector_records(source,utc("2011-09-29 05:00"),utc("2011-09-29 11:00"),"El_VDF.csv")
    comparisons=[]
    vwindows = [("2011-09-29 06:00", "Elenin geometry"), ("2011-09-26 18:00","ICME interval"),
                ("2011-04-30 03:00","Third-ranked 2011 bin"),("2006-06-02 06:00","73P control"),
                ("2006-06-02 15:00","73P control"),("2017-01-31 12:00","2017-01-31 CIR context")]
    for stamp, label in vwindows:
        begin = utc(stamp)
        end = begin+pd.Timedelta(hours=3)
        if begin==event:
            tm = pd.to_datetime(ev.time_utc)
            hits = ev[(tm>=begin)&(tm<end)]
            nrecords = int(((times>=begin)&(times<end)).sum())
        else:
            hits, epochs = vector_records(source, begin, end, begin.strftime("%Y_%m_%d_%H")+"-VDF.csv")
            nrecords = len(epochs)
        for mode in ("TC","DC"):
            rows = hits[hits["mode"]==mode]
            count = int(rows.counts.sum())
            energy = np.repeat(rows.energy_kev_per_e,rows.counts)
            comparisons.append(dict(start_utc=begin.isoformat(),label=label,mode=mode,
                                    scans=nrecords,counts=count,sector7_12_counts=int(rows.loc[
                                    rows.sector_0based.between(7,12),"counts"].sum()),
                                    central_telescope_counts=int(rows.loc[rows.telescope_1based==2,"counts"].sum()),
                                    median_kev_per_e=float(np.median(energy)) if count else np.nan))
    pd.DataFrame(comparisons).to_csv(output/"wind_vdf_controls.csv",index=False)
    ev.to_csv(output/"wind_oxygen_events.csv",index=False)
    report["vdf_3min_records_05_11"]=int(len(times))
    report["vdf_3min_records_06_09"]=int(((times>=event)&(times<event+pd.Timedelta(hours=3))).sum())
    control = moments(source, "2006-05-17", "2006-06-20")
    window_counts(control,"1D").to_csv(output/"wind_73p_control.csv",index_label="start_utc")
    (output/"wind_metrics.json").write_text(json.dumps(report,indent=2)+"\n")
    pd.DataFrame(source.manifest).to_csv(output/"wind_sources.csv",index=False)


if __name__ == "__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cache",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    run(args.cache,args.output)
