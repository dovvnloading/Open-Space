// Turn the published table rows into native observations and display coordinates.
// The optional Gaia comparison changes plotted coordinates only.
const DAY = 86400000;
const toDay = (isoDate) => Date.parse(isoDate + "T00:00:00Z") / DAY;
const number = new Intl.NumberFormat("en-US", {maximumFractionDigits: 9});
const utcDate = new Intl.DateTimeFormat("en-GB", {day: "numeric", month: "short", year: "numeric", timeZone: "UTC"});
const shortDate = (isoDate) => utcDate.format(new Date(toDay(isoDate) * DAY));
const point = (series, value, spread = 0) => ({series, value: Number(value), spread: Number(spread)});

const highestFirstSeries = rows => rows.reduce((best, row, index) =>
  row.points[0].value > rows[best].points[0].value ? index : best, 0);

// Pixel bounds of the four published 16-arcsecond panels in four-epochs.png.
// The viewer regression test pins the source image digest alongside these coordinates.
export const radioCutouts = {
  "2017-10-11": 90, "2020-08-25": 467, "2023-02-09": 844, "2025-09-06": 1220,
};

const profiles = {
  mars: {
    title: "Particle flux in Mars orbit",
    note: "Each point is a daily maximum of hourly median flux. Shading marks days with hours that meet the survey’s event criteria.",
    series: ["Electrons", "Ions"],
    xLabel: "August 2025 · UTC day",
    yLabel: "Flux · (cm² s sr)⁻¹ · log scale",
    log: true,
    row([date, observed, selected, ions, electrons]) {
      return {
        x: toDay(date), stamp: shortDate(date), axisTick: date.slice(8), eventDay: Number(selected) > 0,
        points: [point(0, electrons), point(1, ions)],
        fields: [["Electrons", number.format(electrons)], ["Ions", number.format(ions)], ["Selected hours", `${selected} / ${observed}`]],
      };
    },
    initial: highestFirstSeries,
  },
  sun: {
    title: "The hard X-ray flare",
    note: "One science accumulation per point, recorded aboard Solar Orbiter on 13 July 2021.",
    series: ["15–25 keV"],
    xLabel: "13 July 2021 · spacecraft UTC",
    yLabel: "Count rate · counts/s",
    connect: true,
    zero: true,
    row([time, rate]) {
      const [h, m, s] = time.split(":").map(Number);
      return {
        x: h * 3600 + m * 60 + s, stamp: `${time} UTC`,
        points: [point(0, rate)], fields: [["Count rate", `${rate} counts/s`]],
      };
    },
    initial: highestFirstSeries,
  },
  quasar: {
    title: "One quasar in two images",
    note: "Gaia G-band nights at their observed epochs. Smaller magnitudes are brighter. Bars show median catalog errors.",
    alignedNote: "Image B is shifted by +210 days and −1.94 mag, as in the published figure. This tests a proposed delay.",
    series: ["Image A", "Image B"],
    xLabel: "Epoch · BJD TCB − 2,400,000.5",
    yLabel: "Gaia G · mag · brighter upward",
    inverted: true,
    row([image, epoch, magnitude, error], aligned) {
      const shift = aligned && image === "B";
      return {
        x: Number(epoch) + (shift ? 210 : 0), stamp: `Image ${image}`, axisTick: number.format(epoch),
        points: [point(image === "A" ? 0 : 1, Number(magnitude) - (shift ? 1.94 : 0), error)],
        fields: [["Native epoch", epoch], ["Observed G (mag)", magnitude], ["Median error (mag)", error]],
      };
    },
    initial: rows => {
      const imageB = rows.filter(r => r.points[0].series === 1);
      return rows.indexOf(imageB.reduce((best, row) => best.points[0].value < row.points[0].value ? best : row));
    },
  },
  dwarf: {
    title: "Four visits in radio light",
    note: "Beam fits at the predicted position. Bars span the local control scatter. The 2020 fit is consistent with noise.",
    series: ["Forced flux"],
    xLabel: "VLASS observation · UTC",
    yLabel: "Flux density · mJy/beam",
    zero: true,
    row([date, flux, scatter, centroid]) {
      return {
        x: toDay(date), stamp: shortDate(date), axisTick: date.slice(0, 4),
        points: [point(0, flux, scatter)],
        fields: [["Flux (mJy/beam)", flux], ["Local scatter", `${scatter} mJy/beam`]],
        noCentroid: centroid === "No centroid recovered", cutout: radioCutouts[date],
      };
    },
    initial: rows => rows.findIndex(row => row.noCentroid),
  },
};

export function observationSeries(id, rawRows, aligned = false) {
  const profile = profiles[id];
  if (!profile || !rawRows.length) return null;
  const rows = rawRows.map(raw => profile.row(raw, aligned));
  if (rows.some(r => !Number.isFinite(r.x) || r.points.some(p =>
    !Number.isFinite(p.value) || !Number.isFinite(p.spread) || (profile.log && p.value <= 0)))) return null;
  return {...profile, rows, id, initialRow: Math.max(0, profile.initial(rows))};
}

const extent = numbers => [Math.min(...numbers), Math.max(...numbers)];

function niceStep(span, count) {
  const target = span / count;
  const power = 10 ** Math.floor(Math.log10(target || 1));
  return ([1, 2, 2.5, 5, 10].find(v => v >= target / power) ?? 10) * power;
}

function ticksBetween(lo, hi, step) {
  return Array.from({length: Math.floor(hi / step + 1e-9) - Math.ceil(lo / step - 1e-9) + 1},
    (_, i) => Number(((Math.ceil(lo / step - 1e-9) + i) * step).toPrecision(12)));
}

export function axesFor(data, narrow) {
  const {rows, id} = data;
  let [x0, x1] = extent(rows.map(r => r.x));
  let ticksX;
  if (id === "sun") {
    x0 = Math.floor(x0 / 300) * 300;
    x1 = Math.ceil(x1 / 300) * 300;
    ticksX = ticksBetween(x0, x1, narrow ? 600 : 300).map(t => [t,
      `${Math.floor(t / 3600)}:${String(Math.floor(t / 60) % 60).padStart(2, "0")}`]);
  } else if (id === "quasar") {
    ticksX = ticksBetween(x0, x1, niceStep(x1 - x0, narrow ? 4 : 6)).map(t => [t, number.format(t)]);
  } else {
    const samples = id === "dwarf" ? rows : Array.from({length: narrow ? 4 : 6},
      (_, i) => rows[Math.round(i * (rows.length - 1) / (narrow ? 3 : 5))]);
    ticksX = samples.map(r => [r.x, r.axisTick]);
  }
  const points = rows.flatMap(r => r.points);
  let lo = Math.min(...points.map(p => p.value - p.spread));
  let hi = Math.max(...points.map(p => p.value + p.spread));
  if (data.zero) {
    lo = Math.min(0, lo);
    hi = Math.max(0, hi);
  }
  let ticksY;
  if (data.log) {
    lo = Math.floor(Math.log10(lo));
    hi = Math.ceil(Math.log10(hi));
    ticksY = ticksBetween(lo, hi, 1).map(power => 10 ** power);
    lo = 10 ** lo;
    hi = 10 ** hi;
  } else {
    const step = niceStep(hi - lo, narrow ? 4 : 5);
    lo = Math.floor(lo / step) * step;
    hi = Math.ceil(hi / step) * step;
    if (lo === hi) hi += step;
    ticksY = ticksBetween(lo, hi, step);
  }
  return {x0, x1, lo, hi, ticksX, ticksY};
}
