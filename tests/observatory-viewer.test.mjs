import assert from "node:assert/strict";
import { createHash } from "node:crypto";
import { readFileSync } from "node:fs";
import test from "node:test";
import { axesFor, observationSeries, radioCutouts } from "../observatory/observation-data.js";

const asset = path => new URL(`../observatory/${path}`, import.meta.url);
const catalog = JSON.parse(readFileSync(asset("data/catalog.json"), "utf8"));
const source = Object.fromEntries(catalog.spotlights.map(target => [target.id, target.table.rows]));

test("MAVEN marks the measured days and preserves each channel on a positive log axis", () => {
  const data = observationSeries("mars", source.mars);
  for (const [i, row] of data.rows.entries()) {
    const [date, , selected, ion, electron] = source.mars[i];
    assert.equal(row.x, Date.parse(`${date}T00:00:00Z`) / 86400000);
    assert.equal(row.eventDay, Number(selected) > 0);
    assert.deepEqual(row.points.map(p => p.value), [Number(electron), Number(ion)]);
  }
  assert.equal(source.mars[data.initialRow][0], "2025-08-24");
  const {ticksY} = axesFor(data, false);
  assert.ok(ticksY.every((v, i) => v > 0 && (!i || v / ticksY[i - 1] === 10)));
});

test("STIX times stay in spacecraft UTC, with subsecond precision in the cursor", () => {
  const data = observationSeries("sun", source.sun);
  for (const [i, row] of data.rows.entries()) {
    const [time, rate] = source.sun[i];
    const [hh, mm, ss] = time.split(":").map(Number);
    assert.equal(row.x, hh * 3600 + mm * 60 + ss);
    assert.equal(row.points[0].value, Number(rate));
    assert.equal(row.stamp, `${time} UTC`);
  }
  assert.equal(data.rows[data.initialRow].stamp, "18:27:38.220 UTC");
  assert.equal(data.rows[data.initialRow].points[0].value, 286.35);
});

test("Gaia comparison shifts only B coordinates and retains every published error and native readout", () => {
  const raw = observationSeries("quasar", source.quasar);
  const aligned = observationSeries("quasar", source.quasar, true);
  assert.equal(raw.inverted, true);
  for (const [i, [image, time, magnitude, error]] of source.quasar.entries()) {
    const a = raw.rows[i], b = aligned.rows[i];
    assert.deepEqual(a.fields, b.fields);
    assert.equal(a.x, Number(time));
    assert.equal(b.x, Number(time) + (image === "B" ? 210 : 0));
    assert.equal(a.points[0].value, Number(magnitude));
    assert.equal(b.points[0].value, Number(magnitude) - (image === "B" ? 1.94 : 0));
    assert.equal(b.points[0].spread, Number(error));
  }
});

test("VLASS renders the signed fits, local control scatters and corresponding panel of the source raster", () => {
  const png = readFileSync(asset("data/figures/dwarf.png"));
  assert.equal(createHash("sha256").update(png).digest("hex"),
    "c4a6297aa20953c55199cc97810e61e54219059615ad69ae1d69e63e8a8b19be",
    "Review cutout coordinates against the revised original figure");
  assert.deepEqual([png.readUInt32BE(16), png.readUInt32BE(20)], [1782, 561]);
  const data = observationSeries("dwarf", source.dwarf);
  for (const [i, [date, flux, scatter]] of source.dwarf.entries()) {
    const row = data.rows[i];
    assert.equal(row.points[0].value, Number(flux));
    assert.equal(row.points[0].spread, Number(scatter));
    assert.equal(row.cutout, radioCutouts[date]);
  }
  const selected = data.rows[data.initialRow];
  assert.equal(selected.noCentroid, true);
  assert.equal(selected.points[0].value, -0.150);
  assert.equal(selected.points[0].spread, 0.148);
  assert.ok(axesFor(data, true).ticksY.includes(0));
});
