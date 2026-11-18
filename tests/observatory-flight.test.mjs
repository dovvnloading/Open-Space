import assert from "node:assert/strict";
import test from "node:test";
import { pointingArc } from "../observatory/flight-path.js";

const point = (raDeg, decDeg) => ({ raDeg, decDeg });
test("sky courses show angular separation on the sphere", () => {
  const polar = pointingArc(point(30, 80), point(210, 80));
  assert.ok(Math.abs(polar.degrees - 20) < 1e-8);
  assert.match(polar.path, /^M/);
  const seam = pointingArc(point(352, 0), point(4, 0));
  assert.ok(Math.abs(seam.degrees - 12) < 1e-8);
  assert.equal((seam.path.match(/M/g) || []).length, 2);
});
test("coincident and antipodal pointings produce finite distances", () => {
  for (const [ra, expected] of [
    [0, 0],
    [180, 180],
  ]) {
    const slew = pointingArc(point(0, 0), point(ra, 0));
    assert.equal(slew.degrees, expected);
    assert.equal(slew.path, "");
  }
});
