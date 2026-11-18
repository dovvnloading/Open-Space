import assert from "node:assert/strict";
import test from "node:test";
import {
  address,
  makeSystem,
  orbitRadius,
  parseSystemAddress,
  scaleBar,
  systemSeed,
  universeSize,
} from "../observatory/universe.js";

const { galaxies, systemsPerGalaxy } = universeSize;

test("every displayed system address round-trips and selects a distinct seed", () => {
  const seeds = new Uint32Array(galaxies * systemsPerGalaxy);
  assert.equal(seeds.length, 4_194_304);
  for (let g = 0; g < galaxies; g++) {
    for (let s = 0; s < systemsPerGalaxy; s++) {
      const parsed = parseSystemAddress(address(g, s));
      if (!parsed || parsed.g !== g || parsed.s !== s) assert.fail("Invalid address " + address(g, s));
      seeds[g * systemsPerGalaxy + s] = systemSeed(g, s);
    }
  }
  seeds.sort();
  for (let i = 1; i < seeds.length; i++) {
    if (seeds[i] === seeds[i - 1]) assert.fail("Repeated system seed " + seeds[i]);
  }
  assert.deepEqual(parseSystemAddress(" g1-s1 "), { g: 0, s: 0 });
  for (const invalid of ["G000-S1", "G129-S1", "G1-S00000", "G128–S32769", "G001", "G001-S00001!"]) {
    assert.equal(parseSystemAddress(invalid), null);
  }
});

test("all generated orbits have positive ordered radii and periods consistent with their host mass", () => {
  let minPlanets = Infinity,
    maxPlanets = 0;
  for (let g = 0; g < galaxies; g++) {
    for (let s = 0; s < systemsPerGalaxy; s++) {
      const { star, planets } = makeSystem(g, s);
      minPlanets = Math.min(minPlanets, planets.length);
      maxPlanets = Math.max(maxPlanets, planets.length);
      let previous = 0;
      for (const planet of planets) {
        const massFromPeriod = planet.au ** 3 / (planet.period / 365.25) ** 2;
        if (!(planet.au > previous && planet.period > 0 && Math.abs(massFromPeriod - star.mass) < 1e-12)) {
          assert.fail("Inconsistent orbit at " + address(g, s));
        }
        previous = planet.au;
      }
    }
  }
  assert.equal(minPlanets, 2);
  assert.equal(maxPlanets, 8);
});

test("the opening and final system keep repeatable orbital solutions", () => {
  const last = makeSystem(127, 32767);
  const opening = makeSystem(0, 228);
  assert.equal(opening.star.mass, 1.94);
  assert.deepEqual(
    opening.planets.map((p) => p.au.toFixed(4)),
    ["0.4620", "0.7138", "1.3439", "2.1004", "3.7644", "7.6021", "13.5969"],
  );
  assert.deepEqual(
    last.planets.map((p) => p.au.toFixed(4)),
    ["0.3468", "0.5820"],
  );
  assert.deepEqual(makeSystem(127, 32767), last);
});

test("AU spacing is linear in orbital radius and the ruler spans its labeled distance", () => {
  const { planets } = makeSystem(0, 228);
  const last = planets.length - 1;
  const outer = orbitRadius(planets, last, true);
  assert.equal(outer, orbitRadius(planets, last, false));
  for (const [i, planet] of planets.entries()) {
    assert.ok(Math.abs(orbitRadius(planets, i, true) / outer - planet.au / planets[last].au) < 1e-12);
  }
  for (const unitsPerPixel of [0.0005, 0.01, 0.2, 1, 11, 30, 91.2, 1_000, 90_000]) {
    const ruler = scaleBar(unitsPerPixel);
    assert.ok(ruler.pixels >= 28 && ruler.pixels <= 72);
    assert.ok(Math.abs(ruler.pixels * unitsPerPixel - ruler.units) <= ruler.units * 1e-12);
  }
});
