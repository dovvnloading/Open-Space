/* A deterministic, illustrative sky. The input seed and integer address define
   every galaxy, star and planet independently of the camera or frame rate. */
const TAU = Math.PI * 2;
const GALAXIES = 128;
const SYSTEMS = 32768;
export const universeSize = Object.freeze({ galaxies: GALAXIES, systemsPerGalaxy: SYSTEMS });
const MAP_EXTENT = 1.29;
const SECTOR_MLY_PER_UNIT = 22;
const DAYS_PER_SECOND = 12;
const OUTER_ORBIT_RADIUS = 0.98;
const GALAXY_WORDS = [
  "Vesper",
  "Lyric",
  "Aster",
  "Sable",
  "Halcyon",
  "Lumen",
  "Ember",
  "Aurelia",
  "Cinder",
  "Noctis",
  "Velum",
  "Opal",
  "Argent",
  "Saffron",
  "Cerulean",
  "Serein",
];
const GALAXY_TYPES = ["Spiral", "Reach", "Crown", "Arc", "Drift", "Veil", "Helix", "Sea"];
const STAR_TYPES = [
  { name: "M dwarf", color: "#ff9b73", mass: 0.45 },
  { name: "K dwarf", color: "#ffbc80", mass: 0.76 },
  { name: "G dwarf", color: "#ffe3a3", mass: 0.99 },
  { name: "F dwarf", color: "#fcf0d6", mass: 1.27 },
  { name: "A star", color: "#d6ecff", mass: 1.94 },
  { name: "B star", color: "#a4cbff", mass: 4.38 },
];
const PLANET_TYPES = [
  { name: "Rocky planet", colors: ["#b68767", "#ead0a0", "#424a52"] },
  { name: "Ice world", colors: ["#7fb7c6", "#d7eff4", "#405776"] },
  { name: "Gas giant", colors: ["#cf9868", "#f2d9a9", "#664b5d"] },
  { name: "Ocean world", colors: ["#5386ae", "#a0d9dc", "#29375d"] },
];
export const atlasDestinations = [
  { g: 0, s: 228, level: "galaxy" },
  { g: 22, s: 1912, level: "galaxy" },
  { g: 34, s: 416, level: "galaxy" },
  { g: 0, s: 228, planet: 4, level: "world" },
];
const ROMAN = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"];
const limit = (n, min, max) => Math.max(min, Math.min(max, n));
const pad = (n, size) => String(n + 1).padStart(size, "0");
const galaxyName = (id) => GALAXY_WORDS[id % 16] + " " + GALAXY_TYPES[Math.floor(id / 16)];
const starName = (id) => "System " + pad(id, 5);
export const address = (g, s) => "G" + pad(g, 3) + (s === undefined ? "" : "–S" + pad(s, 5));

export function parseSystemAddress(input) {
  const match = /^g(\d{1,3})[\s–-]+s(\d{1,5})$/i.exec(input.trim());
  if (!match) return null;
  const g = Number(match[1]) - 1,
    s = Number(match[2]) - 1;
  return g >= 0 && g < GALAXIES && s >= 0 && s < SYSTEMS ? { g, s } : null;
}

const openingLocation = Object.freeze({ level: "galaxy", g: 0, s: 228, p: 0 });
const levels = new Set(["sector", "galaxy", "system", "world"]);

export function universeLocationHash({ level, g, s, p }) {
  const at = address(g, s).replace("–", "-");
  const params = new URLSearchParams({ view: level, at });
  if (level === "system" || level === "world") params.set("planet", String(p + 1));
  return "#universe?" + params.toString();
}

export function parseUniverseLocation(search = "") {
  if (!search) return { ...openingLocation };
  const params = new URLSearchParams(search);
  const allowed = new Set(["view", "at", "planet"]);
  for (const name of params.keys()) {
    if (!allowed.delete(name)) return null; // Reject repeated and unknown fields.
  }
  const level = params.get("view");
  const point = parseSystemAddress(params.get("at") ?? "");
  if (!levels.has(level) || !point) return null;
  const ordinal = params.get("planet");
  if (ordinal !== null && (level === "sector" || level === "galaxy" || !/^[1-8]$/.test(ordinal))) return null;
  if (level === "world" && ordinal === null) return null;
  const p = ordinal === null ? 0 : Number(ordinal) - 1;
  if (p >= makeSystem(point.g, point.s).planets.length) return null;
  return { level, ...point, p };
}

function seedNumber(a, b = 0) {
  let h = Math.imul(a + 29, 0x85ebca6b) ^ Math.imul(b + 317, 0xc2b2ae35);
  h = Math.imul(h ^ (h >>> 16), 0x27d4eb2d);
  return (h ^ (h >>> 15)) >>> 0;
}

function random(seed) {
  return () => {
    seed = (seed + 0x6d2b79f5) | 0;
    let h = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    h ^= h + Math.imul(h ^ (h >>> 7), 61 | h);
    return ((h ^ (h >>> 14)) >>> 0) / 4294967296;
  };
}

function canvasSurface(size) {
  const surface = document.createElement("canvas");
  surface.width = surface.height = size;
  return surface;
}

function glow(ctx, x, y, radius, rgb, strength = 1) {
  const grad = ctx.createRadialGradient(x, y, 0, x, y, radius);
  grad.addColorStop(0, "rgba(" + rgb + "," + strength + ")");
  grad.addColorStop(0.16, "rgba(" + rgb + "," + strength * 0.62 + ")");
  grad.addColorStop(0.55, "rgba(" + rgb + "," + strength * 0.15 + ")");
  grad.addColorStop(1, "rgba(" + rgb + ",0)");
  ctx.fillStyle = grad;
  ctx.fillRect(x - radius, y - radius, 2 * radius, 2 * radius);
}

const families = {
  spiral: "Spiral galaxy",
  barred: "Barred spiral",
  edge: "Edge-on spiral",
  ring: "Ring galaxy",
  elliptical: "Elliptical galaxy",
};
export function galaxyProfile(id) {
  const rng = random(seedNumber(id));
  const profile = {
    arms: 2 + Math.floor(rng() * 3),
    winding: 3.7 + rng() * 2.7,
    tilt: 0.49 + 0.27 * rng(),
    angle: -0.4 + rng() * 0.8,
    diameter: (81 + Math.floor(rng() * 56)) * 1000,
    hue: rng(),
    phase: rng() * TAU,
  };
  // The first group has traditional spiral names. Other groups also contain
  // ring and elliptical galaxies. Shape and selectable star density agree.
  const family =
    id < 16
      ? ["spiral", "barred", "edge"][id % 3]
      : ["spiral", "barred", "ring", "edge", "elliptical"][id % 5];
  if (family === "edge") {
    profile.tilt = 0.19 + rng() * 0.08;
    profile.angle = -0.27 + rng() * 0.5;
  }
  if (family === "barred") profile.arms = 2;
  if (family === "ring" || family === "elliptical") {
    profile.arms = 0;
    profile.tilt = 0.6 + rng() * 0.24;
  }
  return { ...profile, family, label: families[family] };
}

function pointAt(radius, angle, thickness, config) {
  const x = radius * Math.cos(angle);
  const y = radius * Math.sin(angle) * config.tilt + thickness;
  return [
    x * Math.cos(config.angle) - y * Math.sin(config.angle),
    x * Math.sin(config.angle) + y * Math.cos(config.angle),
  ];
}

function armAngle(r, arm, config) {
  return (
    (arm * TAU) / config.arms +
    config.phase +
    config.winding * Math.pow(r + 0.015, 0.73) +
    0.12 * Math.sin(r * 11 + arm * 8) +
    0.046 * Math.cos(r * 28 + arm * 2)
  );
}

// Sample a shared density field in the plane of the disk. Broad tidal spurs,
// dusty inner edges and small emission knots break up the logarithmic arms.
function makeDisk(id, config, size = 680) {
  const tile = canvasSurface(size);
  const ctx = tile.getContext("2d");
  const pixels = ctx.createImageData(tile.width, tile.height);
  const rng = random(seedNumber(id, 9650));
  const noise = Float32Array.from({ length: 65536 }, () => rng());
  function sample(x, y) {
    const ix = Math.floor(x),
      iy = Math.floor(y);
    let fx = x - ix,
      fy = y - iy;
    fx = fx * fx * (3 - 2 * fx);
    fy = fy * fy * (3 - 2 * fy);
    const row = (iy & 255) << 8,
      next = ((iy + 1) & 255) << 8;
    const a = ix & 255,
      b = (ix + 1) & 255;
    return (
      (noise[row + a] * (1 - fx) + noise[row + b] * fx) * (1 - fy) +
      (noise[next + a] * (1 - fx) + noise[next + b] * fx) * fy
    );
  }
  const c = Math.cos(config.angle),
    s = Math.sin(config.angle);
  const unit = tile.width / (2 * MAP_EXTENT);
  const barAxis = config.family === "barred" ? armAngle(0.32, 0, config) : 0;
  const barCos = Math.cos(barAxis),
    barSin = Math.sin(barAxis);
  const warm = config.hue > 0.63 ? 1.17 : 1;
  for (let j = 0; j < tile.height; j++) {
    const py = (j - tile.height / 2) / unit;
    for (let i = 0; i < tile.width; i++) {
      const px = (i - tile.width / 2) / unit;
      const x = px * c + py * s,
        y = (py * c - px * s) / config.tilt;
      const r = Math.hypot(x, y);
      if (r > 1.17) continue;
      const broad = sample(x * 7 + 36, y * 7 + 78);
      const cloud =
        0.57 * sample(x * 20 + broad * 4, y * 20 + broad * 3) +
        0.27 * sample(x * 57, y * 57) +
        0.16 * sample(x * 140, y * 140);
      const o = 4 * (j * tile.width + i);
      const fadeEdge = Math.min(1, Math.max(0, (1.17 - r) * 8));
      if (!config.arms) {
        const asymmetry = 0.028 * Math.sin(3 * Math.atan2(y, x) + config.phase) + (broad - 0.5) * 0.055;
        const ring = Math.exp(-(((r - 0.73 + asymmetry) / 0.098) ** 2));
        const nucleus = Math.exp(-r * (config.family === "ring" ? 12 : 5.2));
        const envelope = 52 * Math.exp(-2.5 * r) * (0.87 + 0.23 * cloud);
        const clumps = ring * (38 + 250 * Math.max(0, cloud - 0.21) ** 1.6);
        const emission = config.family === "ring" ? clumps : 0;
        const red = config.family === "ring" ? 0.48 : 1.18;
        pixels.data[o] = red * envelope + nucleus * 255 + emission * 0.68;
        pixels.data[o + 1] = envelope + nucleus * 177 + emission * 1.27;
        pixels.data[o + 2] = envelope * 1.13 + nucleus * 104 + emission * 1.63;
        pixels.data[o + 3] = 255 * fadeEdge;
        continue;
      }
      const phase = Math.atan2(y, x) - config.phase - config.winding * Math.pow(r + 0.015, 0.73);
      const angle = phase * config.arms + (broad - 0.5) * 1.31 + (cloud - 0.5) * 0.7;
      const arm = 0.5 + 0.5 * Math.cos(angle);
      const ridge = arm ** 5.5,
        spur = (0.5 + 0.5 * Math.cos(angle * 1.99 + 2.3 + broad * 4)) ** 6;
      const fade = Math.max(0, 1 - (r / 1.15) ** 3) * Math.exp(-r * 1.18);
      const light = Math.max(0.04, (cloud - 0.23) * 2.3);
      const lane = 1 - 0.83 * (0.5 + 0.5 * Math.cos(angle + 0.48)) ** 26 * Math.min(1, r * 8);
      const disk = (18 + cloud * 16 + 26 * spur * light) * fade * lane;
      const gas = (20 + 162 * light * light) * ridge * fade * lane;
      const core = Math.exp(-r * 8) * 91 * (0.65 + cloud);
      const knot = Math.max(0, (cloud - 0.62) * 6) * ridge * fade * 47;
      const along = x * barCos + y * barSin;
      const across = y * barCos - x * barSin;
      const bar =
        config.family === "barred" ? 106 * Math.exp(-((along / 0.31) ** 4 + (across / 0.055) ** 2)) : 0;
      pixels.data[o] = (bar * 1.15 + disk * 1.07 + gas * 0.72 + core * 1.6 + knot * 2.25) * warm;
      pixels.data[o + 1] = bar * 0.89 + disk * 1.05 + gas * 1.02 + core * 1.1 + knot * 0.53;
      pixels.data[o + 2] = bar * 0.62 + disk * 1.46 + gas * 1.48 + core * 0.6 + knot * 1.28;
      pixels.data[o + 3] = Math.min(255, 255 * (1.17 - r) * 12);
    }
  }
  ctx.putImageData(pixels, 0, 0);
  return tile;
}

let firstGalaxy;
function makeGalaxy(id, paint = true) {
  if (paint && id === 0 && firstGalaxy) return firstGalaxy;
  const rng = random(seedNumber(id, 270));
  const config = galaxyProfile(id);
  const position = new Float32Array(SYSTEMS * 2);
  const magnitude = new Float32Array(SYSTEMS);
  const tint = new Uint8Array(SYSTEMS);
  const inks = ["#fff6d8", "#f6ce9c", "#94aee1", "#b0e5ff", "#d0d2ff", "#e5edff", "#f0a59f"];
  const texture = paint ? canvasSurface(1440) : null;
  const ctx = texture?.getContext("2d");
  const half = 720,
    scale = half / MAP_EXTENT;
  if (ctx) {
    ctx.drawImage(makeDisk(id, config), 0, 0, texture.width, texture.height);
    ctx.globalCompositeOperation = "screen";
    ctx.save();
    ctx.translate(half, half);
    ctx.rotate(config.angle);
    ctx.scale(1, config.tilt);
    glow(ctx, 0, 0, scale * 1.12, "67,76,120", 0.38);
    glow(ctx, 0, 0, scale * 0.53, "146,108,99", 0.7);
    glow(ctx, 0, 0, scale * 0.3, "241,168,98", 0.84);
    ctx.restore();
  }
  if (ctx && config.arms) {
    const cloudRng = random(seedNumber(id, 880));
    const cloudColors = ["139,150,199", "79,123,174", "180,133,129", "132,113,169"];
    const cloudSprites = cloudColors.map((rgb) => {
      const tile = canvasSurface(96);
      const brush = tile.getContext("2d");
      glow(brush, 48, 48, 48, rgb);
      return tile;
    });
    for (let i = 0; i < 4700; i++) {
      const r = 0.017 + cloudRng() * 1.07;
      const arm = i % config.arms;
      const blur = cloudRng() + cloudRng() + cloudRng() - 1.5;
      const a = armAngle(r, arm, config) + blur * (0.11 + 0.11 * r);
      const [x, y] = pointAt(r, a, (cloudRng() - 0.5) * 0.024, config);
      const size = scale * (0.016 + cloudRng() * 0.087) * (1.35 - r * 0.48);
      const tone = r < 0.28 && cloudRng() > 0.4 ? 2 : Math.floor(cloudRng() * 4);
      ctx.globalAlpha = (0.025 + cloudRng() * 0.1) * Math.max(0, 1.12 - r) ** 0.61;
      ctx.drawImage(cloudSprites[tone], half + x * scale - size, half + y * scale - size, size * 2, size * 2);
    }
    const dust = canvasSurface(96);
    ctx.globalCompositeOperation = "source-over";
    glow(dust.getContext("2d"), 48, 48, 48, "8,10,17");
    for (let i = 0; i < 2500; i++) {
      const r = 0.026 + cloudRng() * 0.94;
      const arm = i % config.arms;
      const a = armAngle(r, arm, config) - 0.12 + (cloudRng() - 0.5) * 0.12;
      const [x, y] = pointAt(r, a, (cloudRng() - 0.5) * 0.008, config);
      const w = scale * (0.002 + cloudRng() * 0.031);
      ctx.globalAlpha = 0.035 + cloudRng() * 0.16;
      ctx.drawImage(dust, half + x * scale - w, half + y * scale - w, w * 2, w * 2);
    }
    ctx.globalAlpha = 1;
  }
  const haze = ctx ? canvasSurface(740) : null;
  if (haze) haze.getContext("2d").drawImage(texture, 0, 0, haze.width, haze.height);
  for (let i = 0; i < SYSTEMS; i++) {
    const bulge = rng() < 0.23;
    let r = bulge ? 0.42 * Math.pow(rng(), 1.77) : 0.022 + Math.pow(rng(), 1.12) * 1.065;
    let a =
      bulge || !config.arms || rng() < 0.2
        ? rng() * TAU
        : armAngle(r, i % config.arms, config) + (rng() + rng() + rng() - 1.5) * (0.25 + r * 0.24);
    if (config.family === "ring" && !bulge && rng() < 0.82) {
      r = 0.73 - 0.028 * Math.sin(3 * a + config.phase) + (rng() + rng() + rng() - 1.5) * 0.105;
    } else if (config.family === "elliptical") {
      r = 1.05 * rng() ** 1.85;
    } else if (config.family === "barred" && r < 0.34 && rng() < 0.76) {
      a = armAngle(0.32, 0, config) + (i % 2) * Math.PI + (rng() + rng() - 1) * 0.28;
    }
    const [x, y] = pointAt(r, a, (rng() + rng() - 1) * (bulge ? 0.082 : 0.024), config);
    position[i * 2] = x;
    position[i * 2 + 1] = y;
    magnitude[i] = 0.22 + Math.pow(rng(), 6) * 1.55;
    tint[i] =
      (r < 0.24 || config.family === "elliptical") && rng() < 0.72
        ? Math.floor(rng() * 2)
        : Math.floor(rng() * inks.length);
    if (ctx) {
      ctx.fillStyle = inks[tint[i]];
      ctx.globalAlpha = (0.18 + magnitude[i] * 0.34) * Math.min(1, (1.14 - r) * 3.9);
      const w = 0.42 + magnitude[i];
      ctx.fillRect(half + x * scale, half + y * scale, w, w);
    }
  }
  if (ctx) {
    ctx.globalAlpha = 1;
    ctx.save();
    ctx.translate(half, half);
    ctx.rotate(config.angle);
    ctx.scale(1, config.tilt + 0.09);
    ctx.globalCompositeOperation = "screen";
    glow(ctx, 0, 0, scale * 0.26, "242,190,127", 0.88);
    glow(ctx, 0, 0, scale * 0.102, "255,222,167", 0.98);
    glow(ctx, 0, 0, scale * 0.044, "255,246,212", 0.97);
    ctx.restore();
  }
  const result = { config, position, magnitude, tint, texture, haze, inks };
  if (paint && id === 0) firstGalaxy = result;
  return result;
}

export function systemSeed(g, s) {
  // Flatten the two coordinates before applying the invertible 32-bit mix.
  // This retains the original G001 sequence and gives every address its own seed.
  return seedNumber(742, g * SYSTEMS + s);
}

export function makeSystem(g, s) {
  const rng = random(systemSeed(g, s));
  const star = STAR_TYPES[Math.floor(rng() * STAR_TYPES.length)];
  const count = 2 + Math.floor(rng() * 7);
  let distance = 0.12 + rng() * 0.22;
  const planets = Array.from({ length: count }, (_, i) => {
    distance *= 1.4 + rng() * 0.65;
    let type = Math.floor(rng() * 4);
    if (!i) type = 0;
    return {
      type: PLANET_TYPES[type],
      au: distance,
      phase: rng() * TAU,
      period: 365.25 * Math.sqrt(distance ** 3 / star.mass),
      ring: type === 2 && rng() > 0.42,
      size: (type === 2 ? 0.025 : 0.011) + rng() * 0.008,
    };
  });
  return { star, planets };
}

const ringTilt = -0.29;
const ringAspect = 0.29;
const ringDepth = Math.sqrt(1 - ringAspect ** 2);
const ringCos = Math.cos(ringTilt),
  ringSin = Math.sin(ringTilt);
const light = [-0.82, -0.33, Math.sqrt(1 - 0.82 ** 2 - 0.33 ** 2)];
const lightAcross = -light[0] * ringSin + light[1] * ringCos;
const lightNormal = light[2] * ringAspect - lightAcross * ringDepth;
const ringInner = 1.18,
  ringOuter = 1.79;

function ringDensity(radius) {
  if (radius < ringInner || radius > ringOuter || (radius > 1.565 && radius < 1.594)) return 0;
  const edge = Math.min(1, (radius - ringInner) * 9, (ringOuter - radius) * 16);
  const band = radius < 1.32 ? 0.44 : radius < 1.565 ? 0.92 : 0.69;
  return band * edge * (0.87 + 0.13 * Math.cos(radius * 268));
}

function paintRings(ctx, cx, cy, scale, front) {
  ctx.save();
  ctx.translate(cx, cy);
  ctx.rotate(ringTilt);
  const steps = 72;
  for (let j = 0; j < steps; j++) {
    const radius = ringInner + ((j + 0.5) * (ringOuter - ringInner)) / steps;
    const density = ringDensity(radius);
    if (!density) continue;
    const shade = Math.round(158 + density * 76);
    const ringColor = "rgba(" + shade + "," + (shade - 10) + "," + (shade - 28) + "," + density + ")";
    const shadowColor = "rgba(50,52,67," + density + ")";
    // Segment only where a ray from the ring to the illuminating star meets the globe.
    const first = front ? 0 : Math.PI;
    let start = first,
      previous;
    for (let k = 0; k <= 96; k++) {
      const angle = first + (k * Math.PI) / 96;
      const x = radius * Math.cos(angle),
        y = radius * Math.sin(angle) * ringAspect;
      const z = radius * Math.sin(angle) * ringDepth;
      const dot =
        (x * ringCos - y * ringSin) * light[0] + (y * ringCos + x * ringSin) * light[1] + z * light[2];
      const shadow = dot < 0 && radius * radius - dot * dot < 1;
      if (k && (previous !== shadow || k === 96)) {
        ctx.strokeStyle = previous ? shadowColor : ringColor;
        ctx.lineWidth = (scale * (ringOuter - ringInner)) / steps + 0.24;
        ctx.beginPath();
        ctx.ellipse(0, 0, radius * scale, radius * scale * ringAspect, 0, start, angle);
        ctx.stroke();
        start = angle;
      }
      previous = shadow;
    }
  }
  ctx.restore();
}

function drawWorldPortrait(ctx, cx, cy, radius, planet, globe) {
  if (planet.ring) paintRings(ctx, cx, cy, radius, false);
  ctx.drawImage(globe, cx - radius, cy - radius, radius * 2, radius * 2);
  if (planet.ring) paintRings(ctx, cx, cy, radius, true);
}

// Bucket projected crater centers so each sample touches only nearby impacts.
// Distances and incidence are evaluated in 3D on the unit sphere.
function craterField(rng) {
  const count = 24;
  const bins = Array.from({ length: count * count }, () => []);
  const bin = (x) => limit(Math.floor(((x + 1) * count) / 2), 0, count - 1);
  for (let i = 0; i < 440; i++) {
    const z = rng(),
      angle = rng() * TAU,
      radial = Math.sqrt(1 - z * z);
    const x = radial * Math.cos(angle),
      y = radial * Math.sin(angle);
    const r = 0.007 + 0.118 * rng() ** 5;
    const crater = { x, y, z, r, age: 0.55 + 0.45 * rng() };
    for (let iy = bin(y - r * 1.3); iy <= bin(y + r * 1.3); iy++) {
      for (let ix = bin(x - r * 1.3); ix <= bin(x + r * 1.3); ix++) bins[iy * count + ix].push(crater);
    }
  }
  return (x, y, z) => {
    let relief = 0;
    for (const crater of bins[bin(y) * count + bin(x)]) {
      const dx = x - crater.x,
        dy = y - crater.y,
        dz = z - crater.z;
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz) / crater.r;
      if (d > 1.3) continue;
      const sunward = (dx * light[0] + dy * light[1] + dz * light[2]) / crater.r;
      const rim = Math.exp(-(((d - 0.88) * 9) ** 2));
      const slope = Math.exp(-(((d - 0.66) * 3.7) ** 2));
      relief += crater.age * (rim * 34 + sunward * slope * 83 - 15 * Math.max(0, 1 - d * d));
    }
    return relief;
  };
}

// Evaluate the surface on a sphere so details follow the limb. Each world has
// a fixed terrain field independent of the screen resolution and viewpoint.
function makeWorld(g, s, index, planet, size = 640) {
  const image = canvasSurface(size),
    ctx = image.getContext("2d");
  const pixels = ctx.createImageData(image.width, image.height);
  const rng = random(seedNumber(936, (g * SYSTEMS + s) * ROMAN.length + index));
  const table = Float32Array.from({ length: 32768 }, () => rng());
  const offset = 20 * rng(),
    yaw = rng() * TAU,
    axis = 0.23 + rng() * 0.52;
  const cy = Math.cos(yaw),
    sy = Math.sin(yaw),
    ca = Math.cos(axis),
    sa = Math.sin(axis);
  const contrast = 0.8 + 0.35 * rng(),
    hue = rng();
  function noise(x, y, z) {
    const i = Math.floor(x),
      j = Math.floor(y),
      k = Math.floor(z);
    let a = x - i,
      b = y - j,
      c = z - k;
    a = a * a * (3 - 2 * a);
    b = b * b * (3 - 2 * b);
    c = c * c * (3 - 2 * c);
    const ix = i & 31,
      ix1 = (i + 1) & 31,
      row = (j & 31) << 5,
      row1 = ((j + 1) & 31) << 5;
    const z0 = (k & 31) << 10,
      z1 = ((k + 1) & 31) << 10;
    return (
      (1 - c) *
        ((1 - b) * (table[z0 + row + ix] * (1 - a) + table[z0 + row + ix1] * a) +
          b * (table[z0 + row1 + ix] * (1 - a) + table[z0 + row1 + ix1] * a)) +
      c *
        ((1 - b) * (table[z1 + row + ix] * (1 - a) + table[z1 + row + ix1] * a) +
          b * (table[z1 + row1 + ix] * (1 - a) + table[z1 + row1 + ix1] * a))
    );
  }
  const type = PLANET_TYPES.indexOf(planet.type);
  const half = image.width / 2;
  const craters = type === 0 ? craterField(rng) : null;
  for (let j = 0; j < image.height; j++) {
    const y = (j + 0.5 - half) / (half - 1);
    for (let i = 0; i < image.width; i++) {
      const x = (i + 0.5 - half) / (half - 1);
      const rr = x * x + y * y;
      if (rr >= 1) continue;
      const z = Math.sqrt(1 - rr);
      const ux = x * cy + z * sy,
        uz = z * cy - x * sy;
      const v = y * ca - ux * sa,
        u = ux * ca + y * sa;
      const low = noise(u * 2.7 + offset, v * 2.7, uz * 2.7);
      const mid = noise(u * 7 + 2 * low + offset, v * 7 + low, uz * 7);
      const fine = noise(u * 21 + mid * 3, v * 21 + mid, uz * 21);
      const grain = noise(u * 67, v * 67 + fine, uz * 67);
      const rock = 0.49 * low + 0.26 * mid + 0.17 * fine + 0.08 * grain;
      let r, b, green;
      if (type === 2) {
        const latitude = Math.asin(limit(v, -1, 1));
        const wave = Math.sin(latitude * 24 + (mid - 0.5) * 1.6 + 0.37 * Math.sin(latitude * 91 + fine * 3));
        const swirl = 0.5 + 0.5 * Math.sin(latitude * 71 + (low - 0.5) * 7 + (mid - 0.5) * 8 + fine * 2);
        const belt = 0.49 + wave * 0.23 + swirl * 0.12 + rock * 0.25;
        r = 132 + 111 * belt;
        green = 74 + 130 * belt;
        b = 49 + (103 + hue * 27) * belt;
      } else if (type === 3) {
        const land = limit((rock - 0.54) * 42, 0, 1);
        const cloud = limit((0.55 * mid + 0.3 * fine + grain * 0.15 - 0.58) * 4.5, 0, 0.9);
        r = (26 + 104 * land) * (1 - cloud) + 222 * cloud;
        green = (66 + 69 * land) * (1 - cloud) + 235 * cloud;
        b = (109 + 39 * low - 63 * land) * (1 - cloud) + 244 * cloud;
      } else if (type === 1) {
        const cracks = Math.max(0, 1 - Math.abs(fine * 0.22 + mid * 0.35 + low * 0.43 - 0.5) * 44);
        r = 106 + rock * 123 - cracks * 47;
        green = 150 + rock * 95 - cracks * 32;
        b = 176 + rock * 71 - cracks * 13;
      } else {
        const ridges = (1 - Math.abs(fine - 0.5) * 2) ** 6;
        const strata = Math.max(0, Math.sin(rock * 33 + hue * 5)) * 17;
        const relief = craters(x, y, z);
        const elevation = relief + (rock - 0.45) * 155 * contrast + fine * 25 + grain * 12;
        r = 142 + elevation + ridges * 20 + strata;
        green = 94 + elevation * 0.87 + ridges * 13 + strata;
        b = 74 + elevation * 0.7 + hue * 18 + ridges * 7 + strata;
      }
      const incidence = x * light[0] + y * light[1] + z * light[2];
      const scatter = type === 3 ? 0.15 : 0.08;
      let brightness = 0.028 + 1.12 * Math.max(0, incidence) ** 0.64 + scatter * Math.max(0, z) ** 2;
      if (planet.ring && incidence > 0) {
        const across = -x * ringSin + y * ringCos;
        const t = (across * ringDepth - z * ringAspect) / lightNormal;
        if (t > 0) {
          const ringRadius = Math.hypot(x + t * light[0], y + t * light[1], z + t * light[2]);
          brightness *= 1 - 0.67 * ringDensity(ringRadius);
        }
      }
      const limb = Math.pow(1 - z, 5) * Math.max(0, incidence + 0.44) * (type === 0 ? 0.3 : 0.86);
      const k = 4 * (j * image.width + i);
      pixels.data[k] = r * brightness + limb * 114;
      pixels.data[k + 1] = green * brightness + limb * 151;
      pixels.data[k + 2] = b * brightness + limb * 223;
      pixels.data[k + 3] = Math.min(255, (1 - Math.sqrt(rr)) * half * 255);
    }
  }
  ctx.putImageData(pixels, 0, 0);
  return image;
}

const sector = (() => {
  const rng = random(759831);
  return Array.from({ length: GALAXIES }, (_, i) => {
    const angle = rng() * TAU,
      r = Math.sqrt(rng()) * 1.08;
    return {
      x: Math.cos(angle) * r,
      y: Math.sin(angle) * r * 0.78,
      tilt: rng() * Math.PI,
      color: i % 4,
      size: 2.8 + rng() * 3.9,
    };
  });
})();

// A distant, illustrated dust field gives each address a quiet depth cue.
function makeBackdrop(id) {
  const image = canvasSurface(420);
  image.height = 246;
  const ctx = image.getContext("2d");
  const pixels = ctx.createImageData(image.width, image.height);
  const rng = random(seedNumber(id, 914));
  const grid = Float32Array.from({ length: 4096 }, () => rng());
  function noise(x, y) {
    const ix = Math.floor(x),
      iy = Math.floor(y);
    const tx = x - ix,
      ty = y - iy;
    const u = tx * tx * (3 - 2 * tx),
      v = ty * ty * (3 - 2 * ty);
    const a = (iy & 63) << 6,
      b = ((iy + 1) & 63) << 6;
    const left = ix & 63,
      right = (ix + 1) & 63;
    return (
      (grid[a + left] * (1 - u) + grid[a + right] * u) * (1 - v) +
      (grid[b + left] * (1 - u) + grid[b + right] * u) * v
    );
  }
  for (let y = 0; y < image.height; y++) {
    const v = y / image.height;
    for (let x = 0; x < image.width; x++) {
      const u = x / image.width;
      const broad = noise(u * 4.6 + 5, v * 3 + 16);
      const detail =
        noise(u * 24 + broad * 4, v * 15 + broad * 5) * 0.65 + noise(u * 58 + broad * 4, v * 38) * 0.35;
      const band = Math.exp(-(((v - 0.78 + u * 0.59 + (broad - 0.5) * 0.32) / 0.12) ** 2));
      const filament = band * Math.max(0, 0.65 * broad + 0.7 * detail - 0.46) ** 1.2;
      const edge = Math.max(0, u - 0.35);
      const o = (y * image.width + x) * 4;
      pixels.data[o] = 4 + broad * 3 + filament * (19 + edge * 25);
      pixels.data[o + 1] = 7 + broad * 4 + filament * (31 - edge * 13);
      pixels.data[o + 2] = 13 + broad * 6 + filament * 45;
      pixels.data[o + 3] = 255;
    }
  }
  ctx.putImageData(pixels, 0, 0);
  return image;
}

const sky = (() => {
  const rng = random(43153);
  return Array.from({ length: 510 }, () => [rng(), rng(), rng(), rng()]);
})();

/** Draw the opening galaxy in a dashboard canvas with an explicit CSS size. */
export function drawUniversePreview(canvas) {
  const ctx = canvas.getContext("2d", { alpha: false });
  if (!ctx) return () => {};
  let frame = 0,
    visible = true;
  function paint() {
    frame = 0;
    if (!visible || !canvas.isConnected) return;
    const { width, height } = canvas.getBoundingClientRect();
    if (!width || !height) return;
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.fillStyle = "#05080d";
    ctx.fillRect(0, 0, width, height);
    for (const [x, y, magnitude] of sky) {
      ctx.globalAlpha = 0.14 + magnitude * 0.57;
      ctx.fillStyle = magnitude > 0.82 ? "#f1ddc7" : "#bfd4f1";
      const size = 0.4 + magnitude * 0.78;
      ctx.fillRect(x * width, y * height, size, size);
    }
    ctx.globalAlpha = 1;
    const galaxy = makeGalaxy(0);
    const extent = MAP_EXTENT * Math.min(width * 0.48, height * 0.63);
    ctx.globalCompositeOperation = "screen";
    ctx.drawImage(galaxy.texture, width * 0.5 - extent, height * 0.51 - extent, extent * 2, extent * 2);
    ctx.globalCompositeOperation = "source-over";
  }
  function schedule() {
    if (!frame && visible) frame = requestAnimationFrame(paint);
  }
  const resize = typeof ResizeObserver === "function" ? new ResizeObserver(schedule) : null;
  const intersection =
    typeof IntersectionObserver === "function"
      ? new IntersectionObserver(
          (entries) => {
            visible = entries[0].isIntersecting;
            if (visible) schedule();
          },
          { rootMargin: "200px" },
        )
      : null;
  resize?.observe(canvas);
  if (intersection) {
    visible = false;
    intersection.observe(canvas);
  } else schedule();
  if (!resize) window.addEventListener("resize", schedule);
  return () => {
    visible = false;
    resize?.disconnect();
    intersection?.disconnect();
    cancelAnimationFrame(frame);
    if (!resize) window.removeEventListener("resize", schedule);
  };
}

const icons = {
  close: '<path d="m6 6 12 12M18 6 6 18"/>',
  back: '<path d="m14 6-6 6 6 6"/>',
  forward: '<path d="m10 6 6 6-6 6"/>',
  zoomIn: '<path d="M12 5v14M5 12h14"/>',
  zoomOut: '<path d="M5 12h14"/>',
  center: '<circle cx="12" cy="12" r="5"/><path d="M12 2v5m0 10v5M2 12h5m10 0h5"/>',
  help: '<circle cx="12" cy="12" r="9"/><path d="M9.8 9a2.4 2.4 0 1 1 4.2 1.7c-1.1.8-2 1.4-2 3M12 16.5h.01"/>',
  explore: '<path d="m3 12 17-8-8 17-1.5-7.5L3 12Z"/>',
  link: '<path d="M10 13a4 4 0 0 0 5.6 0l3.5-3.5a4 4 0 0 0-5.7-5.6l-2 2M14 11a4 4 0 0 0-5.6 0l-3.5 3.5a4 4 0 0 0 5.7 5.6l2-2"/>',
  layers: '<path d="m3 9 9-5 9 5-9 5-9-5Zm0 5 9 5 9-5"/>',
  play: '<path d="m8 5 11 7-11 7z"/>',
  pause: '<path d="M8 5v14M16 5v14"/>',
};
function icon(which) {
  return (
    '<svg viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.55" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">' +
    icons[which] +
    "</svg>"
  );
}

const format = new Intl.NumberFormat("en-US");
const rulerFormat = new Intl.NumberFormat("en-US", { maximumSignificantDigits: 3 });

export function orbitRadius(planets, index, proportional) {
  if (proportional) return (OUTER_ORBIT_RADIUS * planets[index].au) / planets.at(-1).au;
  return 0.18 + ((index + 1) / planets.length) * 0.8;
}

// The ruler shows a round 1, 2 or 5 in the same projected plane as the map.
export function scaleBar(unitsPerPixel) {
  const maximum = unitsPerPixel * 72;
  const decade = 10 ** Math.floor(Math.log10(maximum));
  const mantissa = maximum / decade;
  const units = (mantissa >= 5 ? 5 : mantissa >= 2 ? 2 : 1) * decade;
  return { units, pixels: units / unitsPerPixel };
}

/** Mount the modal simulator and return an idempotent teardown function. */
export function mountUniverse(container, { onClose, onResearch, locationSearch = "" } = {}) {
  const previousFocus = document.activeElement;
  const modal = document.createElement("dialog");
  modal.className = "os-universe";
  modal.setAttribute("aria-label", "Universe simulator");
  modal.innerHTML = `
    <div class="osu-shell">
      <header class="osu-bar">
        <div class="osu-brand"><span class="osu-logo">${icon("explore")}</span><span>Universe simulator</span></div>
        <nav class="osu-levels" aria-label="Simulation scale">
          <button type="button" data-level="sector">Sector</button><span>/</span>
          <button type="button" data-level="galaxy">Galaxy</button><span>/</span>
          <button type="button" data-level="system">System</button><span>/</span>
          <button type="button" data-level="world">World</button>
        </nav>
        <div class="osu-baractions">
          ${onResearch ? '<button class="osu-atlas-toggle osu-research-switch" data-act="researchFlights" type="button">Research flights</button>' : ""}
          <button class="osu-atlas-toggle" data-act="atlas" type="button" aria-expanded="false" aria-controls="osu-destination-atlas">${icon("layers")}<span>Atlas</span></button>
          <button class="osu-icon" data-act="help" type="button" aria-expanded="false" aria-label="Navigation help" title="Navigation help">${icon("help")}</button>
          <button class="osu-icon" data-act="exit" type="button" aria-label="Close universe simulator" title="Back to observatory">${icon("close")}</button>
        </div>
      </header>
      <div class="osu-stage">
        <canvas tabindex="0" aria-label="Interactive simulated galaxy. Arrow keys pan, plus and minus zoom. Use navigation buttons below to select and explore destinations."></canvas>
        <div class="osu-heading"><span class="osu-eyebrow" data-ref="eyebrow"></span><h2 data-ref="name"></h2><p data-ref="description"></p><p class="osu-location-error" data-ref="locationError" role="status" hidden>This location is unavailable. Showing the opening galaxy.</p>
          <div class="osu-spacing" role="group" aria-label="Orbit spacing" hidden>
            <button type="button" data-spacing="schematic" aria-pressed="true">Even spacing</button>
            <button type="button" data-spacing="au" aria-pressed="false">AU scale</button>
          </div>
        </div>
        <div class="osu-readout" aria-hidden="true"><span class="osu-signal"></span> PROCEDURAL SKY</div>
        <div class="osu-camera" aria-label="Camera controls">
          <button class="osu-icon" data-act="motion" type="button" aria-label="Pause orbits" title="Pause orbits" hidden></button>
          <button class="osu-step" data-act="step" type="button" aria-label="Advance one quarter of the selected orbit" title="Advance one quarter of the selected orbit" hidden>¼ orbit</button>
          <button class="osu-icon" data-act="zoomIn" type="button" aria-label="Zoom in" title="Zoom in (+)">${icon("zoomIn")}</button>
          <button class="osu-icon" data-act="zoomOut" type="button" aria-label="Zoom out" title="Zoom out (−)">${icon("zoomOut")}</button>
          <button class="osu-icon" data-act="reset" type="button" aria-label="Reset camera" title="Reset camera (Home)">${icon("center")}</button>
        </div>
        <aside class="osu-atlas" id="osu-destination-atlas" hidden aria-label="Destination atlas">
          <div class="osu-atlas-head"><strong>Destination atlas</strong><button class="osu-icon" data-act="closeAtlas" type="button" aria-label="Close destination atlas">${icon("close")}</button></div>
          <div class="osu-atlas-places">${atlasDestinations
            .map((d, i) => {
              const world = d.level === "world";
              const title = world ? "Planet " + ROMAN[d.planet] : galaxyName(d.g);
              const type = world ? "Ringed gas giant" : galaxyProfile(d.g).label;
              return `<button type="button" class="osu-place" data-atlas-index="${i}"><canvas aria-hidden="true"></canvas><span><strong>${title}</strong><small>${type} · ${address(d.g, world ? d.s : undefined)}</small></span>${icon("forward")}</button>`;
            })
            .join("")}</div>
          <div class="osu-share">
            <button type="button" class="osu-copy" data-act="copyLocation">${icon("link")}Copy location link</button>
            <span data-ref="shareStatus" role="status"></span>
            <input data-ref="shareUrl" aria-label="Location link to copy" readonly hidden>
          </div>
        </aside>
        <aside class="osu-help" hidden aria-label="Navigation help">
          <strong>Navigate the simulation</strong>
          <p>Drag to pan. Scroll or pinch to zoom. Select a point to inspect it.</p>
          <p>Focus the map to pan with the arrow keys. Use [ and ] to select a destination, Enter to explore, and Home to center.</p>
          <p>Each address generates a fictional system. Its planets follow circular orbits around one star. AU scale preserves their relative orbital radii. Body sizes are enlarged for visibility.</p>
          <details class="osu-model"><summary>How the model works</summary>
            <p>Orbital periods follow P² = a³/M, with P in years, orbital radius a in astronomical units and stellar mass M in solar masses. The simulation uses 365.25 days per year and treats planets as test particles.</p>
            <p>The galaxy ruler measures projected distances in the plane of the screen. Star colors, climates and cloud patterns are illustrative.</p>
          </details>
          <form class="osu-jump"><label>Go to a system address<input name="address" placeholder="G001–S00001" maxlength="14" required autocomplete="off" spellcheck="false"></label><button type="submit">Go</button></form>
          <small data-ref="jumpError" role="status"></small>
        </aside>
        <div class="osu-scale"><i aria-hidden="true"></i><span data-ref="scale"></span></div>
        <div class="osu-nocanvas" hidden>Use the destination controls to explore this simulated sky.</div>
      </div>
      <div class="osu-dock">
        <div class="osu-object"><span class="osu-eyebrow" data-ref="objectKind"></span><h3 data-ref="objectName"></h3><p data-ref="objectDetail"></p></div>
        <div class="osu-worlds" role="group" aria-label="Worlds in orbital order" hidden></div>
        <div class="osu-dockactions">
          <div class="osu-cycle" role="group" aria-label="Choose a destination">
            <button class="osu-icon" data-act="prev" type="button" aria-label="Previous destination" title="Previous destination">${icon("back")}</button>
            <button class="osu-icon" data-act="next" type="button" aria-label="Next destination" title="Next destination">${icon("forward")}</button>
          </div>
          <button class="osu-primary" data-act="enter" type="button"></button>
        </div>
      </div>
      <div class="osu-note"><span data-ref="modelNote">Generated sky</span><div class="osu-time"><span data-ref="elapsed" hidden></span><span data-ref="hint">Drag to pan · Scroll to zoom</span></div></div>
      <span class="osu-sr" data-ref="announce" aria-live="polite"></span>
    </div>`;
  container.append(modal);
  if (modal.showModal) modal.showModal();
  else modal.setAttribute("open", "");
  const $ = (selector) => modal.querySelector(selector);
  const ref = (name) => $('[data-ref="' + name + '"]');
  const act = (name) => $('[data-act="' + name + '"]');
  const canvas = $(".osu-stage > canvas");
  const ctx = canvas.getContext("2d", { alpha: false });
  const help = $(".osu-help");
  const atlas = $(".osu-atlas");
  let paintedAtlas = false;
  const motion = window.matchMedia("(prefers-reduced-motion: reduce)");
  const touch = window.matchMedia("(pointer: coarse)");
  const abort = new AbortController();
  const listen = (node, event, fn, options = {}) =>
    node.addEventListener(event, fn, { ...options, signal: abort.signal });
  const initial = parseUniverseLocation(locationSearch);
  let { level, g: galaxyId, s: systemId, p: planetId } = initial ?? openingLocation;
  ref("locationError").hidden = !!initial;
  let galaxy = level === "galaxy" ? makeGalaxy(galaxyId, !!ctx) : { config: galaxyProfile(galaxyId) },
    system = makeSystem(galaxyId, systemId);
  let backdrop = null,
    backdropId = -1;
  let photosphere = null,
    photosphereSeed = -1;
  let globe = null,
    globeKey = "";
  let pixelRatio = 1;
  let width = 1,
    height = 1,
    frame = 0,
    lastTime = 0,
    days = 0,
    dirty = true;
  let playing = !motion.matches,
    disposed = false;
  let lineupAddress = "";
  let proportional = false;
  let zoom = 1,
    panX = 0,
    panY = 0;
  const pointers = new Map();
  let drag = null,
    pinch = null;
  let hits = [];
  const stage = $(".osu-stage");

  function cleanup() {
    if (disposed) return;
    disposed = true;
    abort.abort();
    observer?.disconnect();
    cancelAnimationFrame(frame);
    if (modal.open && modal.close) modal.close();
    modal.remove();
    if (previousFocus?.isConnected) previousFocus.focus({ preventScroll: true });
  }
  function close() {
    cleanup();
    onClose?.();
  }
  function redraw() {
    dirty = true;
    if (!frame && !disposed) frame = requestAnimationFrame(render);
  }
  function reset() {
    zoom = 1;
    panX = 0;
    panY = 0;
    redraw();
  }
  function baseScale() {
    if (level === "system") return Math.min(width * 0.44, height * 0.55);
    return Math.min(width * 0.46, height * (level === "galaxy" ? 0.57 : 0.49));
  }
  function center() {
    return [width / 2 + panX, height * 0.55 + panY];
  }
  function project(x, y) {
    const [cx, cy] = center(),
      k = baseScale() * zoom;
    return [cx + x * k, cy + y * k];
  }
  function zoomAt(factor, x = width / 2, y = height * 0.55) {
    const next = limit(zoom * factor, 0.55, level === "galaxy" ? 28 : 6);
    const [cx, cy] = center();
    panX += (x - cx) * (1 - next / zoom);
    panY += (y - cy) * (1 - next / zoom);
    zoom = next;
    redraw();
  }
  function update(announce = false) {
    const isSector = level === "sector",
      isSystem = level === "system",
      isWorld = level === "world";
    const body = system.planets[planetId];
    const title = isSector
      ? "Simulated sector"
      : isWorld
        ? "Planet " + ROMAN[planetId]
        : isSystem
          ? starName(systemId)
          : galaxyName(galaxyId);
    ref("name").textContent = title;
    ref("eyebrow").textContent = isSector
      ? "Sector overview"
      : isSystem || isWorld
        ? address(galaxyId, systemId)
        : address(galaxyId) + " / " + galaxy.config.label;
    ref("description").textContent = isSector
      ? format.format(GALAXIES) +
        " galaxies · " +
        format.format(GALAXIES * SYSTEMS) +
        " simulated star systems"
      : isWorld
        ? body.type.name + " · " + starName(systemId)
        : isSystem
          ? system.star.name +
            " · " +
            system.star.mass +
            " solar masses · " +
            system.planets.length +
            " orbiting worlds"
          : format.format(SYSTEMS) +
            " systems · " +
            format.format(galaxy.config.diameter) +
            " ly reference diameter";
    ref("objectKind").textContent = isSector
      ? "Selected galaxy"
      : isWorld
        ? body.type === PLANET_TYPES[2]
          ? "Illustrated cloud tops"
          : "Illustrated surface"
        : isSystem
          ? "Selected world"
          : "Selected star system";
    ref("objectName").textContent = isSector
      ? galaxyName(galaxyId)
      : isWorld
        ? body.type.name
        : isSystem
          ? "Planet " + ROMAN[planetId]
          : starName(systemId);
    ref("objectDetail").textContent = isSector
      ? address(galaxyId) + " · " + galaxy.config.label
      : isWorld
        ? body.au.toFixed(2) +
          " AU from its host star · " +
          format.format(Math.round(body.period)) +
          "-day orbit"
        : isSystem
          ? body.type.name +
            " · " +
            body.au.toFixed(2) +
            " AU · " +
            format.format(Math.round(body.period)) +
            "-day orbit"
          : system.star.name + " · " + system.planets.length + " planets · " + address(galaxyId, systemId);
    act("enter").innerHTML = isWorld
      ? icon("back") + "Return to system"
      : (isSector ? "Enter galaxy" : isSystem ? "Inspect world" : "Enter system") + icon("forward");
    act("motion").hidden = !isSystem;
    act("step").hidden = !isSystem;
    ref("elapsed").hidden = !isSystem || !ctx;
    act("step").title =
      "Advance by " + format.format(Math.round(body.period * 0.25 * 10) / 10) + " model days";
    act("motion").innerHTML = icon(playing ? "pause" : "play");
    act("motion").setAttribute("aria-label", playing ? "Pause orbital motion" : "Play orbital motion");
    act("motion").title = playing ? "Pause orbits" : "Play orbits";
    act("motion").setAttribute("aria-pressed", String(playing));
    $(".osu-spacing").hidden = !isSystem || !ctx;
    for (const button of modal.querySelectorAll("[data-spacing]")) {
      button.setAttribute("aria-pressed", String((button.dataset.spacing === "au") === proportional));
    }
    const lineup = $(".osu-worlds");
    lineup.hidden = !(isSystem || isWorld);
    if (!lineup.hidden) {
      const key = address(galaxyId, systemId);
      if (key !== lineupAddress) {
        lineupAddress = key;
        lineup.replaceChildren(
          ...system.planets.map((planet, i) => {
            const button = document.createElement("button");
            button.type = "button";
            button.dataset.world = i;
            button.dataset.type = PLANET_TYPES.indexOf(planet.type);
            if (planet.ring) button.dataset.ring = "";
            button.setAttribute("aria-label", "Select planet " + ROMAN[i] + ", " + planet.type.name);
            button.title = "Planet " + ROMAN[i] + " · " + planet.type.name;
            const swatch = document.createElement("i");
            swatch.setAttribute("aria-hidden", "true");
            button.append(swatch, ROMAN[i]);
            return button;
          }),
        );
      }
      for (const button of lineup.children)
        button.setAttribute("aria-pressed", String(Number(button.dataset.world) === planetId));
    }
    ref("modelNote").textContent = isSystem
      ? "Generated orbits · body sizes enlarged"
      : isWorld
        ? "Generated globe · illustrative lighting"
        : isSector
          ? "Galaxy markers enlarged"
          : "Generated sky · projected distances";
    ref("hint").textContent = isSystem
      ? playing
        ? DAYS_PER_SECOND + " model days / second"
        : "Orbital motion paused"
      : touch.matches
        ? "Drag to pan · Pinch to zoom"
        : "Drag to pan · Scroll to zoom";
    if (!ctx) ref("hint").textContent = "Use the destination controls";
    for (const button of modal.querySelectorAll("[data-level]")) {
      if (button.dataset.level === level) button.setAttribute("aria-current", "step");
      else button.removeAttribute("aria-current");
    }
    canvas.setAttribute(
      "aria-label",
      title +
        ", interactive simulated map. Arrow keys pan, plus and minus zoom. Destination controls are below the map.",
    );
    if (announce)
      ref("announce").textContent =
        title + ". " + ref("objectName").textContent + ". " + ref("objectDetail").textContent;
    redraw();
  }
  function currentHash() {
    return universeLocationHash({ level, g: galaxyId, s: systemId, p: planetId });
  }
  function reflectLocation() {
    ref("locationError").hidden = true;
    // Navigating inside a visit keeps Back pointed at the preceding page.
    if (location.hash !== currentHash()) history.replaceState(history.state, "", currentHash());
    ref("shareStatus").textContent = "";
    ref("shareUrl").hidden = true;
  }
  async function copyLocation() {
    const url = new URL(location.href);
    url.hash = currentHash();
    const field = ref("shareUrl");
    try {
      await navigator.clipboard.writeText(url.href);
      if (disposed) return;
      field.hidden = true;
      ref("shareStatus").textContent = "Location link copied.";
    } catch {
      if (disposed) return;
      field.value = url.href;
      field.hidden = false;
      field.focus();
      field.select();
      ref("shareStatus").textContent = "Copy the selected link.";
    }
  }
  function select(id) {
    if (level === "sector") {
      if (id !== galaxyId) {
        galaxyId = id;
        // Delay making a galaxy image until the visitor enters it.
        galaxy = { config: galaxyProfile(id) };
        system = makeSystem(galaxyId, systemId);
        planetId = 0;
        days = 0;
      }
    } else if (level === "galaxy") {
      systemId = id;
      system = makeSystem(galaxyId, systemId);
      planetId = 0;
      days = 0;
    } else planetId = id;
    reflectLocation();
    update(true);
  }
  function cycle(delta) {
    const count = level === "sector" ? GALAXIES : level === "galaxy" ? SYSTEMS : system.planets.length;
    const current = level === "sector" ? galaxyId : level === "galaxy" ? systemId : planetId;
    select((current + delta + count) % count);
    // Keep a keyboard-selected point in view, including at high magnification.
    if (level === "sector" || level === "galaxy") {
      const x = level === "sector" ? sector[galaxyId].x : galaxy.position[systemId * 2];
      const y = level === "sector" ? sector[galaxyId].y : galaxy.position[systemId * 2 + 1];
      const [sx, sy] = project(x, y);
      if (sx < 36 || sx > width - 36 || sy < 42 || sy > height - 42) {
        panX -= sx - width * 0.5;
        panY -= sy - height * 0.55;
      }
    }
  }
  function go(next) {
    if (!atlas.hidden) toggleAtlas();
    if (!help.hidden) toggleHelp();
    if (next === "galaxy" && !galaxy.position) galaxy = makeGalaxy(galaxyId, !!ctx);
    level = next;
    reflectLocation();
    reset();
    lastTime = 0;
    update(true);
  }
  function enter() {
    go(level === "galaxy" || level === "world" ? "system" : level === "system" ? "world" : "galaxy");
  }
  function toggleAtlas() {
    atlas.hidden = !atlas.hidden;
    act("atlas").setAttribute("aria-expanded", String(!atlas.hidden));
    if (atlas.hidden) {
      act("atlas").focus();
      return;
    }
    help.hidden = true;
    act("help").setAttribute("aria-expanded", "false");
    if (!paintedAtlas && ctx) {
      for (const [i, thumb] of Array.from(atlas.querySelectorAll("canvas")).entries()) {
        const d = atlasDestinations[i];
        thumb.width = 156;
        thumb.height = 88;
        const brush = thumb.getContext("2d");
        if (!brush) continue;
        if (d.level === "world") {
          const planet = makeSystem(d.g, d.s).planets[d.planet];
          drawWorldPortrait(brush, 78, 44, 34, planet, makeWorld(d.g, d.s, d.planet, planet, 112));
        } else {
          const config = galaxyProfile(d.g);
          const sprite = makeDisk(d.g, config, 152);
          brush.drawImage(sprite, -6, -39, 168, 168);
        }
      }
      paintedAtlas = true;
    }
    atlas.querySelector(".osu-place").focus();
  }
  function visitDestination(index) {
    const d = atlasDestinations[index];
    if (!d) return;
    if (d.g !== galaxyId) galaxy = { config: galaxyProfile(d.g) };
    galaxyId = d.g;
    systemId = d.s;
    planetId = d.planet ?? 0;
    days = 0;
    system = makeSystem(galaxyId, systemId);
    toggleAtlas();
    go(d.level);
    canvas.focus();
  }
  function toggleHelp() {
    atlas.hidden = true;
    act("atlas").setAttribute("aria-expanded", "false");
    help.hidden = !help.hidden;
    act("help").setAttribute("aria-expanded", String(!help.hidden));
    if (!help.hidden) $(".osu-jump input").focus();
    else act("help").focus();
  }
  listen($(".osu-jump"), "submit", (event) => {
    event.preventDefault();
    const parsed = parseSystemAddress($(".osu-jump input").value);
    if (!parsed) {
      ref("jumpError").textContent = "Use G001–G128 and S00001–S32768, for example G013–S02189.";
      return;
    }
    galaxyId = parsed.g;
    systemId = parsed.s;
    planetId = 0;
    days = 0;
    galaxy = { config: galaxyProfile(galaxyId) };
    system = makeSystem(galaxyId, systemId);
    ref("jumpError").textContent = "";
    toggleHelp();
    go("system");
    canvas.focus();
  });
  listen(modal, "click", (event) => {
    const button = event.target.closest("button");
    if (!button) return;
    if (button.dataset.world !== undefined) {
      select(Number(button.dataset.world));
      return;
    }
    if (button.dataset.atlasIndex !== undefined) {
      visitDestination(Number(button.dataset.atlasIndex));
      return;
    }
    if (button.dataset.spacing) {
      proportional = button.dataset.spacing === "au";
      reset();
      update(true);
      return;
    }
    if (button.dataset.level) {
      go(button.dataset.level);
      return;
    }
    switch (button.dataset.act) {
      case "researchFlights":
        onResearch?.(currentHash());
        break;
      case "exit":
        close();
        break;
      case "copyLocation":
        void copyLocation();
        break;
      case "atlas":
      case "closeAtlas":
        toggleAtlas();
        break;
      case "help":
        toggleHelp();
        break;
      case "prev":
        cycle(-1);
        break;
      case "next":
        cycle(1);
        break;
      case "enter":
        enter();
        break;
      case "motion":
        playing = !playing;
        lastTime = 0;
        update(true);
        break;
      case "step":
        days += system.planets[planetId].period / 4;
        lastTime = 0;
        update();
        ref("announce").textContent =
          "Advanced one quarter orbit to model day " + format.format(Math.round(days * 10) / 10) + ".";
        break;
      case "zoomIn":
        zoomAt(1.55);
        break;
      case "zoomOut":
        zoomAt(1 / 1.55);
        break;
      case "reset":
        reset();
        break;
    }
  });
  listen(modal, "cancel", (event) => {
    event.preventDefault();
    if (!atlas.hidden) toggleAtlas();
    else if (!help.hidden) toggleHelp();
    else close();
  });
  listen(canvas, "keydown", (event) => {
    const step = event.shiftKey ? 100 : 40;
    switch (event.key) {
      case "ArrowLeft":
        panX += step;
        break;
      case "ArrowRight":
        panX -= step;
        break;
      case "ArrowUp":
        panY += step;
        break;
      case "ArrowDown":
        panY -= step;
        break;
      case "+":
      case "=":
        zoomAt(1.35);
        break;
      case "-":
        zoomAt(1 / 1.35);
        break;
      case "Home":
        reset();
        break;
      case "[":
        cycle(-1);
        break;
      case "]":
        cycle(1);
        break;
      case "Enter":
        enter();
        break;
      default:
        return;
    }
    event.preventDefault();
    redraw();
  });
  function localPoint(event) {
    const rect = canvas.getBoundingClientRect();
    return { x: event.clientX - rect.left, y: event.clientY - rect.top };
  }
  function twoPointers() {
    const [a, b] = [...pointers.values()];
    return {
      x: (a.x + b.x) / 2,
      y: (a.y + b.y) / 2,
      distance: Math.max(1, Math.hypot(a.x - b.x, a.y - b.y)),
    };
  }
  listen(
    canvas,
    "wheel",
    (event) => {
      event.preventDefault();
      const { x, y } = localPoint(event);
      const units = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? height : 1;
      zoomAt(Math.exp(-limit(event.deltaY * units, -220, 220) * 0.0023), x, y);
    },
    { passive: false },
  );
  listen(canvas, "pointerdown", (event) => {
    if (event.button !== 0 || pointers.size >= 2) return;
    canvas.focus({ preventScroll: true });
    canvas.setPointerCapture(event.pointerId);
    const pos = localPoint(event);
    pointers.set(event.pointerId, pos);
    drag = { ...pos, moved: false };
    if (pointers.size === 2) {
      pinch = twoPointers();
      drag.moved = true;
    }
  });
  listen(canvas, "pointermove", (event) => {
    if (!pointers.has(event.pointerId)) return;
    const previous = pointers.get(event.pointerId);
    const pos = localPoint(event);
    pointers.set(event.pointerId, pos);
    if (pointers.size >= 2) {
      const next = twoPointers();
      if (pinch) {
        panX += next.x - pinch.x;
        panY += next.y - pinch.y;
        zoomAt(next.distance / pinch.distance, next.x, next.y);
      }
      pinch = next;
      drag.moved = true;
    } else {
      if (Math.hypot(pos.x - drag.x, pos.y - drag.y) > 4) drag.moved = true;
      if (drag.moved) {
        panX += pos.x - previous.x;
        panY += pos.y - previous.y;
        redraw();
      }
    }
  });
  function finishPointer(event) {
    if (!pointers.has(event.pointerId)) return;
    if (event.type === "pointerup" && !drag.moved && pointers.size === 1) {
      const pos = localPoint(event);
      pick(pos.x, pos.y);
    }
    pointers.delete(event.pointerId);
    if (pointers.size === 1 && pinch) drag = { ...pointers.values().next().value, moved: true };
    if (pointers.size < 2) pinch = null;
    if (!pointers.size) drag = null;
  }
  listen(canvas, "pointerup", finishPointer);
  listen(canvas, "pointercancel", finishPointer);
  listen(canvas, "lostpointercapture", finishPointer);

  function pick(x, y) {
    if (level === "world") return;
    let nearest = -1,
      best = level === "system" ? 34 : 18;
    if (level === "galaxy") {
      const [cx, cy] = center(),
        k = baseScale() * zoom;
      const mx = (x - cx) / k,
        my = (y - cy) / k;
      let distance = (best / k) ** 2;
      for (let i = 0; i < SYSTEMS; i++) {
        const dx = galaxy.position[i * 2] - mx,
          dy = galaxy.position[i * 2 + 1] - my;
        const d = dx * dx + dy * dy;
        if (d < distance) {
          nearest = i;
          distance = d;
        }
      }
    } else {
      for (const [sx, sy, id] of hits) {
        const d = Math.hypot(sx - x, sy - y);
        if (d < best) {
          best = d;
          nearest = id;
        }
      }
    }
    if (nearest >= 0) select(nearest);
  }
  function marker(x, y, label, radius = 10) {
    if (x < 0 || x > width || y < 0 || y > height) return;
    ctx.save();
    ctx.lineWidth = 1;
    ctx.strokeStyle = "#ffac76";
    ctx.fillStyle = "#fbeaDC";
    ctx.beginPath();
    for (let j = 0; j < 4; j++) {
      const a = ((j + 0.12) * Math.PI) / 2;
      ctx.moveTo(x + Math.cos(a) * radius, y + Math.sin(a) * radius);
      ctx.arc(x, y, radius, a, a + Math.PI * 0.32);
    }
    ctx.stroke();
    ctx.font = "10px ui-monospace, monospace";
    const measured = ctx.measureText(label).width;
    const left = x + radius + 9 + measured > width - 6;
    const tx = left ? x - radius - measured - 8 : x + radius + 8;
    ctx.fillStyle = "rgba(5,8,13,.83)";
    ctx.fillRect(tx - 4, y - 10, measured + 8, 19);
    ctx.fillStyle = "#ffe5d2";
    ctx.fillText(label, tx, y + 3.2);
    ctx.restore();
  }
  function drawSky() {
    if (backdropId !== galaxyId) {
      backdrop = makeBackdrop(galaxyId);
      backdropId = galaxyId;
    }
    const bx = -28 + limit(panX * 0.025, -26, 26);
    const by = -28 + limit(panY * 0.025, -26, 26);
    ctx.drawImage(backdrop, bx, by, width + 56, height + 56);
    for (const [u, v, b, z] of sky) {
      const x = (((u * width + panX * z * 0.025) % width) + width) % width;
      const y = (((v * height + panY * z * 0.025) % height) + height) % height;
      const size = b > 0.982 ? 1.8 : 0.45 + b * 0.8;
      ctx.globalAlpha = 0.17 + b * 0.57;
      ctx.fillStyle = z > 0.67 ? "#f6ddc8" : "#c1d5ea";
      ctx.fillRect(x, y, size, size);
      if (b > 0.982) {
        ctx.globalAlpha = 0.22;
        ctx.fillRect(x - 2.3, y + 0.5, 6, 0.55);
        ctx.fillRect(x + 0.5, y - 2.3, 0.55, 6);
      }
    }
    ctx.globalAlpha = 1;
  }
  function showRuler(unitsPerPixel, unit) {
    const bar = scaleBar(unitsPerPixel);
    const line = $(".osu-scale i");
    line.hidden = false;
    line.style.width = bar.pixels + "px";
    ref("scale").textContent = rulerFormat.format(bar.units) + " " + unit;
  }
  function showScaleNote(note) {
    $(".osu-scale i").hidden = true;
    ref("scale").textContent = note;
  }
  function drawGalaxy() {
    const [cx, cy] = center(),
      k = baseScale() * zoom;
    const extent = k * MAP_EXTENT;
    ctx.globalCompositeOperation = "screen";
    ctx.drawImage(
      zoom <= 2.6 ? galaxy.texture : galaxy.haze,
      cx - extent,
      cy - extent,
      extent * 2,
      extent * 2,
    );
    ctx.globalCompositeOperation = "source-over";
    if (zoom > 1.65) {
      ctx.globalAlpha = Math.min(0.87, (zoom - 1.65) * 0.35);
      for (let i = 0; i < SYSTEMS; i++) {
        const x = cx + galaxy.position[i * 2] * k,
          y = cy + galaxy.position[i * 2 + 1] * k;
        if (x < -2 || x > width || y < -2 || y > height) continue;
        ctx.fillStyle = galaxy.inks[galaxy.tint[i]];
        const w = 0.45 + galaxy.magnitude[i] * Math.min(1.4, zoom * 0.38);
        ctx.fillRect(x, y, w, w);
      }
      ctx.globalAlpha = 1;
    }
    marker(
      cx + galaxy.position[systemId * 2] * k,
      cy + galaxy.position[systemId * 2 + 1] * k,
      "S" + pad(systemId, 5),
    );
    showRuler(galaxy.config.diameter / (2 * k), "ly");
  }
  function drawSector() {
    hits = [];
    const tints = {
      elliptical: "250,194,137",
      ring: "118,209,235",
      edge: "211,191,177",
      barred: "221,177,146",
      spiral: "182,193,239",
    };
    sector.forEach((g, id) => {
      const [x, y] = project(g.x, g.y);
      if (x < -20 || x > width + 20 || y < -20 || y > height + 20) return;
      const config = galaxyProfile(id);
      const radius = g.size * ((width < 600 ? 1.15 : 1.68) + Math.sqrt(zoom) * 0.8);
      hits.push([x, y, id]);
      ctx.save();
      ctx.translate(x, y);
      ctx.rotate(config.angle);
      ctx.scale(1, config.tilt);
      glow(ctx, 0, 0, radius, tints[config.family], 0.64);
      if (config.family === "ring") {
        ctx.strokeStyle = "rgba(153,215,245,.73)";
        ctx.lineWidth = 0.95;
        ctx.beginPath();
        ctx.arc(0, 0, radius * 0.55, 0, TAU);
        ctx.stroke();
      } else if (config.arms) {
        ctx.strokeStyle = "rgba(" + tints[config.family] + ",.70)";
        ctx.lineWidth = Math.max(0.65, radius * 0.065);
        for (let arm = 0; arm < config.arms; arm++) {
          ctx.beginPath();
          for (let j = 0; j < 13; j++) {
            const r = 0.06 + j * 0.062;
            const a = armAngle(r, arm, config);
            const px = radius * r * Math.cos(a),
              py = radius * r * Math.sin(a);
            if (!j) ctx.moveTo(px, py);
            else ctx.lineTo(px, py);
          }
          ctx.stroke();
        }
      }
      glow(ctx, 0, 0, radius * 0.35, "255,224,181", 0.85);
      ctx.fillStyle = "#fff1d4";
      ctx.beginPath();
      ctx.arc(0, 0, Math.max(0.75, radius * 0.085), 0, TAU);
      ctx.fill();
      ctx.restore();
    });
    const g = sector[galaxyId];
    marker(...project(g.x, g.y), address(galaxyId), 12);
    showRuler(SECTOR_MLY_PER_UNIT / (baseScale() * zoom), "Mly");
  }
  function stellarSurface() {
    const seed = systemSeed(galaxyId, systemId);
    if (photosphereSeed === seed) return photosphere;
    photosphereSeed = seed;
    const tile = canvasSurface(128),
      brush = tile.getContext("2d");
    const pixels = brush.createImageData(128, 128),
      rng = random(seed);
    const color = Number.parseInt(system.star.color.slice(1), 16);
    const rgb = [color >> 16, (color >> 8) & 255, color & 255];
    for (let y = 0; y < 128; y++) {
      for (let x = 0; x < 128; x++) {
        const r = Math.hypot((x - 63.5) / 63, (y - 63.5) / 63);
        if (r > 1) continue;
        const z = Math.sqrt(1 - r * r);
        const granulation = (rng() - 0.5) * (0.1 + 0.12 * z);
        const intensity = 0.63 + 0.45 * z + granulation;
        const o = 4 * (y * 128 + x);
        for (let c = 0; c < 3; c++) pixels.data[o + c] = rgb[c] * intensity + 14 * z;
        pixels.data[o + 3] = 255 * Math.min(1, (1 - r) * 64);
      }
    }
    brush.putImageData(pixels, 0, 0);
    photosphere = tile;
    return tile;
  }
  function drawSystem() {
    const clockLabel = "Model day " + format.format(Math.floor(days));
    if (ref("elapsed").textContent !== clockLabel) ref("elapsed").textContent = clockLabel;
    const [cx, cy] = center(),
      k = baseScale() * zoom;
    const tilt = 0.51,
      turn = -0.22;
    const radii = system.planets.map((_, i) => orbitRadius(system.planets, i, proportional));
    const bodies = system.planets.map((planet, i) => {
      const phase = planet.phase + (days / planet.period) * TAU;
      const dx = radii[i] * Math.cos(phase),
        dy = radii[i] * Math.sin(phase) * tilt;
      return {
        id: i,
        phase,
        x: cx + k * (dx * Math.cos(turn) - dy * Math.sin(turn)),
        y: cy + k * (dx * Math.sin(turn) + dy * Math.cos(turn)),
        depth: dy,
        planet,
      };
    });
    ctx.save();
    ctx.translate(cx, cy);
    ctx.rotate(turn);
    radii.forEach((r, i) => {
      ctx.beginPath();
      ctx.ellipse(0, 0, k * r, k * r * tilt, 0, 0, TAU);
      ctx.strokeStyle = i === planetId ? "rgba(239,159,104,.58)" : "rgba(189,203,224,.22)";
      ctx.lineWidth = i === planetId ? 1.1 : 0.75;
      ctx.stroke();
    });
    // Mark the selected planet's approach along the same projected orbit.
    const selectedPhase = bodies[planetId].phase;
    const orbitPixels = k * radii[planetId];
    ctx.lineWidth = 2;
    for (let j = 0; j < 12; j++) {
      const start = selectedPhase - (12 - j) * 0.045;
      ctx.strokeStyle = "rgba(255,174,115," + (j + 1) / 18 + ")";
      ctx.beginPath();
      ctx.ellipse(0, 0, orbitPixels, orbitPixels * tilt, 0, start, start + 0.046);
      ctx.stroke();
    }
    ctx.restore();
    const starRadius = Math.max(5.5, Math.min(35, k * 0.039));
    const starRGB = parseInt(system.star.color.slice(1), 16);
    const rgb = [(starRGB >> 16) & 255, (starRGB >> 8) & 255, starRGB & 255].join(",");
    function drawStar() {
      glow(ctx, cx, cy, starRadius * 10.4, rgb, 0.46);
      glow(ctx, cx, cy, starRadius * 4.1, rgb, 0.88);
      const corona = ctx.createRadialGradient(cx, cy, starRadius * 0.98, cx, cy, starRadius * 1.78);
      corona.addColorStop(0, "rgba(" + rgb + ",.89)");
      corona.addColorStop(0.13, "rgba(" + rgb + ",.15)");
      corona.addColorStop(1, "rgba(" + rgb + ",0)");
      ctx.fillStyle = corona;
      ctx.fillRect(cx - starRadius * 2, cy - starRadius * 2, starRadius * 4, starRadius * 4);
      ctx.drawImage(stellarSurface(), cx - starRadius, cy - starRadius, starRadius * 2, starRadius * 2);
    }
    function drawBody({ x, y, planet: p }) {
      const r = Math.max(3.7, Math.min(45, p.size * k));
      const direction = Math.atan2(cy - y, cx - x);
      const sx = x + Math.cos(direction) * r * 0.76,
        sy = y + Math.sin(direction) * r * 0.76;
      const grad = ctx.createRadialGradient(sx, sy, 0, sx, sy, r * 1.88);
      grad.addColorStop(0, p.type.colors[1]);
      grad.addColorStop(0.48, p.type.colors[0]);
      grad.addColorStop(0.69, p.type.colors[2]);
      grad.addColorStop(1, "#080d17");
      function drawRing(front) {
        ctx.strokeStyle = "rgba(230,203,167,.74)";
        ctx.lineWidth = Math.max(1, r * 0.19);
        ctx.beginPath();
        ctx.ellipse(x, y, r * 1.75, r * 0.38, -0.38, front ? 0 : Math.PI, front ? Math.PI : TAU);
        ctx.stroke();
      }
      if (p.ring) drawRing(false);
      ctx.fillStyle = grad;
      ctx.beginPath();
      ctx.arc(x, y, r, 0, TAU);
      ctx.fill();
      if (p.ring) drawRing(true);
    }
    bodies.sort((a, b) => a.depth - b.depth);
    bodies.filter((b) => b.depth < 0).forEach(drawBody);
    drawStar();
    bodies.filter((b) => b.depth >= 0).forEach(drawBody);
    hits = bodies.map((b) => [b.x, b.y, b.id]);
    const selected = bodies.find((b) => b.id === planetId);
    marker(
      selected.x,
      selected.y,
      "PLANET " + ROMAN[planetId] + (proportional ? " · " + selected.planet.au.toFixed(2) + " AU" : ""),
      Math.max(12, Math.min(45, selected.planet.size * k) + 5),
    );
    if (proportional) showRuler(system.planets.at(-1).au / (OUTER_ORBIT_RADIUS * k), "AU");
    else showScaleNote("SCHEMATIC");
  }
  function drawWorld() {
    const planet = system.planets[planetId],
      key = [galaxyId, systemId, planetId].join(".");
    if (globeKey !== key) {
      globe = makeWorld(galaxyId, systemId, planetId, planet);
      globeKey = key;
    }
    const [cx, cy] = center(),
      radius = Math.min(width * (planet.ring ? 0.27 : 0.36), height * 0.38) * zoom;
    drawWorldPortrait(ctx, cx, cy, radius, planet, globe);
    showScaleNote(planet.type === PLANET_TYPES[2] ? "CLOUD TOPS" : "SURFACE VIEW");
  }
  function render(now) {
    frame = 0;
    if (disposed || !ctx || document.hidden) return;
    const active = level === "system" && playing;
    if (dirty || (active && now - lastTime >= 32)) {
      if (active && lastTime) days += ((now - lastTime) / 1000) * DAYS_PER_SECOND;
      lastTime = now;
      dirty = false;
      // Resizing clears the backing store. Resize and paint in the same frame,
      // especially when the planet chooser changes the stage height on phones.
      const nextWidth = Math.max(1, Math.round(width * pixelRatio));
      const nextHeight = Math.max(1, Math.round(height * pixelRatio));
      if (canvas.width !== nextWidth || canvas.height !== nextHeight) {
        canvas.width = nextWidth;
        canvas.height = nextHeight;
      }
      ctx.setTransform(pixelRatio, 0, 0, pixelRatio, 0, 0);
      drawSky();
      if (level === "sector") drawSector();
      else if (level === "galaxy") drawGalaxy();
      else if (level === "world") drawWorld();
      else drawSystem();
    }
    if (active) frame = requestAnimationFrame(render);
  }
  function sizeCanvas() {
    const rect = stage.getBoundingClientRect();
    pixelRatio = Math.min(window.devicePixelRatio || 1, 2);
    width = Math.max(1, rect.width);
    height = Math.max(1, rect.height);
    redraw();
  }
  const observer = typeof ResizeObserver === "function" ? new ResizeObserver(sizeCanvas) : null;
  observer?.observe(stage);
  if (!observer) listen(window, "resize", sizeCanvas);
  if (!ctx) {
    $(".osu-nocanvas").hidden = false;
    $(".osu-camera").hidden = true;
    $(".osu-scale").hidden = true;
  }
  listen(document, "visibilitychange", () => {
    lastTime = 0;
    redraw();
  });
  listen(motion, "change", () => {
    playing = !motion.matches;
    lastTime = 0;
    update();
  });
  sizeCanvas();
  update();
  act("exit").focus({ preventScroll: true });
  return cleanup;
}
