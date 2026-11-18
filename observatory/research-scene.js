import { icon } from "./icons.js";

const FLIGHT_MS = 1800;
const TAU = 2 * Math.PI;
const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));
const mix = (a, b, t) => a + (b - a) * t;
const smooth = (t) => t * t * (3 - 2 * t);
const palettes = {
  mars: [196, 116, 74],
  sun: [210, 104, 71],
  quasar: [149, 134, 204],
  dwarf: [104, 152, 211],
  moon: [154, 168, 188],
  star: [199, 157, 120],
  galaxy: [145, 133, 198],
  comet: [100, 173, 172],
  transient: [167, 141, 194],
};

function randomSource(value) {
  let seed = 2166136261;
  for (const char of value) seed = Math.imul(seed ^ char.charCodeAt(0), 16777619);
  return () => {
    seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
    return seed / 0x100000000;
  };
}

function element(tag, cls, text) {
  const node = document.createElement(tag);
  if (cls) node.className = cls;
  if (text) node.textContent = text;
  return node;
}

function mediaFor(destination, plate) {
  if (plate?.file)
    return {
      src: "./assets/" + plate.file,
      href: plate.source,
      credit: plate.credit,
      kind: plate.kind,
      alt: plate.description,
      markers: plate.markers,
      fieldArcsec: plate.fieldArcsec,
      mode:
        destination.id === "mars" || destination.id === "sun"
          ? "disc"
          : destination.id === "dwarf"
            ? "illustration"
            : "field",
    };
  const media = destination.image?.src ? destination.image : destination.figure;
  if (!media?.src) return null;
  return {
    src: media.src,
    href: media.url,
    credit: media.credit || "Repository figure",
    kind: media.type === "archive-image" ? "Archive image" : "Research figure",
    alt: media.alt || media.caption || destination.name,
    mode: media.type === "archive-image" ? "field" : "figure",
  };
}

// A deterministic navigation backdrop. The observational images are separate
// DOM layers, so their pixels never inherit the renderer's color or noise.
function dustPlate(random, tint) {
  const dust = document.createElement("canvas");
  dust.width = 416;
  dust.height = 244;
  const ctx = dust.getContext("2d");
  if (!ctx) return dust;
  const grid = Float32Array.from({ length: 64 * 64 }, random);
  const sample = (x, y) => {
    const ix = Math.floor(x),
      iy = Math.floor(y),
      fx = smooth(x - ix),
      fy = smooth(y - iy);
    const a = (iy & 63) << 6,
      b = ((iy + 1) & 63) << 6;
    const u = ix & 63,
      v = (ix + 1) & 63;
    return mix(mix(grid[a + u], grid[a + v], fx), mix(grid[b + u], grid[b + v], fx), fy);
  };
  const pixels = ctx.createImageData(dust.width, dust.height);
  const slope = (random() > 0.5 ? 1 : -1) * 0.33;
  const bend = random() * 2;
  for (let y = 0; y < dust.height; y++) {
    const v = y / dust.height;
    for (let x = 0; x < dust.width; x++) {
      const u = x / dust.width;
      const coarse = sample(u * 4 + 12, v * 4);
      const mid = sample(u * 11 + coarse * 2, v * 11 + coarse * 2);
      const fine = sample(u * 33 + 1.7 * mid, v * 33);
      const lane = v - 0.58 - slope * (u - 0.5) - 0.095 * Math.sin(6 * u + bend) - 0.12 * coarse;
      const body = Math.exp(-lane * lane * 17);
      const light = body * Math.max(0, coarse * 0.46 + mid * 0.36 + fine * 0.18 - 0.27);
      const shadow = 1 - 0.76 * Math.exp(-Math.pow((lane + 0.05 + 0.069 * mid) * 17, 2));
      const i = (y * dust.width + x) * 4;
      const grain = random() * 1.4;
      pixels.data[i] = 5 + grain + light * shadow * tint[0] * 0.45;
      pixels.data[i + 1] = 8 + grain + light * shadow * tint[1] * 0.44;
      pixels.data[i + 2] = 14 + grain + light * shadow * tint[2] * 0.51;
      pixels.data[i + 3] = 255;
    }
  }
  ctx.putImageData(pixels, 0, 0);
  return dust;
}

function starCatalog(random) {
  return Array.from({ length: 1050 }, () => {
    const blue = random();
    return {
      x: random() * 2.8 - 1.4,
      y: random() * 2.8 - 1.4,
      z: random(),
      radius: 0.25 + Math.pow(random(), 4) * 1.08,
      strength: 0.25 + random() * 0.67,
      color: blue > 0.79 ? "#f7d4aa" : blue < 0.35 ? "#a7c3e1" : "#e7e5df",
    };
  });
}

function drawField(ctx, dust, stars, view) {
  const { w, h, focusX, focusY, panX, panY, zoom, t, speed } = view;
  ctx.globalAlpha = 1;
  ctx.fillStyle = "#070a0f";
  ctx.fillRect(0, 0, w, h);
  const edge = Math.max(w / dust.width, h / dust.height) * (1.17 + 0.05 * (zoom - 1));
  const width = dust.width * edge,
    height = dust.height * edge;
  ctx.drawImage(dust, (w - width) / 2 + panX * 0.065, (h - height) / 2 + panY * 0.055, width, height);
  const perspective = Math.max(w * 0.56, h * 0.8) * (1 + 0.13 * (zoom - 1));
  const travel = smooth(t) * 2.4;
  const n = Math.min(stars.length, Math.max(430, Math.round((w * h) / 510)));
  for (let i = 0; i < n; i++) {
    const star = stars[i];
    const z = 0.14 + ((star.z - (travel % 1) + 1) % 1) * 1.56;
    const scale = perspective / z;
    const px = focusX + star.x * scale + (panX * 0.14) / z;
    const py = focusY + star.y * scale + (panY * 0.14) / z;
    if (px < -7 || py < -7 || px > w + 7 || py > h + 7) continue;
    const alpha = star.strength * Math.min(1, (1.7 - z) * 7) * Math.min(1, 0.65 / z + 0.35);
    ctx.fillStyle = ctx.strokeStyle = star.color;
    ctx.globalAlpha = alpha;
    const radius = Math.max(0.34, star.radius / Math.sqrt(z));
    if (speed > 0.012) {
      const trail = Math.min(0.27, (0.055 * speed) / z);
      ctx.lineWidth = Math.min(1.5, radius);
      ctx.beginPath();
      ctx.moveTo(px - (px - focusX) * trail, py - (py - focusY) * trail);
      ctx.lineTo(px, py);
      ctx.stroke();
    }
    ctx.beginPath();
    ctx.arc(px, py, radius, 0, TAU);
    ctx.fill();
    if (radius > 1.4 && alpha > 0.56) {
      ctx.globalAlpha *= 0.26;
      ctx.fillRect(px - 3 * radius, py - 0.4, 6 * radius, 0.8);
      ctx.fillRect(px - 0.4, py - 2 * radius, 0.8, 4 * radius);
    }
  }
  ctx.globalAlpha = 1;
}

function glow(ctx, x, y, r, color, alpha) {
  const g = ctx.createRadialGradient(x, y, 0, x, y, r);
  g.addColorStop(0, `rgba(${color},${alpha})`);
  g.addColorStop(0.48, `rgba(${color},${alpha * 0.26})`);
  g.addColorStop(1, `rgba(${color},0)`);
  ctx.fillStyle = g;
  ctx.fillRect(x - r, y - r, r * 2, r * 2);
}

// Geometrical collecting aperture for reports with no published image.
// Its petals depict the simulator's instrument, not an astronomical object.
function drawAperture(ctx, x, y, r, tint) {
  if (r < 2) return;
  const color = tint.join(",");
  ctx.save();
  ctx.translate(x, y);
  ctx.rotate(-0.16);
  ctx.scale(1, 0.84);
  glow(ctx, 0, 0, r * 1.55, color, 0.25);
  const ring = ctx.createLinearGradient(-r, -r, r, r);
  ring.addColorStop(0, "#8b9aab");
  ring.addColorStop(0.23, "#242e3d");
  ring.addColorStop(0.46, "#11151c");
  ring.addColorStop(0.76, `rgb(${color})`);
  ring.addColorStop(1, "#252f3f");
  ctx.fillStyle = ring;
  ctx.beginPath();
  ctx.arc(0, 0, r, 0, TAU);
  ctx.fill();
  ctx.fillStyle = "#0d1520";
  ctx.beginPath();
  ctx.arc(0, 0, r * 0.94, 0, TAU);
  ctx.fill();
  for (let i = 0; i < 10; i++) {
    const a = (i * TAU) / 10;
    ctx.save();
    ctx.rotate(a);
    const blade = ctx.createLinearGradient(0, -r * 0.8, r * 0.4, 0);
    blade.addColorStop(0, `rgba(${color},${0.24 + 0.19 * (1 + Math.cos(a))})`);
    blade.addColorStop(1, "#0c1521");
    ctx.fillStyle = blade;
    ctx.strokeStyle = "#070c13";
    ctx.lineWidth = Math.max(0.5, r * 0.006);
    ctx.beginPath();
    ctx.moveTo(0, -r * 0.91);
    ctx.arc(0, 0, r * 0.91, -Math.PI / 2, -Math.PI / 2 + TAU / 10);
    ctx.lineTo(r * 0.18, -r * 0.36);
    ctx.lineTo(-r * 0.24, -r * 0.3);
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
    ctx.restore();
  }
  const glass = ctx.createRadialGradient(-r * 0.08, -r * 0.15, 0, 0, 0, r * 0.44);
  glass.addColorStop(0, "#273f58");
  glass.addColorStop(0.52, "#0b1725");
  glass.addColorStop(0.87, "#080d15");
  glass.addColorStop(1, `rgba(${color},.65)`);
  ctx.fillStyle = glass;
  ctx.beginPath();
  ctx.arc(0, 0, r * 0.37, 0, TAU);
  ctx.fill();
  ctx.strokeStyle = "rgba(207,221,237,.48)";
  ctx.lineWidth = Math.max(0.6, r * 0.006);
  ctx.beginPath();
  ctx.arc(0, 0, r * 0.94, 3.28, 4.8);
  ctx.stroke();
  ctx.restore();
}

// Orthographic projection of the published equatorial coordinates. East
// increases to the left, the convention used when looking out at the sky.
function drawPointing(ctx, x, y, r, coords, approach, yaw, pitch) {
  if (r < 3) return;
  const rad = Math.PI / 180;
  const longitude = (coords.raDeg + 28 + yaw + (1 - approach) * 55) * rad;
  const latitude = clamp(coords.decDeg * 0.7 + pitch, -85, 85) * rad;
  const sl = Math.sin(latitude),
    cl = Math.cos(latitude);
  function project(ra, dec) {
    const ha = ra * rad - longitude,
      s = Math.sin(dec * rad),
      c = Math.cos(dec * rad);
    return {
      x: x - r * c * Math.sin(ha),
      y: y - r * (cl * s - sl * c * Math.cos(ha)),
      front: sl * s + cl * c * Math.cos(ha) >= 0,
    };
  }
  ctx.save();
  glow(ctx, x, y, r * 1.32, "106,168,176", 0.18);
  const face = ctx.createRadialGradient(x - 0.28 * r, y - 0.44 * r, 0, x, y, r);
  face.addColorStop(0, "#253c47");
  face.addColorStop(0.67, "#14242f");
  face.addColorStop(0.94, "#0b1621");
  face.addColorStop(1, "#2f4852");
  ctx.fillStyle = face;
  ctx.beginPath();
  ctx.arc(x, y, r, 0, TAU);
  ctx.fill();
  ctx.lineWidth = 1;
  function gridLine(start, end, fn, color) {
    ctx.strokeStyle = color;
    ctx.beginPath();
    let drawing = false;
    for (let angle = start; angle <= end; angle += 3) {
      const p = fn(angle);
      if (p.front) {
        if (drawing) ctx.lineTo(p.x, p.y);
        else ctx.moveTo(p.x, p.y);
      }
      drawing = p.front;
    }
    ctx.stroke();
  }
  for (let ra = 0; ra < 360; ra += 15)
    gridLine(-90, 90, (dec) => project(ra, dec), ra % 60 === 0 ? "#617d85" : "#3d535e");
  for (let dec = -75; dec <= 75; dec += 15)
    gridLine(0, 360, (ra) => project(ra, dec), dec === 0 ? "#9db4ad" : "#495f69");
  const point = project(coords.raDeg, coords.decDeg);
  if (r > 105 && approach > 0.8) {
    ctx.font = '10px "IBM Plex Mono", monospace';
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    ctx.fillStyle = "#a7babd";
    for (let ra = 0; ra < 360; ra += 45) {
      const pos = project(ra, 0);
      if (
        pos.front &&
        Math.cos(ra * rad - longitude) > 0.48 &&
        Math.hypot(pos.x - point.x, pos.y - point.y) > 48
      )
        ctx.fillText(ra / 15 + "h", pos.x, pos.y + 5);
    }
    const pole = project(0, latitude >= 0 ? 90 : -90);
    if (pole.front && Math.abs(latitude) > 0.3) {
      ctx.textBaseline = latitude >= 0 ? "bottom" : "top";
      ctx.fillText(latitude >= 0 ? "N" : "S", pole.x, pole.y + (latitude >= 0 ? -4 : 4));
    }
  }
  if (point.front && approach > 0.55) {
    ctx.strokeStyle = "#f7c18e";
    ctx.fillStyle = "#ffdab0";
    ctx.lineWidth = 1.3;
    ctx.beginPath();
    ctx.arc(point.x, point.y, 7, 0, TAU);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(point.x, point.y, 2.7, 0, TAU);
    ctx.fill();
    if (r > 100) {
      const lineY = point.y + Math.min(30, r * 0.18),
        lineX = point.x - Math.min(52, r * 0.34);
      ctx.beginPath();
      ctx.moveTo(point.x - 5, point.y + 6);
      ctx.lineTo(point.x - 23, lineY);
      ctx.lineTo(lineX, lineY);
      ctx.stroke();
      ctx.font = '10px "IBM Plex Mono", monospace';
      ctx.textAlign = "right";
      ctx.textBaseline = "middle";
      const textX = lineX - 4;
      for (const [label, dy] of [
        ["RA  " + coords.raDeg.toFixed(4) + "°", -7],
        ["Dec " + coords.decDeg.toFixed(4) + "°", 8],
      ]) {
        const width = ctx.measureText(label).width;
        ctx.fillStyle = "#0b111bec";
        ctx.fillRect(textX - width - 3, lineY + dy - 7, width + 6, 14);
        ctx.fillStyle = "#ffe0ba";
        ctx.fillText(label, textX, lineY + dy);
      }
    }
  }
  ctx.restore();
}

function buildScene(name, media, inspectable) {
  const root = element("div", "research-scene");
  const space = element("div", "rs-camera");
  space.tabIndex = 0;
  space.setAttribute("role", "group");
  space.setAttribute(
    "aria-label",
    `View of ${name}. Drag or use arrow keys to look around. Use plus and minus to zoom and Home to center.`,
  );
  const sky = element("canvas", "rs-stars");
  sky.setAttribute("aria-hidden", "true");
  space.append(sky);
  const plateNode = element("div", "rs-plate");
  plateNode.hidden = true;
  plateNode.dataset.mode = media?.mode || "station";
  const img = element("img", "rs-source");
  img.alt = media?.alt || "";
  img.draggable = false;
  plateNode.append(img);
  for (const marker of media?.markers || []) {
    const mark = element(inspectable ? "button" : "span", "rs-target-mark");
    if (inspectable) {
      mark.type = "button";
      mark.dataset.marker = marker.label;
      mark.setAttribute("aria-label", "Inspect lensed image " + marker.label);
      mark.title = "Inspect image " + marker.label;
    }
    mark.style.left = marker.x * 100 + "%";
    mark.style.top = marker.y * 100 + "%";
    mark.append(element("span", "", marker.label));
    plateNode.append(mark);
  }
  const scale = element("span", "rs-image-scale", "5″");
  scale.hidden = !media?.fieldArcsec;
  plateNode.append(scale);
  space.append(plateNode);
  const controls = element("div", "rs-controls");
  controls.setAttribute("role", "group");
  controls.setAttribute("aria-label", "Research camera controls");
  function button(label, content) {
    const b = element("button", "");
    b.type = "button";
    b.title = label;
    b.setAttribute("aria-label", label);
    b.innerHTML = content;
    controls.append(b);
    return b;
  }
  const smaller = button("Zoom out of research view", "−"),
    larger = button("Zoom into research view", "+");
  const reset = button("Center research camera", icon("focus"));
  const note = element("span", "rs-zoom-note", "1×");
  controls.append(note);
  const footer = element("div", "rs-credit");
  const credit = element(
    media?.href ? "a" : "span",
    "rs-origin",
    media ? "Loading source image…" : "Illustrative observatory aperture",
  );
  if (media?.href) {
    credit.href = media.href;
    credit.target = "_blank";
    credit.rel = "noopener";
  }
  footer.append(credit);
  const skip = element("button", "rs-skip", "Skip flight");
  skip.type = "button";
  const actions = element("div", "rs-actions");
  const mapToggle = element("button", "");
  mapToggle.type = "button";
  mapToggle.hidden = true;
  const inspect = element(
    "button",
    "",
    media?.mode === "illustration"
      ? "Inspect radio data"
      : "Inspect " + (media?.mode === "figure" ? "figure" : "image"),
  );
  if (media?.mode === "illustration") inspect.title = "Inspect the VLASS radio observations";
  inspect.type = "button";
  inspect.hidden = !inspectable;
  actions.append(mapToggle, inspect);
  root.append(space, controls, footer, skip, actions);
  return {
    root,
    space,
    sky,
    plateNode,
    img,
    scale,
    smaller,
    larger,
    reset,
    note,
    credit,
    skip,
    inspect,
    mapToggle,
  };
}

/** Mount one guided research flight. All camera offsets use display units. */
export function mountResearchScene(host, { destination, plate, animate = true, onArrive, onInspect } = {}) {
  const media = mediaFor(destination, plate);
  const ui = buildScene(destination.name, media, Boolean(media && onInspect));
  host.append(ui.root);
  const ctx = ui.sky.getContext("2d");
  const random = randomSource(destination.id);
  const tint = palettes[destination.id] || palettes[destination.kind] || [128, 157, 189];
  const dust = ctx && dustPlate(random, tint),
    stars = ctx ? starCatalog(random) : [];
  const events = new AbortController(),
    opts = { signal: events.signal };
  const motion = matchMedia("(prefers-reduced-motion: reduce)");
  const pointers = new Map();
  const pointing = destination.coordinates;
  const hasPointing =
    ctx && media?.mode === "figure" && Number.isFinite(pointing?.raDeg) && Number.isFinite(pointing?.decDeg);
  let mapView = Boolean(hasPointing);
  let tap = null;
  let w = 0,
    h = 0,
    dpr = 1,
    frame = 0,
    disposed = false,
    arrived = false,
    loaded = false,
    failed = false;
  let start = performance.now(),
    zoom = 1,
    panX = 0,
    panY = 0;
  let flying = Boolean(animate && !motion.matches && ctx && !document.hidden);
  ui.root.dataset.phase = flying ? "flight" : "arrival";
  ui.skip.hidden = !flying;
  if (!ctx) {
    ui.root.classList.add("rs-no-canvas");
  }
  function updateView() {
    ui.root.dataset.view = mapView ? "pointing" : failed ? "station" : media?.mode || "station";
    ui.mapToggle.hidden = !hasPointing || failed;
    ui.mapToggle.textContent = mapView ? "View figure" : "View pointing";
    if (mapView || !loaded) ui.plateNode.hidden = true;
    ui.inspect.hidden = !onInspect || !loaded;
    const available = loaded ? [media.kind, media.credit].filter(Boolean).join(" · ") : "Loading source image…";
    ui.credit.textContent = mapView
      ? "Equatorial pointing · " + pointing.frame
      : failed
        ? "Image unavailable · " + (media.credit || "Source record")
        : media
          ? available
          : "Illustrative observatory aperture";
    if (mapView) {
      ui.space.setAttribute(
        "aria-label",
        [
          pointing.label,
          "Right ascension " + pointing.raDeg + " degrees",
          "Declination " + pointing.decDeg + " degrees",
          "Drag or use arrow keys to rotate the celestial globe. Use plus or minus to zoom and Home to center",
        ].join(". "),
      );
      if (pointing.source?.url) ui.credit.href = pointing.source.url;
    } else {
      ui.space.setAttribute(
        "aria-label",
        `View of ${destination.name}. Drag or use arrow keys to look around. Use plus and minus to zoom and Home to center.`,
      );
      if (media?.href) ui.credit.href = media.href;
    }
  }
  updateView();
  function request() {
    if (!disposed && !frame) frame = requestAnimationFrame(draw);
  }
  function finish() {
    if (disposed || arrived) return;
    flying = false;
    arrived = true;
    if (document.activeElement === ui.skip) ui.space.focus({ preventScroll: true });
    ui.skip.hidden = true;
    ui.root.dataset.phase = "arrival";
    queueMicrotask(() => {
      if (!disposed) onArrive?.();
    });
  }
  function fit() {
    const rect = ui.space.getBoundingClientRect();
    if (!rect.width || !rect.height) return;
    w = rect.width;
    h = rect.height;
    dpr = Math.min(devicePixelRatio || 1, 2);
    ui.sky.width = Math.round(w * dpr);
    ui.sky.height = Math.round(h * dpr);
    request();
  }
  const observer = new ResizeObserver(fit);
  observer.observe(ui.space);
  fit();
  function cameraCenter() {
    return {
      x: w * (media?.mode === "figure" ? 0.5 : w < 500 ? 0.53 : 0.59),
      y: h * (w < 500 ? 0.54 : 0.535),
    };
  }
  function draw(now) {
    frame = 0;
    if (disposed) return;
    const t = flying ? clamp((now - start) / FLIGHT_MS, 0, 1) : 1;
    if (t >= 1) finish();
    if (!w || !h) {
      if (flying) request();
      return;
    }
    const travel = smooth(t),
      rush = Math.sin(Math.PI * t) ** 1.45;
    const { x: x0, y: y0 } = cameraCenter();
    const focusX = x0 - rush * w * 0.095,
      focusY = y0 + rush * h * 0.085;
    const approach = 0.021 + 0.979 * smooth(clamp((t - 0.1) / 0.9, 0, 1)) ** 1.65;
    const px = focusX + (mapView ? 0 : panX),
      py = focusY + (mapView ? 0 : panY);
    if (ctx) {
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      drawField(ctx, dust, stars, { w, h, focusX, focusY, panX, panY, zoom, t, speed: rush });
      glow(ctx, px, py, Math.max(28, Math.min(w * 0.51, h * 0.72) * zoom * approach), tint.join(","), 0.23);
      if (t < 0.74) glow(ctx, px, py, Math.max(7, 18 * rush), tint.join(","), 0.65 * rush);
      if (mapView)
        drawPointing(
          ctx,
          px,
          py,
          Math.min(w * 0.34, h * 0.32) * zoom * approach,
          pointing,
          approach,
          (panX * 420) / (w * zoom),
          (panY * 230) / (h * zoom),
        );
      else if (!media || failed)
        drawAperture(ctx, px, py, Math.min(w * 0.25, h * 0.31) * zoom * approach, tint);
    } else if (!media || failed) {
      ui.space.style.setProperty("--rs-x", px + "px");
      ui.space.style.setProperty("--rs-y", py + "px");
      ui.space.style.setProperty("--rs-diameter", 2 * Math.min(w * 0.25, h * 0.31) * zoom + "px");
    }
    if (loaded && !mapView) {
      const aspect = ui.img.naturalWidth / ui.img.naturalHeight;
      const room =
        media.mode === "illustration"
          ? [w * (w < 500 ? 1.44 : 1.2), h * 0.81]
          : media.mode === "figure"
            ? [w * 0.86, h * 0.62]
            : media.mode === "disc"
              ? [w * (w < 500 ? 0.8 : 0.68), h * 0.81]
              : [w * (w < 500 ? 0.81 : 0.7), h * 0.68];
      const iw = Math.min(room[0], room[1] * aspect) * zoom * approach,
        ih = iw / aspect;
      ui.plateNode.style.width = iw + "px";
      ui.plateNode.style.height = ih + "px";
      ui.plateNode.style.transform = `translate3d(${px - iw / 2}px,${py - ih / 2}px,0)`;
      ui.plateNode.style.opacity = String(smooth(clamp((travel - 0.05) / 0.45, 0, 1)));
      ui.plateNode.hidden = false;
      if (media.fieldArcsec) ui.scale.style.width = (iw * 5) / media.fieldArcsec + "px";
    }
    if (flying) request();
  }
  if (media) {
    ui.img.src = media.src;
    ui.img.decode().then(
      () => {
        if (disposed) return;
        loaded = true;
        updateView();
        request();
      },
      () => {
        if (disposed) return;
        failed = true;
        updateView();
        request();
      },
    );
  }
  function setCamera(z, dx, dy) {
    zoom = clamp(z, 0.65, 3.5);
    panX = clamp(dx, -w * 0.36 * zoom, w * 0.36 * zoom);
    panY = clamp(dy, -h * 0.33 * zoom, h * 0.33 * zoom);
    ui.note.textContent = Number(zoom.toFixed(1)) + "×";
    ui.smaller.disabled = zoom <= 0.65;
    ui.larger.disabled = zoom >= 3.5;
    request();
  }
  function magnify(factor, ax, ay) {
    const center = cameraCenter();
    ax ??= center.x;
    ay ??= center.y;
    const z = clamp(zoom * factor, 0.65, 3.5),
      ratio = z / zoom;
    setCamera(z, panX * ratio + (ax - center.x) * (1 - ratio), panY * ratio + (ay - center.y) * (1 - ratio));
  }
  function control(fn) {
    return () => {
      finish();
      fn();
    };
  }
  function inspect(marker) {
    if (!loaded || !onInspect) return;
    finish();
    onInspect({ marker });
  }
  ui.mapToggle.addEventListener(
    "click",
    control(() => {
      mapView = !mapView;
      updateView();
      setCamera(1, 0, 0);
    }),
    opts,
  );
  ui.inspect.addEventListener("click", () => inspect(), opts);
  ui.plateNode.addEventListener(
    "click",
    (e) => {
      const mark = e.target.closest("[data-marker]");
      if (mark) inspect(mark.dataset.marker);
    },
    opts,
  );
  ui.smaller.addEventListener(
    "click",
    control(() => magnify(1 / 1.32)),
    opts,
  );
  ui.larger.addEventListener(
    "click",
    control(() => magnify(1.32)),
    opts,
  );
  ui.reset.addEventListener(
    "click",
    control(() => setCamera(1, 0, 0)),
    opts,
  );
  ui.skip.addEventListener(
    "click",
    () => {
      finish();
      request();
    },
    opts,
  );
  ui.space.addEventListener(
    "keydown",
    (e) => {
      const arrows = { ArrowLeft: [27, 0], ArrowRight: [-27, 0], ArrowUp: [0, 27], ArrowDown: [0, -27] };
      if (arrows[e.key]) {
        finish();
        setCamera(zoom, panX + arrows[e.key][0], panY + arrows[e.key][1]);
      } else if (e.key === "+" || e.key === "=") {
        finish();
        magnify(1.32);
      } else if (e.key === "-") {
        finish();
        magnify(1 / 1.32);
      } else if (e.key === "Home") {
        finish();
        setCamera(1, 0, 0);
      } else return;
      e.preventDefault();
      e.stopPropagation();
    },
    opts,
  );
  function local(e) {
    const r = ui.space.getBoundingClientRect();
    return { x: e.clientX - r.left, y: e.clientY - r.top };
  }
  function gesture() {
    const [a, b] = pointers.values();
    if (!a) return null;
    return b ? { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, d: Math.hypot(a.x - b.x, a.y - b.y) } : { ...a, d: 0 };
  }
  ui.space.addEventListener(
    "pointerdown",
    (e) => {
      if (e.button !== 0 || pointers.size >= 2 || e.target.closest("[data-marker]")) return;
      tap = !pointers.size && !flying ? local(e) : null;
      ui.space.focus({ preventScroll: true });
      ui.space.setPointerCapture(e.pointerId);
      pointers.set(e.pointerId, local(e));
      ui.space.classList.add("rs-dragging");
    },
    opts,
  );
  ui.space.addEventListener(
    "pointermove",
    (e) => {
      if (!pointers.has(e.pointerId)) return;
      const before = gesture();
      pointers.set(e.pointerId, local(e));
      const after = gesture();
      if (tap && Math.hypot(after.x - tap.x, after.y - tap.y) > 4) tap = null;
      if (!before || !after || (!after.d && Math.hypot(after.x - before.x, after.y - before.y) < 1)) return;
      finish();
      if (before.d > 1 && after.d > 1) magnify(after.d / before.d, before.x, before.y);
      setCamera(zoom, panX + after.x - before.x, panY + after.y - before.y);
    },
    opts,
  );
  ui.space.addEventListener(
    "pointerup",
    (e) => {
      if (tap && loaded && !mapView && !e.target.closest("[data-marker]") && onInspect) {
        const rect = ui.plateNode.getBoundingClientRect();
        const u = (2 * (e.clientX - rect.left)) / rect.width - 1;
        const v = (2 * (e.clientY - rect.top)) / rect.height - 1;
        if (media.mode === "disc" ? u * u + v * v <= 1 : Math.max(Math.abs(u), Math.abs(v)) <= 1) inspect();
      }
      tap = null;
    },
    opts,
  );
  ui.space.addEventListener(
    "lostpointercapture",
    (e) => {
      pointers.delete(e.pointerId);
      if (!pointers.size) ui.space.classList.remove("rs-dragging");
    },
    opts,
  );
  ui.space.addEventListener(
    "wheel",
    (e) => {
      if (document.activeElement !== ui.space && !e.ctrlKey) return;
      e.preventDefault();
      e.stopPropagation();
      finish();
      const point = local(e);
      magnify(Math.exp(-clamp(e.deltaY, -120, 120) * 0.003), point.x, point.y);
    },
    { ...opts, passive: false },
  );
  document.addEventListener(
    "visibilitychange",
    () => {
      if (document.hidden) {
        finish();
        request();
      }
    },
    opts,
  );
  motion.addEventListener(
    "change",
    () => {
      if (motion.matches) {
        finish();
        request();
      }
    },
    opts,
  );
  if (!flying) finish();
  request();
  return () => {
    disposed = true;
    events.abort();
    observer.disconnect();
    cancelAnimationFrame(frame);
    pointers.clear();
    ui.root.remove();
  };
}
