import { axesFor, observationSeries } from "./observation-data.js";
import { icon } from "./icons.js";
import { esc } from "./html.js";

const axisNumber = new Intl.NumberFormat("en-US", {notation: "compact", maximumFractionDigits: 2});
const asCoord = n => Number(n.toFixed(2));

function glyph(x, y, radius, series, attributes = "") {
  if (series === 1) {
    return `<path d="M${x},${y - radius}l${radius},${radius} -${radius},${radius} -${radius},-${radius}Z" ${attributes}/>`;
  }
  return `<circle cx="${x}" cy="${y}" r="${radius}" ${attributes}/>`;
}

function frame(data, width) {
  const narrow = width < 465;
  const height = narrow ? 194 : 218;
  const left = narrow ? 36 : 46;
  const right = width - (narrow ? 19 : 29);
  const top = 13;
  const bottom = height - 32;
  const {x0, x1, lo, hi, ticksX, ticksY} = axesFor(data, narrow);
  const x = value => asCoord(left + (value - x0) / (x1 - x0 || 1) * (right - left));
  const project = data.log ? Math.log10 : value => value;
  const y = value => {
    const fraction = (project(value) - project(lo)) / (project(hi) - project(lo) || 1);
    return asCoord(data.inverted ? top + fraction * (bottom - top) : bottom - fraction * (bottom - top));
  };
  const eventBands = data.rows.filter(row => row.eventDay).map(row => {
    const start = Math.max(left, x(row.x - .5));
    const end = Math.min(right, x(row.x + .5));
    return `<rect x="${start}" y="${top}" width="${asCoord(end-start)}" height="${bottom-top}" class="obs-event-band"/>`;
  }).join("");
  const xGrid = ticksX.map(([value, label]) => {
    const px = x(value);
    return `<g class="obs-tick"><path class="obs-grid-v" d="M${px},${top}V${bottom}"/>
      <text x="${px}" y="${bottom + 19}" text-anchor="middle">${esc(label)}</text></g>`;
  }).join("");
  const yGrid = ticksY.map(value => {
    const py = y(value);
    return `<g class="obs-tick"><path class="${value === 0 ? "obs-zero" : "obs-grid-h"}" d="M${left},${py}H${right}"/>
      <text x="${left - 8}" y="${py + 3.6}" text-anchor="end">${esc(axisNumber.format(value))}</text></g>`;
  }).join("");
  const coords = data.rows.map((row, index) => row.points.map(p => ({
    x: x(row.x), y: y(p.value), index, series: p.series, spread: p.spread,
    low: y(p.value - p.spread), high: y(p.value + p.spread),
  })));
  let trace = "";
  if (data.connect) {
    const path = coords.map(([p], i) => `${i ? "L" : "M"}${p.x},${p.y}`).join(" ");
    trace = `<path class="obs-fill" d="${path}V${bottom}H${coords[0][0].x}Z"/>
      <path class="obs-line" d="${path}"/>`;
  }
  const marks = coords.flat().map(p => {
    const bar = p.spread ? `<path class="obs-whisker" d="M${p.x},${p.low}V${p.high}
      M${p.x - 3},${p.low}h6M${p.x - 3},${p.high}h6"/>` : "";
    return `<g class="obs-series-${p.series}">${bar}${glyph(p.x, p.y, data.id === "dwarf" ? 5 : 2.8, p.series,
      `class="obs-mark${data.rows[p.index].noCentroid ? " obs-hollow" : ""}"`)}</g>`;
  }).join("");
  const html = `<defs><linearGradient id="obs-curve-fill" x1="0" x2="0" y1="0" y2="1" gradientUnits="objectBoundingBox">
      <stop offset="0" stop-color="var(--signal-a)" stop-opacity=".31"/><stop offset="1" stop-color="var(--signal-a)" stop-opacity=".025"/>
    </linearGradient></defs>${eventBands}${xGrid}${yGrid}${trace}${marks}<g class="obs-cursor"></g>`;
  return {html, width, height, coords, top, bottom};
}

class ObservationViewer extends HTMLElement {
  connectedCallback() {
    queueMicrotask(() => {
      if (!this.isConnected) return;
      if (!this.data) this.initialize();
      if (this.data) this.observer.observe(this.querySelector('.obs-trace'));
    });
  }
  disconnectedCallback() {
    this.observer?.disconnect();
    cancelAnimationFrame(this.pendingDraw);
  }
  initialize() {
    const table = this.closest("article").querySelector(".table-reader table");
    if (!table) return;
    this.raw = Array.from(table.tBodies[0].rows, row => Array.from(row.cells, cell => cell.textContent));
    this.data = observationSeries(this.dataset.target, this.raw);
    if (!this.data) return;
    this.selected = this.data.initialRow;
    const {title, id, series, rows} = this.data;
    const legend = series.map((name, i) =>
      `<span><i class="obs-key obs-series-${i}" aria-hidden="true"></i>${esc(name)}</span>`).join("");
    this.insertAdjacentHTML("afterbegin", `<div class="obs-head"><h3>${esc(title)}</h3>
      <div class="obs-views" role="group" aria-label="Figure view">
        <button type="button" data-view="interactive" aria-pressed="true">Interactive plot</button>
        <button type="button" data-view="published" aria-pressed="false">Published figure</button></div></div>
      <div class="obs-lab"><div class="obs-chart-head"><div class="obs-legend">${legend}</div>
        ${id === "quasar" ? `<label class="obs-shift"><input type="checkbox"> Apply published offsets</label>` : ""}</div>
        <div class="obs-trace-body"><div class="obs-trace"><div class="obs-y-label"></div>
        <svg class="obs-plot" role="img" aria-label="${esc(title)}" aria-describedby="observation-chart-note"></svg>
        <div class="obs-x-label"></div></div>
        ${id === "dwarf" && rows.every(r => r.cutout !== undefined) ? `<aside class="obs-sky-panel">
          <div class="obs-epoch-caption"></div>
          <div class="obs-cutout-box" data-state="loading" aria-busy="true">
            <svg class="obs-cutout" role="img" aria-label="Selected VLASS survey cutout, north up">
              <image width="1782" height="561"/></svg>
            <div class="obs-cutout-status" role="status"><span>Loading cutout…</span>
              <a hidden href="${esc(this.querySelector(".science-plate figcaption a").href)}" target="_blank" rel="noopener">Source figure</a>
            </div></div>
          <span class="obs-cutout-key"><i aria-hidden="true"></i> Predicted position</span>
          <small>16″ field · north up</small></aside>` : ""}
        </div>
        <div class="obs-inspector"><div class="obs-readout" aria-live="off"></div>
          <div class="obs-scrub">
          <button type="button" data-move="-1" aria-label="Previous observation">${icon("left")}</button>
          <input type="range" min="1" max="${rows.length}" value="${this.selected + 1}" step="1"
            aria-label="Inspect observation">
          <button type="button" data-move="1" aria-label="Next observation">${icon("right")}</button>
          <span class="obs-count" aria-hidden="true"></span></div></div>
        <p class="obs-note" id="observation-chart-note"></p>
      </div>`);
    this.svg = this.querySelector(".obs-plot");
    this.slider = this.querySelector('input[type="range"]');
    this.plate = this.querySelector(".science-plate");
    this.plate.hidden = true;
    this.loadCutout();
    this.classList.add("obs-ready");
    this.addEventListener("click", event => {
      const view = event.target.closest("[data-view]")?.dataset.view;
      if (view) this.switchView(view);
      const move = event.target.closest("[data-move]")?.dataset.move;
      if (move) this.select(this.selected + Number(move), true);
    });
    this.slider.addEventListener("input", () => this.select(Number(this.slider.value) - 1));
    this.querySelector(".obs-shift input")?.addEventListener("change", event => {
      this.aligned = event.target.checked;
      this.data = observationSeries(id, this.raw, this.aligned);
      this.draw();
    });
    this.svg.addEventListener("pointerdown", event => this.inspectAt(event));
    this.svg.addEventListener("pointermove", event => {
      if (event.pointerType === "mouse" || event.buttons) this.inspectAt(event);
    });
    this.observer = new ResizeObserver(() => {
      if (!this.svg.clientWidth || this.width === this.svg.clientWidth) return;
      cancelAnimationFrame(this.pendingDraw);
      this.pendingDraw = requestAnimationFrame(() => this.draw());
    });
    this.draw();
  }
  loadCutout() {
    const box = this.querySelector(".obs-cutout-box");
    if (!box) return;
    const raster = box.querySelector("image");
    const status = box.querySelector(".obs-cutout-status");
    const fail = () => {
      box.dataset.state = "error";
      box.setAttribute("aria-busy", "false");
      status.querySelector("span").textContent = "Cutout unavailable";
      status.querySelector("a").hidden = false;
    };
    raster.addEventListener("error", fail, {once: true});
    raster.addEventListener("load", async () => {
      try {
        if (raster.decode) await raster.decode();
        box.dataset.state = "ready";
        box.setAttribute("aria-busy", "false");
        status.hidden = true;
      } catch {
        fail();
      }
    }, {once: true});
    // Register listeners first so cached images follow the same ready state.
    raster.setAttribute("href", this.plate.querySelector("img").src);
  }
  switchView(view) {
    this.querySelector(".obs-lab").hidden = view !== "interactive";
    this.plate.hidden = view !== "published";
    this.querySelectorAll("[data-view]").forEach(button => {
      button.setAttribute("aria-pressed", String(button.dataset.view === view));
    });
    if (!this.plate.hidden) return;
    this.draw();
  }
  draw() {
    const width = this.svg.clientWidth;
    if (!width) return;
    this.width = width;
    this.frame = frame(this.data, width);
    this.svg.setAttribute("viewBox", `0 0 ${width} ${this.frame.height}`);
    this.svg.innerHTML = this.frame.html;
    this.querySelector(".obs-y-label").textContent = this.data.yLabel;
    this.querySelector(".obs-x-label").textContent = this.data.xLabel + (this.aligned ? " · B shifted" : "");
    this.querySelector(".obs-note").textContent = this.aligned ? this.data.alignedNote : this.data.note;
    this.select(this.selected);
  }
  inspectAt(event) {
    if (!this.frame) return;
    const rect = this.svg.getBoundingClientRect();
    const px = (event.clientX - rect.left) * this.frame.width / rect.width;
    const py = (event.clientY - rect.top) * this.frame.height / rect.height;
    let best = Infinity;
    let nearest = this.selected;
    for (const coord of this.frame.coords.flat()) {
      // Dense, nearly simultaneous Gaia pairs must be selectable by image.
      const distance = (coord.x - px) ** 2 + (this.data.id === "quasar" ? (coord.y - py) ** 2 : 0);
      if (distance < best) {
        best = distance;
        nearest = coord.index;
      }
    }
    if (nearest !== this.selected) this.select(nearest);
  }
  select(next, announce = false) {
    if (!this.frame) return;
    const {rows} = this.data;
    this.selected = Math.max(0, Math.min(rows.length - 1, next));
    const row = rows[this.selected];
    const marks = this.frame.coords[this.selected];
    const {top, bottom} = this.frame;
    const cursor = this.svg.querySelector(".obs-cursor");
    cursor.innerHTML = `<path class="obs-cross" d="M${marks[0].x},${top}V${bottom}"/>` +
      marks.map(p => `<g class="obs-series-${p.series}"><circle class="obs-halo" cx="${p.x}" cy="${p.y}" r="10"/>
        ${glyph(p.x, p.y, 4.4, p.series, `class="obs-selected${row.noCentroid ? " obs-hollow" : ""}"`)}</g>`).join("");
    const readout = this.querySelector(".obs-readout");
    readout.setAttribute("aria-live", announce ? "polite" : "off");
    const fields = row.fields.map(([name, value]) =>
      `<div><dt>${esc(name)}</dt><dd>${esc(value)}</dd></div>`).join("");
    readout.innerHTML = `<div class="obs-stamp">${esc(row.stamp)}</div><dl>${fields}</dl>`;
    const cutout = this.querySelector(".obs-cutout");
    if (cutout) {
      cutout.setAttribute("viewBox", `${row.cutout} 108 353 353`);
      this.querySelector(".obs-epoch-caption").textContent = row.stamp;
    }
    this.slider.value = String(this.selected + 1);
    this.slider.setAttribute("aria-valuetext", `${row.stamp}. ${row.fields.map(pair => pair.join(" ")).join(", ")}.`);
    this.slider.style.setProperty("--obs-progress", `${100 * this.selected / Math.max(1, rows.length - 1)}%`);
    this.querySelector(".obs-count").textContent = `${this.selected + 1} / ${rows.length}`;
    this.querySelector('[data-move="-1"]').setAttribute("aria-disabled", String(this.selected === 0));
    this.querySelector('[data-move="1"]').setAttribute("aria-disabled", String(this.selected === rows.length - 1));
  }
}

customElements.define("observation-viewer", ObservationViewer);
