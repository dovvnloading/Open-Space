import { esc } from "./html.js";
import { icon } from "./icons.js";

const levels = [1, 2, 4];
const titles = {
  mars: "Valles Marineris hemisphere",
  sun: "Seven hours in the solar corona",
  quasar: "Two images of one quasar",
  dwarf: "A brown dwarf, imagined",
};
const explanations = {
  sun: "Follow the outer corona through eight hourly STEREO-A exposures on 13 July 2021. Two eruptions overlap in this field. The source-direction test uses Solar Orbiter’s separate X-ray measurements.",
  quasar:
    "A and B mark the catalog positions of the two images in this 40-arcsecond Pan-STARRS1 field. Their separation is 5.83 arcseconds. The timing investigation uses ZTF and Gaia.",
};

function utc(iso) {
  return iso.slice(11, 19) + " UTC";
}

export function openImageDesk(id, plates) {
  const desk = document.createElement("image-desk");
  desk.targetId = id;
  desk.plates = plates;
  return desk;
}

class ImageDesk extends HTMLElement {
  connectedCallback() {
    const id = this.targetId;
    const plate = this.plates[id];
    this.frames = plate.sequence || [plate];
    this.index = 0;
    this.level = 0;
    this.offset = { x: 0, y: 0 };
    const sequence = this.frames.length > 1;
    this.innerHTML = `<article class="desk-layout">
      <header class="desk-head"><div><p class="desk-kicker">${esc(plate.kind)} · ${esc(plate.name)}</p>
      <h2 id="dossier-title">${esc(titles[id])}</h2></div></header>
      <div class="desk-window" role="group" tabindex="0" aria-label="Image inspection area. At higher magnification use the arrow keys to pan.">
        <div class="desk-sheet"><img draggable="false" alt="${esc(plate.kind + " of " + plate.name)}"/>
        ${(plate.markers || []).map((m) => `<span class="desk-pin" style="left:${100 * m.x}%;top:${100 * m.y}%"><b>${esc(m.label)}</b></span>`).join("")}</div>
        <span class="desk-scale" hidden>5″</span>
      </div>
      <div class="desk-tools"><div class="desk-zoom" role="group" aria-label="Image magnification">
        ${levels.map((n, i) => `<button type="button" data-zoom="${i}" aria-pressed="${i === 0}" aria-label="Magnify ${n} times">${n}×</button>`).join("")}
      </div><small data-hint>Fit to window</small>
      <a class="desk-original" target="_blank" rel="noopener">Original image ${icon("external")}</a></div>
      ${
        sequence
          ? `<div class="desk-film"><button type="button" data-play>Play sequence</button>
        <label><span>13 JUL 2021</span><output data-clock></output>
        <input type="range" min="0" max="${this.frames.length - 1}" step="1" value="${this.index}" aria-label="Solar frame time"></label></div>`
          : ""
      }
      <div class="desk-caption"><p>${esc(explanations[id] || plate.description)}</p>
        <a href="${esc(plate.source)}" target="_blank" rel="noopener">${esc(plate.credit)} ${icon("external")}</a>
        <a href="#target/${esc(id)}">Read investigation ${icon("arrow")}</a></div>
      <nav class="desk-plates" aria-label="Image collection">${Object.entries(this.plates)
        .map(
          ([key, p]) =>
            `<a href="#plate/${esc(key)}" ${key === id ? 'aria-current="page"' : ""}><img src="./assets/${esc(p.file)}" alt=""/><span>${esc(p.name)}</span></a>`,
        )
        .join("")}</nav>
      </article>`;
    const find = (selector) => this.querySelector(selector);
    this.stage = find(".desk-window");
    this.sheet = find(".desk-sheet");
    this.photo = find(".desk-sheet img");
    this.clock = find("[data-clock]");
    this.range = find("input[type=range]");
    this.playButton = find("[data-play]");
    this.observer = new ResizeObserver(() => this.layout());
    this.observer.observe(this.stage);
    this.photo.onload = () => this.layout();
    this.querySelectorAll("[data-zoom]").forEach(
      (button) =>
        (button.onclick = () => {
          this.level = Number(button.dataset.zoom);
          this.querySelectorAll("[data-zoom]").forEach((b) =>
            b.setAttribute("aria-pressed", String(b === button)),
          );
          if (!this.level) this.offset = { x: 0, y: 0 };
          find("[data-hint]").textContent = this.level ? "Drag or use arrow keys to pan" : "Fit to window";
          this.layout();
        }),
    );
    this.stage.onkeydown = (e) => {
      const delta = { ArrowLeft: [40, 0], ArrowRight: [-40, 0], ArrowUp: [0, 40], ArrowDown: [0, -40] }[
        e.key
      ];
      if (!delta || !this.level) return;
      e.preventDefault();
      this.pan(...delta);
    };
    this.stage.onpointerdown = (e) => {
      if (!this.level || !e.isPrimary || e.button !== 0) return;
      this.drag = { x: e.clientX, y: e.clientY };
      this.stage.setPointerCapture(e.pointerId);
    };
    this.stage.onpointermove = (e) => {
      if (!this.drag || !e.isPrimary) return;
      this.pan(e.clientX - this.drag.x, e.clientY - this.drag.y);
      this.drag = { x: e.clientX, y: e.clientY };
    };
    this.stage.onlostpointercapture = () => (this.drag = null);
    if (sequence) {
      this.range.oninput = () => {
        this.stop();
        this.showFrame(Number(this.range.value));
      };
      this.playButton.onclick = () => {
        if (this.timer) {
          this.stop();
          return;
        }
        if (this.index === this.frames.length - 1) this.showFrame(0);
        this.playButton.textContent = "Pause sequence";
        this.timer = setInterval(() => {
          if (this.index === this.frames.length - 1) this.stop();
          else this.showFrame(this.index + 1);
        }, 760);
      };
      this.visibility = () => {
        if (document.hidden) this.stop();
      };
      document.addEventListener("visibilitychange", this.visibility);
    }
    this.showFrame(this.index);
  }

  showFrame(index) {
    this.index = index;
    const frame = this.frames[index];
    this.photo.src = "./assets/" + frame.file;
    this.querySelector(".desk-original").href = frame.archiveUrl || "./assets/" + frame.file;
    if (this.range) {
      const time = utc(frame.observedAt);
      this.range.value = String(index);
      this.range.setAttribute("aria-valuetext", time);
      this.clock.textContent = time;
      this.photo.alt = "STEREO-A COR2 at " + time + " on 13 July 2021";
      const next = this.frames[index + 1];
      if (next) {
        const preload = new Image();
        preload.src = "./assets/" + next.file;
      }
    }
    this.layout();
  }

  pan(dx, dy) {
    this.offset.x += dx;
    this.offset.y += dy;
    this.layout();
  }

  layout() {
    const w = this.stage.clientWidth,
      h = this.stage.clientHeight;
    const iw = this.photo.naturalWidth,
      ih = this.photo.naturalHeight;
    if (!iw || !w || !h) return;
    const fit = Math.min(w / iw, h / ih) * levels[this.level];
    const width = iw * fit,
      height = ih * fit;
    const clamp = (v, bound) => Math.min(bound, Math.max(-bound, v));
    this.offset.x = clamp(this.offset.x, Math.max(0, (width - w) / 2));
    this.offset.y = clamp(this.offset.y, Math.max(0, (height - h) / 2));
    Object.assign(this.sheet.style, {
      width: width + "px",
      height: height + "px",
      left: (w - width) / 2 + this.offset.x + "px",
      top: (h - height) / 2 + this.offset.y + "px",
    });
    this.stage.classList.toggle("pannable", this.level > 0);
    const scale = this.querySelector(".desk-scale");
    const field = this.plates[this.targetId].fieldArcsec;
    scale.hidden = !field;
    if (field) scale.style.width = (width * 5) / field + "px";
  }

  stop() {
    clearInterval(this.timer);
    this.timer = null;
    if (this.playButton) this.playButton.textContent = "Play sequence";
  }

  disconnectedCallback() {
    this.stop();
    this.observer?.disconnect();
    if (this.visibility) document.removeEventListener("visibilitychange", this.visibility);
  }
}

customElements.define("image-desk", ImageDesk);
