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

export function openImageDesk(id, plates, { embedded = false } = {}) {
  const desk = document.createElement("image-desk");
  desk.targetId = id;
  desk.plates = plates;
  desk.embedded = embedded;
  return desk;
}

class ImageDesk extends HTMLElement {
  connectedCallback() {
    const id = this.targetId;
    const plate = this.plates[id];
    this.frames = plate.sequence || [plate];
    this.index = null;
    this.requestId ??= 0;
    this.playing = false;
    this.level = 0;
    this.offset = { x: 0, y: 0 };
    this.centerMarker = null;
    const sequence = this.frames.length > 1;
    this.innerHTML = `<article class="desk-layout">
      <header class="desk-head"><div><p class="desk-kicker">${esc(plate.kind)} · ${esc(plate.name)}</p>
      ${this.embedded ? `<h3>${esc(titles[id])}</h3>` : `<h2 id="dossier-title">${esc(titles[id])}</h2>`}</div></header>
      <div class="desk-window" role="group" tabindex="0" aria-label="Image inspection area. At higher magnification use the arrow keys to pan.">
        <div class="desk-sheet" hidden><img draggable="false" alt="${esc(plate.kind + " of " + plate.name)}"/>
        ${(plate.markers || []).map((m) => `<span class="desk-pin" data-pin="${esc(m.label)}" style="left:${100 * m.x}%;top:${100 * m.y}%"><b>${esc(m.label)}</b></span>`).join("")}</div>
        <span class="desk-scale" hidden>5″</span>
        <div class="desk-load-state" hidden><span role="status"></span><button type="button" hidden>Retry image</button></div>
      </div>
      <div class="desk-tools"><div class="desk-zoom" role="group" aria-label="Image magnification">
        ${levels.map((n, i) => `<button type="button" data-zoom="${i}" aria-pressed="${i === 0}" aria-label="Magnify ${n} times">${n}×</button>`).join("")}
      </div><small data-hint>Fit to window</small>
      <a class="desk-original" hidden target="_blank" rel="noopener">Original image ${icon("external")}</a></div>
      ${this.embedded && plate.markers?.length ? `<div class="desk-marker-tools" role="group" aria-label="Catalog image positions">
        ${plate.markers.map(marker => `<button type="button" data-marker="${esc(marker.label)}" aria-pressed="false">Image ${esc(marker.label)}</button>`).join("")}
        <span role="status" data-marker-status></span></div>` : ""}
      ${
        sequence
          ? `<div class="desk-film"><button type="button" data-play disabled>Play sequence</button>
        <label><span>13 JUL 2021</span><output data-clock></output>
        <input type="range" min="0" max="${this.frames.length - 1}" step="1" value="0" aria-label="Solar frame time"></label></div>`
          : ""
      }
      <div class="desk-caption"><p>${esc(explanations[id] || plate.description)}</p>
        <a href="${esc(plate.source)}" target="_blank" rel="noopener">${esc(plate.credit)} ${icon("external")}</a>
        ${this.embedded ? "" : `<a href="#target/${esc(id)}">Read investigation ${icon("arrow")}</a>`}</div>
      ${this.embedded ? "" : `<nav class="desk-plates" aria-label="Image collection">${Object.entries(this.plates)
        .map(
          ([key, p]) =>
            `<a href="#plate/${esc(key)}" ${key === id ? 'aria-current="page"' : ""}><img src="./assets/${esc(p.file)}" alt=""/><span>${esc(p.name)}</span></a>`,
        )
        .join("")}</nav>`}
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
    this.notice = find(".desk-load-state");
    this.message = find(".desk-load-state span");
    this.retry = find(".desk-load-state button");
    this.original = find(".desk-original");
    this.retry.onclick = () => {
      this.stage.focus();
      this.showFrame(this.requestedIndex);
    };
    this.querySelectorAll("[data-zoom]").forEach(button => {
      button.onclick = () => this.zoom(Number(button.dataset.zoom));
    });
    this.querySelectorAll("[data-marker]").forEach(button => {
      button.onclick = () => this.focusMarker(button.dataset.marker);
    });
    this.stage.onkeydown = (e) => {
      const delta = { ArrowLeft: [40, 0], ArrowRight: [-40, 0], ArrowUp: [0, 40], ArrowDown: [0, -40] }[
        e.key
      ];
      if (!delta || !this.level) return;
      e.preventDefault();
      this.pan(...delta);
    };
    this.stage.onpointerdown = (e) => {
      if (!this.level || !e.isPrimary || e.button !== 0 || e.target.closest("button")) return;
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
        if (this.playing) {
          this.stop();
          if (this.pendingImage) this.showFrame(this.index);
          return;
        }
        this.playing = true;
        this.playButton.textContent = "Pause sequence";
        if (this.index === this.frames.length - 1) this.showFrame(0);
        else if (!this.pendingImage) this.scheduleNext();
      };
      this.visibility = () => {
        if (document.hidden) this.stop();
      };
      document.addEventListener("visibilitychange", this.visibility);
    }
    this.showFrame(0);
  }

  cancelPending() {
    this.requestId++;
    clearTimeout(this.loadingTimer);
    if (this.pendingImage) this.pendingImage.removeAttribute("src");
    this.pendingImage = null;
  }

  scheduleNext() {
    clearTimeout(this.timer);
    if (!this.playing) return;
    if (this.index === this.frames.length - 1) this.stop();
    else this.timer = setTimeout(() => this.showFrame(this.index + 1), 760);
  }

  async showFrame(index) {
    this.cancelPending();
    this.requestedIndex = index;
    const request = this.requestId;
    const frame = this.frames[index];
    const time = this.range ? utc(frame.observedAt) : "";
    this.retry.hidden = true;
    if (this.range) {
      this.range.value = String(index);
      this.range.setAttribute("aria-valuetext", time);
    }
    if (index === this.index) {
      this.sheet.removeAttribute("aria-busy");
      this.notice.hidden = true;
      this.scheduleNext();
      return;
    }
    this.sheet.setAttribute("aria-busy", "true");
    const pending = this.preload?.index === index ? this.preload.image : new Image();
    this.preload = null;
    this.pendingImage = pending;
    const file = "./assets/" + frame.file;
    if (!pending.hasAttribute("src")) pending.src = file;
    this.notice.hidden = true;
    this.loadingTimer = setTimeout(() => {
      this.message.textContent = time ? `Loading ${time}…` : "Loading image…";
      this.notice.hidden = false;
    }, 150);
    try {
      await pending.decode();
    } catch {
      if (request !== this.requestId) return;
      clearTimeout(this.loadingTimer);
      this.pendingImage = null;
      this.sheet.removeAttribute("aria-busy");
      this.stop();
      this.message.textContent = time ? `Could not load ${time}.` : "Could not load this image.";
      this.notice.hidden = false;
      this.retry.hidden = false;
      return;
    }
    if (request !== this.requestId || !this.isConnected) return;
    clearTimeout(this.loadingTimer);
    this.pendingImage = null;
    this.index = index;
    this.notice.hidden = true;
    pending.draggable = false;
    pending.alt = time ? `STEREO-A COR2 at ${time} on 13 July 2021` : this.photo.alt;
    this.photo.replaceWith(pending);
    this.photo = pending;
    this.sheet.hidden = false;
    this.sheet.removeAttribute("aria-busy");
    this.original.href = frame.archiveUrl || file;
    this.original.hidden = false;
    if (time) {
      this.clock.textContent = time;
      this.playButton.disabled = false;
      const next = this.frames[index + 1];
      if (next) {
        const image = new Image();
        image.src = "./assets/" + next.file;
        this.preload = { index: index + 1, image };
      }
    }
    this.layout();
    this.scheduleNext();
  }

  zoom(level) {
    this.level = level;
    if (!level) this.offset = { x: 0, y: 0 };
    this.querySelectorAll("[data-zoom]").forEach(button =>
      button.setAttribute("aria-pressed", String(Number(button.dataset.zoom) === level))
    );
    this.querySelector("[data-hint]").textContent = level ? "Drag or use arrow keys to pan" : "Fit to window";
    this.layout();
  }

  focusMarker(label) {
    const marker = this.plates[this.targetId]?.markers?.find(item => item.label === label);
    if (!marker) return;
    this.centerMarker = marker;
    this.querySelectorAll("[data-pin]").forEach(pin => pin.classList.toggle("desk-pin-selected", pin.dataset.pin === label));
    this.querySelectorAll("[data-marker]").forEach(button => button.setAttribute("aria-pressed", String(button.dataset.marker === label)));
    const status = this.querySelector("[data-marker-status]");
    if (status) status.textContent = `Image ${label} selected`;
    this.zoom(levels.length - 1);
  }

  pan(dx, dy) {
    this.centerMarker = null;
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
    if (this.centerMarker) {
      this.offset.x = (.5 - this.centerMarker.x) * width;
      this.offset.y = (.5 - this.centerMarker.y) * height;
    }
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
    clearTimeout(this.timer);
    this.playing = false;
    this.timer = null;
    if (this.playButton) this.playButton.textContent = "Play sequence";
  }

  disconnectedCallback() {
    this.cancelPending();
    this.preload?.image.removeAttribute("src");
    this.preload = null;
    this.stop();
    this.observer?.disconnect();
    if (this.visibility) document.removeEventListener("visibilitychange", this.visibility);
  }
}

customElements.define("image-desk", ImageDesk);
